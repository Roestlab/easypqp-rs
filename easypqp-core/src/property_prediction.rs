use anyhow::Result;
use chrono::Utc;
use core::f32;
use itertools::multiunzip;
use rayon::iter::ParallelIterator;
use rayon::prelude::*;
use redeem_properties::models::model_interface::PredictionValue;
use redeem_properties::utils::logging::Progress;
use redeem_properties::{
    models::model_interface::DLModels,
    utils::peptdeep_utils::{
        ccs_to_mobility_bruker, get_modification_string, remove_mass_shift, ModificationMap,
    },
};
use rustyms::fragment::FragmentKind;
use rustyms::system::e;
use rustyms::system::usize::Charge;
use rustyms::FragmentationModel;
use sage_core::peptide::Peptide;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use crate::{
    is_allowed_fragment, select_model, InsilicoPQPSettings, PeptideProperties, PrecursorProperties,
    ProductProperties,
};


#[derive(Clone)]
struct PeptideFeatures {
    naked_sequence: String,
    mod_string: String,
    mod_sites: String,
}

pub struct PropertyPrediction<'db, 'params> {
    pub peptides: &'db [Peptide],
    pub insilico_settings: &'params InsilicoPQPSettings,
    pub dl_models: DLModels,
    pub modifications: HashMap<(String, Option<char>), ModificationMap>,
    pub batch_size: usize,
    pub predictions:
        Arc<Mutex<HashMap<(u32, u8), (Option<f32>, Option<f32>, Option<PredictionValue>)>>>,
}

impl<'db, 'params> PropertyPrediction<'db, 'params> {
    pub fn new(
        peptides: &'db [Peptide],
        insilico_settings: &'params InsilicoPQPSettings,
        dl_models: DLModels,
        modifications: HashMap<(String, Option<char>), ModificationMap>,
        batch_size: usize,
    ) -> Self
    {
        PropertyPrediction {
            peptides,
            insilico_settings,
            dl_models,
            modifications,
            batch_size,
            predictions: Arc::new(Mutex::new(HashMap::new())),
        }
    }

    pub fn predict_properties(&mut self) -> Result<Vec<PeptideProperties>> {
        let rt_model = self.dl_models.rt_model.as_ref().cloned();
        let ccs_model = self.dl_models.ccs_model.as_ref().cloned();
        let ms2_model = self.dl_models.ms2_model.as_ref().cloned();

        // Step 1: Extract features once per peptide
        let total_batches = self.peptides.len();
        let timestamp = Utc::now().format("[%Y-%m-%dT%H:%M:%SZ INFO  easypqp_core::property_prediction]").to_string();
        let description = format!("{} Preparing features...", timestamp);
        let progress = Progress::new(total_batches, &description);
        let (send, recv) = crossbeam_channel::bounded(1024);
        let progress_thread = std::thread::spawn(move || {
            while recv.recv().is_ok() {
                progress.inc();
            }
            progress.finish();
        });
        let peptide_features: Vec<PeptideFeatures> = self.peptides
            .par_iter()
            .map_init(
                || send.clone(),
                |sender, peptide| {
                let peptide_string = peptide.to_string();
                let _ = sender.send(());
                PeptideFeatures {
                    naked_sequence: remove_mass_shift(&peptide_string)
                        .trim_start_matches('-')
                        .to_string(),
                    mod_string: get_modification_string(&peptide_string, &self.modifications),
                    mod_sites: peptide
                        .modifications
                        .iter()
                        .enumerate()
                        .filter_map(|(i, &v)| (v > 0.0).then(|| i.to_string()))
                        .collect::<Vec<_>>()
                        .join(";"),
                }
            })
            .collect();
        drop(send);
        progress_thread.join().unwrap();

        // Step 2: Predict RT in batches
        let batched_inputs: Vec<(Vec<String>, Vec<String>, Vec<String>)> = peptide_features
            .chunks(self.batch_size)
            .map(|chunk| {
                (
                    chunk.iter().map(|f| f.naked_sequence.clone()).collect(),
                    chunk.iter().map(|f| f.mod_string.clone()).collect(),
                    chunk.iter().map(|f| f.mod_sites.clone()).collect(),
                )
            })
            .collect();

        let total_batches = batched_inputs.len();
        let timestamp = Utc::now().format("[%Y-%m-%dT%H:%M:%SZ INFO  easypqp_core::property_prediction]").to_string();
        let description = format!("{} Predicting RT...", timestamp);
        let progress = Progress::new(total_batches, &description);
        let (send, recv) = crossbeam_channel::bounded(1024);
        let progress_thread = std::thread::spawn(move || {
            while recv.recv().is_ok() {
                progress.inc();
            }
            progress.finish();
        });

        let rt_predictions: Arc<Vec<Option<f32>>> = Arc::new(
            batched_inputs
                .par_iter()
                .map_init(
                    || send.clone(),
                    |sender, (seqs, mods, sites)| {
                        let preds = rt_model
                            .as_ref()
                            .and_then(|m| m.predict(seqs, mods, sites).ok());
                        let results = (0..seqs.len())
                            .map(|i| preds.as_ref().map(|p| p.get_prediction_entry(i)[0]))
                            .collect::<Vec<_>>();
                        let _ = sender.send(());
                        results
                    },
                )
                .flatten()
                .collect(),
        );
        drop(send);
        progress_thread.join().unwrap();

        // Step 3: Expand peptides by charge
        let timestamp = Utc::now().format("[%Y-%m-%dT%H:%M:%SZ INFO  easypqp_core::property_prediction]").to_string();
        let description = format!("{} Expanding features...", timestamp);
        let progress = Progress::new(self.peptides.len(), &description);
        let (send, recv) = crossbeam_channel::bounded(1024);
        let progress_thread = std::thread::spawn(move || {
            while recv.recv().is_ok() {
                progress.inc();
            }
            progress.finish();
        });

        let expanded: Vec<_> = peptide_features
            .par_iter()
            .enumerate()
            .map_init(
                || send.clone(),
                |sender, (idx, f)| {
                    let _ = sender.send(());
                    self.insilico_settings
                        .precursor_charge
                        .iter()
                        .map(move |&charge| {
                            (
                                f.naked_sequence.clone(),
                                f.mod_string.clone(),
                                f.mod_sites.clone(),
                                charge as i32,
                                (idx as u32, charge as u8),
                            )
                        })
                        .collect::<Vec<_>>()
                },
            )
            .flatten()
            .collect();
        drop(send);
        progress_thread.join().unwrap();

        let (seqs, mods, sites, charges, peptide_indices): (
            Vec<String>,
            Vec<String>,
            Vec<String>,
            Vec<i32>,
            Vec<(u32, u8)>,
        ) = multiunzip(expanded);

        // Step 4: Predict CCS + MS2
        let param = self
            .dl_models
            .params
            .as_ref()
            .expect("DL model parameters must be present")
            .clone();

        let batches: Vec<_> = (0..seqs.len())
            .step_by(self.batch_size)
            .map(|start| {
                let end = (start + self.batch_size).min(seqs.len());
                (
                    seqs[start..end].to_vec(),
                    mods[start..end].to_vec(),
                    sites[start..end].to_vec(),
                    charges[start..end].to_vec(),
                    peptide_indices[start..end].to_vec(),
                )
            })
            .collect();

        let total_batches = batches.len();
        let timestamp = Utc::now().format("[%Y-%m-%dT%H:%M:%SZ INFO  easypqp_core::property_prediction]").to_string();
        let description = format!("{} Predicting CCS and MS2...", timestamp);
        let progress = Progress::new(total_batches, &description);
        let (send, recv) = crossbeam_channel::bounded(1024);
        let progress_thread = std::thread::spawn(move || {
            while recv.recv().is_ok() {
                progress.inc();
            }
            progress.finish();
        });

        let predictions: Vec<_> = batches
            .par_iter()
            .map_init(
                || (send.clone(), Arc::clone(&rt_predictions)),
                |(sender, rt_preds), (seqs, mods, sites, chgs, indices)| {
                    let ccs_preds = ccs_model
                        .as_ref()
                        .and_then(|m| m.predict(seqs, mods, sites, chgs.clone()).ok());
                    let ms2_preds = ms2_model.as_ref().and_then(|m| {
                        m.predict(
                            seqs,
                            mods,
                            sites,
                            chgs.clone(),
                            vec![param.nce as i32; seqs.len()],
                            vec![param.instrument.clone(); seqs.len()],
                        )
                        .ok()
                    });

                    let result: Vec<_> = indices
                        .iter()
                        .enumerate()
                        .map(|(i, (peptide_idx, charge))| {
                            (
                                (*peptide_idx, *charge),
                                (
                                    rt_preds.get(*peptide_idx as usize).copied().flatten(),
                                    ccs_preds.as_ref().map(|p| p.get_prediction_entry(i)[0]),
                                    ms2_preds.as_ref().map(|p| p.get_prediction_entry(i).clone()),
                                ),
                            )
                        })
                        .collect();
                    let _ = sender.send(());
                    result
                },
            )
            .flatten()
            .collect();

        drop(send);
        progress_thread.join().unwrap();

        {
            let mut final_predictions = self.predictions.lock().unwrap();
            for (key, value) in predictions {
                final_predictions.insert(key, value);
            }
        }

        // Step 4: Build assays from predictions
        let preds_clone = {
            let preds = self.predictions.lock().unwrap();
            preds.clone()
        };

        let total_assays = preds_clone.len();
        let timestamp = Utc::now().format("[%Y-%m-%dT%H:%M:%SZ INFO  easypqp_core::property_prediction]").to_string();
        let description = format!("{} Creating PQP assays...", timestamp);
        let progress = Progress::new(total_assays, &description);
        let (send, recv) = crossbeam_channel::bounded(1024);

        let progress_thread = std::thread::spawn(move || {
            while recv.recv().is_ok() {
                progress.inc();
            }
            progress.finish();
        });

        let model = select_model(
            &self.insilico_settings.fragmentation_model,
            FragmentationModel::cid_hcd(),
        );

        let assays: Vec<_> = preds_clone
            .par_iter()
            .map_init(
                || (model, send.clone()),
                |(model, sender), ((peptide_idx, charge), (rt_pred, ccs_pred, ms2_pred))| {
                    let peptide = &self.peptides[*peptide_idx as usize];
                    let precursor_mz = peptide.monoisotopic / (*charge as f32);

                    let peptidoform = rustyms::peptidoform::PeptidoformIon::pro_forma(
                        &peptide.to_string(),
                        None,
                    )?;

                    let fragments = peptidoform.generate_theoretical_fragments(
                        Charge::new::<e>(self.insilico_settings.max_fragment_charge),
                        *model,
                    );

                    let mut product = ProductProperties {
                        peptide_index: *peptide_idx,
                        ion_type: Vec::new(),
                        ion_ordinal: Vec::new(),
                        charge: Vec::new(),
                        product_mz: Vec::new(),
                        intensity: Vec::new(),
                    };

                    for fragment in fragments {
                        if !is_allowed_fragment(
                            fragment.ion.kind(),
                            &self.insilico_settings.allowed_fragment_types,
                        ) || !fragment.neutral_loss.is_empty()
                        {
                            continue;
                        }

                        let ion = &fragment.ion;
                        let mz = fragment.mz(rustyms::MassMode::Monoisotopic).unwrap().value;

                        let (row, col) = match (ion.kind(), fragment.charge.value) {
                            (FragmentKind::b, 1) => {
                                (ion.position().unwrap().series_number as u8 - 1, 0)
                            }
                            (FragmentKind::b, 2) => {
                                (ion.position().unwrap().series_number as u8 - 1, 1)
                            }
                            (FragmentKind::y, 1) => {
                                (ion.position().unwrap().series_number as u8 - 1, 2)
                            }
                            (FragmentKind::y, 2) => {
                                (ion.position().unwrap().series_number as u8 - 1, 3)
                            }
                            _ => continue,
                        };

                        let intensity = ms2_pred
                            .as_ref()
                            .and_then(|pred| pred.get(row as usize, col as usize))
                            .copied()
                            .unwrap_or(0.0);

                        product.ion_type.push(ion.kind());
                        product
                            .ion_ordinal
                            .push(ion.position().unwrap().series_number as u8);
                        product.charge.push(fragment.charge.value as u8);
                        product.product_mz.push(mz);
                        product.intensity.push(intensity);
                    }

                    let _ = sender.send(());

                    Ok(PeptideProperties {
                        peptide_index: *peptide_idx,
                        retention_time: rt_pred.unwrap_or_default(),
                        precursor: PrecursorProperties {
                            peptide_index: *peptide_idx,
                            charge: *charge,
                            precursor_mz,
                            ion_mobility: ccs_pred
                                .map(|ccs| {
                                    ccs_to_mobility_bruker(
                                        ccs as f64,
                                        *charge as f64,
                                        precursor_mz as f64,
                                    )
                                })
                                .unwrap_or(f64::NAN),
                        },
                        product,
                    })
                },
            )
            .collect::<Result<Vec<_>>>()?;

        drop(send);
        progress_thread.join().unwrap();

        Ok(assays)
    }
}
