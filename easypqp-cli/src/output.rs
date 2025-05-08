use std::fs::{OpenOptions, metadata};
use std::io::{Write, BufWriter};
use std::path::Path;

use easypqp_core::{InsilicoPQPSettings, PeptideProperties};
use redeem_properties::utils::peptdeep_utils::remove_mass_shift;
use sage_core::peptide::Peptide;


pub fn write_assays_to_tsv<P: AsRef<Path>>(
    assays: &[PeptideProperties],
    peptides: &[Peptide],
    path: P,
    insilico_settings: &InsilicoPQPSettings,
) -> std::io::Result<()> {
    let path = path.as_ref();
    let write_header = !path.exists() || metadata(path)?.len() == 0;

    let file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)?;
    let mut writer = BufWriter::new(file);

    if write_header {
        writeln!(
            writer,
            "PrecursorMz\t\
            ProductMz\t\
            PrecursorCharge\t\
            ProductCharge\t\
            LibraryIntensity\t\
            NormalizedRetentionTime\t\
            PeptideSequence\t\
            ModifiedPeptideSequence\t\
            PeptideGroupLabel\t\
            LabelType\t\
            CompoundName\t\
            SumFormula\t\
            SMILES\t\
            Adducts\t\
            ProteinId\t\
            UniprotId\t\
            GeneName\t\
            FragmentType\t\
            FragmentSeriesNumber\t\
            Annotation\t\
            CollisionEnergy\t\
            PrecursorIonMobility\t\
            TransitionGroupId\t\
            TransitionId\t\
            Decoy\t\
            DetectingTransition\t\
            IdentifyingTransition\t\
            QuantifyingTransition\t\
            Peptidoforms"
        )?;
    }

    for assay in assays {
        let peptide_idx = assay.peptide_index as usize;
        let decoy = peptides[peptide_idx].decoy;
        let modified_peptide = peptides[peptide_idx].to_string();
        let naked_peptide = remove_mass_shift(&modified_peptide);
        let protein = &peptides[peptide_idx].proteins;

        let non_zero_indices: Vec<usize> = assay.product.intensity
            .iter()
            .enumerate()
            .filter(|(_, &intensity)| intensity > 0.0)
            .map(|(i, _)| i)
            .collect();

        if non_zero_indices.len() < insilico_settings.min_transitions as usize {
            continue;
        }

        let product_indices = if insilico_settings.max_transitions > 0 {
            let mut sorted = non_zero_indices.clone();
            sorted.sort_by(|&a, &b| {
                assay.product.intensity[b]
                    .partial_cmp(&assay.product.intensity[a])
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            sorted.into_iter()
                .take(insilico_settings.max_transitions as usize)
                .collect()
        } else {
            non_zero_indices
        };

        for &i in &product_indices {
            let fragment_type = &assay.product.ion_type[i];
            let series_number = assay.product.ion_ordinal[i];
            let product_charge = assay.product.charge[i];
            let annotation = format!("{}{}^{}", fragment_type, series_number, product_charge);

            writeln!(
                writer,
                "{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                assay.precursor.precursor_mz,
                assay.product.product_mz[i],
                assay.precursor.charge,
                product_charge,
                assay.product.intensity[i] * 10_000.0,
                assay.retention_time,
                naked_peptide,
                modified_peptide,
                "", "", "", "", "", "", // Placeholder columns
                protein.join(";"),
                "", "", // UniprotId, GeneName
                fragment_type,
                series_number,
                annotation,
                "", // CollisionEnergy
                assay.precursor.ion_mobility,
                assay.peptide_index,
                assay.peptide_index * 1000 + i as u32,
                decoy as i8,
                1, 0, 1,
                "" // Peptidoforms
            )?;
        }
    }

    Ok(())
}