//! UniMod database parser and mass-bracket-to-UniMod re-annotation.
//!
//! Provides [`UnimodDb`] for parsing the [UniMod XML](https://www.unimod.org) database
//! and a public helper [`reannotate_modified_sequence`] that converts Sage-style mass
//! bracket notation (e.g. `[+42.0106]-MAGIR`, `C[+57.0215]`) into UniMod accession
//! notation (e.g. `.(UniMod:1)MAGIR`, `C(UniMod:4)`).
//!
//! # Embedded database
//! The standard `unimod.xml` is compiled into the binary via [`include_bytes!`] so no
//! external file is required at runtime. A custom XML file can also be loaded from disk.

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};
use quick_xml::events::Event;
use quick_xml::reader::Reader;
use regex::Regex;

// ---------------------------------------------------------------------------
// Embedded default unimod.xml
// ---------------------------------------------------------------------------
static EMBEDDED_UNIMOD_XML: &[u8] = include_bytes!("../data/unimod.xml");

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// A single UniMod modification entry: record_id -> monoisotopic delta mass,
/// keyed by (site, position).
#[derive(Debug, Clone)]
struct UnimodEntry {
    record_id: i32,
    mono_mass: f64,
}

/// Lookup table: `site` -> `position` -> `Vec<UnimodEntry>`.
///
/// Sites are single amino acid characters (`A`..`Z`), plus `"N-term"` and `"C-term"`.
/// Positions are one of `"Anywhere"`, `"Any N-term"`, `"Any C-term"`,
/// `"Protein N-term"`, or `"Protein C-term"`.
type PtmMap = HashMap<String, HashMap<String, Vec<UnimodEntry>>>;

/// Parsed UniMod database used for mass-to-accession lookups.
#[derive(Debug, Clone)]
pub struct UnimodDb {
    ptms: PtmMap,
    max_delta: f64,
}

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

impl UnimodDb {
    /// Parse the embedded default `unimod.xml` bundled with the crate.
    pub fn from_embedded(max_delta: f64) -> Result<Self> {
        Self::from_bytes(EMBEDDED_UNIMOD_XML, max_delta)
    }

    /// Parse a `unimod.xml` file from the given path.
    pub fn from_file<P: AsRef<Path>>(path: P, max_delta: f64) -> Result<Self> {
        let bytes = std::fs::read(path.as_ref())
            .with_context(|| format!("reading unimod XML from {}", path.as_ref().display()))?;
        Self::from_bytes(&bytes, max_delta)
    }

    /// Parse raw XML bytes into a [`UnimodDb`].
    pub fn from_bytes(xml_bytes: &[u8], max_delta: f64) -> Result<Self> {
        let ptms = parse_unimod_xml(xml_bytes)?;
        Ok(Self { ptms, max_delta })
    }

    // -----------------------------------------------------------------
    // Lookup
    // -----------------------------------------------------------------

    /// Find the best-matching UniMod record_id for a given `site`, set of
    /// candidate `positions`, and observed `delta_mass`.
    ///
    /// Returns `Some(record_id)` for the closest match within `max_delta`,
    /// preferring the entry with the **smallest record_id** among ties.
    /// Returns `None` if no match is found.
    pub fn get_id(&self, site: &str, positions: &[&str], delta_mass: f64) -> Option<i32> {
        let site_map = self.ptms.get(site)?;

        let mut best: Option<(i32, f64)> = None; // (record_id, delta_diff)

        for &pos in positions {
            if let Some(entries) = site_map.get(pos) {
                for entry in entries {
                    let diff = (entry.mono_mass - delta_mass).abs();
                    if diff < self.max_delta {
                        let dominated = match best {
                            Some((_, best_diff)) => {
                                diff < best_diff
                                    || (diff == best_diff && entry.record_id < best.unwrap().0)
                            }
                            None => true,
                        };
                        if dominated {
                            best = Some((entry.record_id, diff));
                        }
                    }
                }
            }
        }

        best.map(|(id, _)| id)
    }
}

// ---------------------------------------------------------------------------
// XML parsing (quick-xml streaming)
// ---------------------------------------------------------------------------

/// Low-level XML parser that builds the `PtmMap` from raw XML bytes.
fn parse_unimod_xml(xml_bytes: &[u8]) -> Result<PtmMap> {
    let ns = b"umod:";

    let mut ptms: PtmMap = HashMap::new();

    // Pre-populate all known sites and positions so lookups never fail on missing keys.
    let sites: &[&str] = &[
        "A", "R", "N", "D", "C", "E", "Q", "G", "H", "O", "I", "L", "K", "M", "F", "P", "U",
        "S", "T", "W", "Y", "V", "N-term", "C-term", "B", "J", "X", "Z",
    ];
    let positions: &[&str] = &[
        "Anywhere",
        "Any N-term",
        "Any C-term",
        "Protein N-term",
        "Protein C-term",
    ];
    for &site in sites {
        let pos_map = ptms.entry(site.to_string()).or_default();
        for &pos in positions {
            pos_map.entry(pos.to_string()).or_default();
        }
    }

    // State machine variables while scanning
    let mut reader = Reader::from_reader(xml_bytes);
    reader.config_mut().trim_text(true);

    let mut buf = Vec::new();

    // Current <umod:mod> being parsed
    let mut current_record_id: Option<i32> = None;
    let mut current_mono_mass: Option<f64> = None;
    // Collected specificities for the current mod
    let mut current_specificities: Vec<(String, String)> = Vec::new(); // (site, position)

    loop {
        match reader.read_event_into(&mut buf) {
            Ok(Event::Start(ref e)) | Ok(Event::Empty(ref e)) => {
                let local_name = e.local_name();
                let name = local_name.as_ref();

                if name == b"mod" || name == strip_ns(e.name().as_ref(), ns) && strip_ns(e.name().as_ref(), ns) == b"mod" {
                    // <umod:mod ... record_id="1" ...>
                    // We detect <umod:mod> by checking the local name == "mod"
                    let mut rid: Option<i32> = None;
                    for attr in e.attributes().flatten() {
                        if attr.key.as_ref() == b"record_id" {
                            rid = std::str::from_utf8(&attr.value).ok().and_then(|s| s.parse().ok());
                        }
                    }
                    current_record_id = rid;
                    current_mono_mass = None;
                    current_specificities.clear();
                }

                if name == b"specificity" {
                    // <umod:specificity site="C" position="Anywhere" .../>
                    let mut site: Option<String> = None;
                    let mut position: Option<String> = None;
                    for attr in e.attributes().flatten() {
                        match attr.key.as_ref() {
                            b"site" => {
                                site = std::str::from_utf8(&attr.value).ok().map(|s| s.to_string());
                            }
                            b"position" => {
                                position = std::str::from_utf8(&attr.value).ok().map(|s| s.to_string());
                            }
                            _ => {}
                        }
                    }
                    if let (Some(s), Some(p)) = (site, position) {
                        current_specificities.push((s, p));
                    }
                }

                if name == b"delta" {
                    // <umod:delta mono_mass="42.010565" ...>
                    for attr in e.attributes().flatten() {
                        if attr.key.as_ref() == b"mono_mass" {
                            current_mono_mass = std::str::from_utf8(&attr.value)
                                .ok()
                                .and_then(|s| s.parse().ok());
                        }
                    }
                }
            }
            Ok(Event::End(ref e)) => {
                let local_name = e.local_name();
                let name = local_name.as_ref();

                if name == b"mod" {
                    // Closing </umod:mod> — commit the collected data
                    if let (Some(rid), Some(mm)) = (current_record_id, current_mono_mass) {
                        for (site, position) in &current_specificities {
                            let pos_map = ptms.entry(site.clone()).or_default();
                            let entries = pos_map.entry(position.clone()).or_default();
                            entries.push(UnimodEntry {
                                record_id: rid,
                                mono_mass: mm,
                            });
                        }
                    }
                    current_record_id = None;
                    current_mono_mass = None;
                    current_specificities.clear();
                }
            }
            Ok(Event::Eof) => break,
            Err(e) => return Err(anyhow::anyhow!("Error parsing unimod XML: {}", e)),
            _ => {}
        }
        buf.clear();
    }

    Ok(ptms)
}

/// Helper: strip a namespace prefix from an element name.
fn strip_ns<'a>(name: &'a [u8], ns: &[u8]) -> &'a [u8] {
    if name.starts_with(ns) {
        &name[ns.len()..]
    } else {
        name
    }
}

// ---------------------------------------------------------------------------
// Public re-annotation API
// ---------------------------------------------------------------------------

/// Re-annotate mass-bracket modifications in a Sage-style modified peptide
/// string to UniMod accession notation.
///
/// # Transformation rules
///
/// | Input (Sage)                | Output (UniMod)               |
/// |-----------------------------|-------------------------------|
/// | `[+42.0106]-MAGIR`          | `.(UniMod:1)MAGIR`            |
/// | `M[+15.9949]AGIR`           | `M(UniMod:35)AGIR`            |
/// | `C[+57.0215]PEPTIDE`        | `C(UniMod:4)PEPTIDE`          |
///
/// N-terminal modifications (indicated by `[+mass]-` or `[+mass]` at position 0)
/// are written as `.(UniMod:X)` with a leading dot. Residue modifications are
/// written directly after the amino acid.
///
/// If `enable_unannotated` is true, mass brackets that cannot be matched are
/// kept (normalised) in the output. If false, an error is returned.
pub fn reannotate_modified_sequence(
    modified_peptide: &str,
    db: &UnimodDb,
    enable_unannotated: bool,
) -> Result<String> {
    lazy_static_regex();

    let matches: Vec<MassBracketMatch> = find_mass_brackets(modified_peptide);
    if matches.is_empty() {
        return Ok(modified_peptide.to_string());
    }

    // Process from right to left to keep character positions stable.
    let mut result = modified_peptide.to_string();
    for m in matches.into_iter().rev() {
        let is_nterm = m.is_nterm;

        let (site, positions): (&str, Vec<&str>) = if is_nterm {
            ("N-term", vec!["Any N-term", "Protein N-term"])
        } else {
            // Determine the amino acid just before the bracket.
            let aa = find_preceding_aa(&result, m.start);
            match aa {
                Some((aa_char, aa_count, total_aa)) => {
                    let aa_str = &result[aa_char..aa_char + 1];
                    if aa_count == 1 {
                        (aa_str, vec!["Anywhere", "Any N-term", "Protein N-term"])
                    } else if aa_count == total_aa {
                        (aa_str, vec!["Anywhere", "Any C-term", "Protein C-term"])
                    } else {
                        (aa_str, vec!["Anywhere"])
                    }
                }
                None => {
                    // Nothing before the bracket -> treat as N-term
                    ("N-term", vec!["Any N-term", "Protein N-term"])
                }
            }
        };

        let record_id = db.get_id(site, &positions, m.delta_mass);

        if let Some(rid) = record_id {
            let annotation = if is_nterm {
                format!(".(UniMod:{})", rid)
            } else {
                format!("(UniMod:{})", rid)
            };
            result.replace_range(m.start..m.end, &annotation);
        } else if enable_unannotated {
            // Normalise: N-term `[+42]-` -> `.[+42]`
            if is_nterm && m.has_hyphen {
                let sign = if m.delta_mass >= 0.0 { "+" } else { "" };
                let normalized = format!(".[{}{}]", sign, m.delta_mass);
                result.replace_range(m.start..m.end, &normalized);
            }
            // else: keep as-is
        } else {
            anyhow::bail!(
                "Could not annotate modification {} at site {} with delta mass {} to UniMod",
                &modified_peptide[m.start..m.end],
                site,
                m.delta_mass,
            );
        }
    }

    Ok(result)
}

// ---------------------------------------------------------------------------
// Internal regex helpers
// ---------------------------------------------------------------------------

/// Parsed mass bracket match.
struct MassBracketMatch {
    start: usize,
    end: usize,
    delta_mass: f64,
    has_hyphen: bool,
    is_nterm: bool,
}

/// Thread-local compiled regex for mass brackets.
fn lazy_static_regex() {
    // warm up the thread-local regex (no-op after first call)
    MASS_BRACKET_RE.with(|_| {});
}

thread_local! {
    /// `[+57.0215]` or `[-18.0153]`, optionally followed by `-` for N-term.
    static MASS_BRACKET_RE: Regex = Regex::new(r"\[([+-]?\d+\.?\d*)\](-)?").unwrap();
    /// Strips annotation brackets/parens for counting amino acids.
    static STRIP_ANNOT_RE: Regex = Regex::new(r"\([^)]*\)|\[[^\]]*\]").unwrap();
}

/// Find all mass-bracket matches in a modified peptide string.
fn find_mass_brackets(s: &str) -> Vec<MassBracketMatch> {
    MASS_BRACKET_RE.with(|re| {
        re.find_iter(s)
            .map(|_| ()) // just to get count; re-search with captures below
            .count(); // no-op

        re.captures_iter(s)
            .map(|cap| {
                let full_match = cap.get(0).unwrap();
                let delta_mass: f64 = cap[1].parse().unwrap_or(0.0);
                let has_hyphen = cap.get(2).is_some();
                let start = full_match.start();
                let is_nterm = start == 0 || has_hyphen;

                MassBracketMatch {
                    start,
                    end: full_match.end(),
                    delta_mass,
                    has_hyphen,
                    is_nterm,
                }
            })
            .collect()
    })
}

/// Walk backwards from `bracket_start` to find the preceding amino acid.
///
/// Returns `(byte_pos_of_aa, 1-based_aa_count, total_aa_count)` or `None`.
fn find_preceding_aa(s: &str, bracket_start: usize) -> Option<(usize, usize, usize)> {
    // Walk backwards skipping over existing annotations like (…) or […]
    let chars: Vec<(usize, char)> = s[..bracket_start].char_indices().collect();
    let mut idx = chars.len().checked_sub(1)?;

    loop {
        let ch = chars[idx].1;
        if ch == ')' || ch == ']' {
            // skip backwards over matching opener
            let close = ch;
            let open = if close == ')' { '(' } else { '[' };
            let mut depth = 1u32;
            if idx == 0 {
                return None;
            }
            idx -= 1;
            while depth > 0 {
                let c = chars[idx].1;
                if c == close {
                    depth += 1;
                } else if c == open {
                    depth -= 1;
                }
                if depth > 0 {
                    if idx == 0 {
                        return None;
                    }
                    idx -= 1;
                }
            }
            // idx now points at the opener; move one before
            if idx == 0 {
                return None;
            }
            idx -= 1;
        } else {
            break;
        }
    }

    let aa_byte_pos = chars[idx].0;

    // Count amino acids in the prefix up to (but not including) the bracket
    let aa_count = STRIP_ANNOT_RE.with(|re| {
        let prefix = &s[..bracket_start];
        let stripped = re.replace_all(prefix, "");
        stripped.chars().filter(|c| c.is_ascii_alphabetic()).count()
    });

    // Total amino acids in the full string
    let total_aa = STRIP_ANNOT_RE.with(|re| {
        let stripped = re.replace_all(s, "");
        stripped.chars().filter(|c| c.is_ascii_alphabetic()).count()
    });

    Some((aa_byte_pos, aa_count, total_aa))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn db() -> UnimodDb {
        UnimodDb::from_embedded(0.02).unwrap()
    }

    #[test]
    fn test_parse_embedded_unimod() {
        let db = db();
        // Acetyl (record_id=1) should be present for N-term
        let entries = &db.ptms["N-term"]["Any N-term"];
        assert!(
            entries.iter().any(|e| e.record_id == 1),
            "Acetyl (id=1) should exist for N-term / Any N-term"
        );
    }

    #[test]
    fn test_get_id_acetyl_nterm() {
        let db = db();
        let id = db.get_id("N-term", &["Any N-term", "Protein N-term"], 42.0106);
        assert_eq!(id, Some(1), "Acetyl should match record_id 1");
    }

    #[test]
    fn test_get_id_carbamidomethyl_cys() {
        let db = db();
        let id = db.get_id("C", &["Anywhere"], 57.0215);
        assert_eq!(id, Some(4), "Carbamidomethyl on C should match record_id 4");
    }

    #[test]
    fn test_get_id_oxidation_met() {
        let db = db();
        let id = db.get_id("M", &["Anywhere"], 15.9949);
        assert_eq!(id, Some(35), "Oxidation on M should match record_id 35");
    }

    #[test]
    fn test_get_id_no_match() {
        let db = db();
        let id = db.get_id("A", &["Anywhere"], 999.999);
        assert_eq!(id, None);
    }

    #[test]
    fn test_reannotate_nterm_acetyl() {
        let db = db();
        let result = reannotate_modified_sequence("[+42.0106]-MAGIR", &db, true).unwrap();
        assert_eq!(result, ".(UniMod:1)MAGIR");
    }

    #[test]
    fn test_reannotate_carbamidomethyl() {
        let db = db();
        let result = reannotate_modified_sequence("MAAC[+57.0215]GRVR", &db, true).unwrap();
        assert_eq!(result, "MAAC(UniMod:4)GRVR");
    }

    #[test]
    fn test_reannotate_oxidation() {
        let db = db();
        let result = reannotate_modified_sequence("M[+15.9949]AGIR", &db, true).unwrap();
        assert_eq!(result, "M(UniMod:35)AGIR");
    }

    #[test]
    fn test_reannotate_nterm_plus_internal() {
        let db = db();
        // Acetyl N-term + Carbamidomethyl on C
        let result =
            reannotate_modified_sequence("[+42.0106]-MAAC[+57.0215]GRVR", &db, true).unwrap();
        assert_eq!(result, ".(UniMod:1)MAAC(UniMod:4)GRVR");
    }

    #[test]
    fn test_reannotate_no_mods() {
        let db = db();
        let result = reannotate_modified_sequence("MAGIR", &db, true).unwrap();
        assert_eq!(result, "MAGIR");
    }

    #[test]
    fn test_reannotate_unannotated_kept() {
        let db = db();
        // Some absurd mass that won't match anything
        let result = reannotate_modified_sequence("M[+999.999]AGIR", &db, true).unwrap();
        // Should be kept as-is
        assert_eq!(result, "M[+999.999]AGIR");
    }

    #[test]
    fn test_reannotate_unannotated_error() {
        let db = db();
        let result = reannotate_modified_sequence("M[+999.999]AGIR", &db, false);
        assert!(result.is_err());
    }
}
