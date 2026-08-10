//! Parse-path regression tests against committed GTF/GFF3 fixtures.
//!
//! The GTF/GFF reader is the part of grangers most likely to change behaviour
//! silently across a noodles upgrade: strand and phase are decoded from library
//! enums into the `+`/`-` and `0`/`1`/`2` strings the rest of the crate expects,
//! and for a while that decoding was done by matching on `format!("{:?}", ..)`
//! with a catch-all arm. These tests pin down the decoding so that class of
//! regression cannot land unnoticed.

use grangers::Grangers;
use std::path::Path;

fn gtf() -> Grangers {
    Grangers::from_gtf(Path::new("tests/data/strand_phase.gtf"), true)
        .expect("the GTF fixture should parse")
}

fn gff() -> Grangers {
    Grangers::from_gff(Path::new("tests/data/strand_phase.gff3"), true)
        .expect("the GFF3 fixture should parse")
}

/// Collects a string column as `Option<String>`, preserving row order and nulls.
fn col_opt(gr: &Grangers, name: &str) -> Vec<Option<String>> {
    gr.df()
        .column(name)
        .unwrap_or_else(|e| panic!("column {name} should exist: {e}"))
        .str()
        .expect("column should be of string type")
        .into_iter()
        .map(|v| v.map(|s| s.to_string()))
        .collect()
}

#[test]
fn gtf_row_count_and_columns() {
    let gr = gtf();
    // 15 data lines, 2 comment lines
    assert_eq!(gr.df().height(), 15);
    for c in ["seqname", "start", "end", "strand", "phase", "feature_type"] {
        assert!(
            gr.df().column(c).is_ok(),
            "expected the {c} column to be present"
        );
    }
}

#[test]
fn gtf_strand_is_decoded_to_plus_and_minus() {
    let strand = col_opt(&gtf(), "strand");
    assert!(
        strand
            .iter()
            .all(|s| matches!(s.as_deref(), Some("+") | Some("-"))),
        "every strand should decode to + or -, got {strand:?}"
    );
    // chr1 records are all +, chr2 records all -
    assert_eq!(strand[0].as_deref(), Some("+"));
    assert_eq!(strand[8].as_deref(), Some("-"));
    assert_eq!(
        strand.iter().filter(|s| s.as_deref() == Some("+")).count(),
        8
    );
    assert_eq!(
        strand.iter().filter(|s| s.as_deref() == Some("-")).count(),
        7
    );
}

/// Phase must map Zero/One/Two to "0"/"1"/"2" and stay null where absent.
/// A catch-all arm here previously turned every unrecognised phase into "0".
#[test]
fn gtf_phase_is_decoded_exactly() {
    let gr = gtf();
    let phase = col_opt(&gr, "phase");
    let ftype = col_opt(&gr, "feature_type");

    let cds: Vec<Option<String>> = ftype
        .iter()
        .zip(phase.iter())
        .filter(|(f, _)| f.as_deref() == Some("CDS"))
        .map(|(_, p)| p.clone())
        .collect();
    assert_eq!(
        cds,
        vec![
            Some("0".to_string()),
            Some("1".to_string()),
            Some("2".to_string())
        ]
    );

    // every non-CDS record in the fixture has a "." phase, which must be null
    assert!(ftype
        .iter()
        .zip(phase.iter())
        .filter(|(f, _)| f.as_deref() != Some("CDS"))
        .all(|(_, p)| p.is_none()));
}

#[test]
fn gtf_coordinates_and_attributes_survive_parsing() {
    let gr = gtf();
    let start = gr.df().column("start").unwrap().i64().unwrap();
    let end = gr.df().column("end").unwrap().i64().unwrap();
    assert_eq!(start.get(0), Some(100));
    assert_eq!(end.get(0), Some(2000));

    let gene_ids = col_opt(&gr, "gene_id");
    assert_eq!(gene_ids[0].as_deref(), Some("g1"));
    assert_eq!(gene_ids[14].as_deref(), Some("g3"));

    let names = col_opt(&gr, "gene_name");
    assert_eq!(names[0].as_deref(), Some("GENE1"));
}

#[test]
fn gff_parses_and_decodes_the_same_way() {
    let gr = gff();
    assert_eq!(gr.df().height(), 14);

    let ftype = col_opt(&gr, "feature_type");
    let phase = col_opt(&gr, "phase");
    let cds: Vec<Option<String>> = ftype
        .iter()
        .zip(phase.iter())
        .filter(|(f, _)| f.as_deref() == Some("CDS"))
        .map(|(_, p)| p.clone())
        .collect();
    assert_eq!(
        cds,
        vec![
            Some("0".to_string()),
            Some("1".to_string()),
            Some("2".to_string())
        ]
    );
}

/// GFF3 permits `.` (no strand) and `?` (unknown strand). grangers models strand
/// as `+`/`-` only, so both are coerced to `+`. The coercion is intentional; what
/// matters is that it stays total -- no nulls, nothing else leaking through --
/// and, in the reader, that it is counted and warned about rather than silent.
#[test]
fn gff_coerces_unstranded_and_unknown_strand_to_plus() {
    let gr = gff();
    let strand = col_opt(&gr, "strand");
    assert!(
        strand
            .iter()
            .all(|s| matches!(s.as_deref(), Some("+") | Some("-"))),
        "no strand should be null or exotic, got {strand:?}"
    );
    // the last four fixture rows are the `.` and `?` chr3 records
    for (i, s) in strand.iter().enumerate().skip(10) {
        assert_eq!(
            s.as_deref(),
            Some("+"),
            "row {i} (a . or ? strand record) should be coerced to +"
        );
    }
}

/// `Grangers::exons` runs the polars group-by/apply pipeline that the 0.53
/// upgrade rewrote (`explode`, `drop_nulls`, `exclude_cols`, `Expr::apply`).
#[test]
fn exons_and_introns_round_trip_through_the_polars_pipeline() {
    let gr = gtf();
    let exons = gr.exons(None, true).expect("exons should be derivable");
    assert_eq!(exons.df().height(), 6, "the fixture has six exon records");

    let introns = gr
        .introns(None, None, None, true)
        .expect("introns should be derivable");
    // t1 has three exons -> two introns; t2 has two exons -> one intron;
    // t3 has a single exon -> none.
    assert_eq!(introns.df().height(), 3);
    let start = introns.df().column("start").unwrap().i64().unwrap();
    let end = introns.df().column("end").unwrap().i64().unwrap();
    let mut spans: Vec<(i64, i64)> = (0..introns.df().height())
        .map(|i| (start.get(i).unwrap(), end.get(i).unwrap()))
        .collect();
    spans.sort();
    assert_eq!(spans, vec![(301, 799), (701, 1199), (1001, 1799)]);
}
