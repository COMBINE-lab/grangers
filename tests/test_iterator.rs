use anyhow;
use grangers::grangers_utils::IntervalType;
use grangers::options::FieldColumns;
use grangers::{options, Grangers};
use noodles::fasta::record::Record;
use polars::prelude::*;
use std::pin::Pin;

#[test]
fn test_iterator() -> anyhow::Result<()> {
    let df = df!(
        "seqname" => ["chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1"],
        "feature_type" => ["gene", "transcript", "exon", "exon", "gene", "transcript", "exon", "exon"],
        "start" => [1i64, 1, 1, 21, 1, 1, 1, 41],
        "end" => [30i64, 30, 10, 30, 50, 50, 10, 50],
        "strand"=> ["+", "+", "+", "+", "-", "-", "-", "-"],
        "gene_id" => ["g1", "g1", "g1", "g1", "g2", "g2", "g2", "g2"],
        "transcript_id" => [None, Some("t1"), Some("t1"), Some("t1"), None, Some("t2"), Some("t2"), Some("t2")],
        "exon_id" => [None, None, Some("e1"), Some("e2"), None, None, Some("e3"), Some("e4")],
    )
    .unwrap();

    let mut gr = Grangers::new(
        df,
        None,
        None,
        IntervalType::Inclusive(1),
        FieldColumns::default(),
        false,
    )
    .unwrap();

    let gene_id_s = gr.get_column_name("gene_id", false)?;

    let fa_string = ">chr1\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";

    let siter = gr.iter_sequences_from_reader(
        std::io::Cursor::new(fa_string),
        false,
        Some(&gene_id_s),
        options::OOBOption::Truncate,
    )?;

    // Is this safe?
    let mut svec: Vec<Record> = Pin::into_inner(siter).map(|x| x.1).collect();

    // could also do this
    //while let Some(r) = siter.next() {
    //svec.push(r);
    //}

    println!("total records collected by iterator = {}", svec.len());

    let grangers::GrangersSequenceCollection {
        records: seq_vec,
        signature: _,
    } = gr.get_sequences_from_read(
        std::io::Cursor::new(fa_string),
        false,
        Some(&gene_id_s),
        options::OOBOption::Truncate,
    )?;

    let mut rvec: Vec<Record> = seq_vec.into_iter().map(|x| x.1).collect();

    println!("total records collected by get_sequences = {}", rvec.len());

    // the iterator version and the eager function call should return
    // equivalent sequences modulo order. Thus, after sorting, they
    // should compare as equal.
    svec.sort_unstable_by_key(|x| x.name().to_owned());
    rvec.sort_unstable_by_key(|x| x.name().to_owned());

    assert_eq!(svec, rvec);

    Ok(())
}
