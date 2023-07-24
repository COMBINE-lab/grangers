use anyhow;
use grangers::{options, Grangers};
use noodles::fasta::record::Record;
use std::path::Path;
use std::pin::Pin;

#[test]
fn test_iterator() -> anyhow::Result<()> {
    let gtf = Path::new("data/refdata-gex-GRCh38-2020-A/genes/genes_small.gtf");
    let fa = Path::new("data/refdata-gex-GRCh38-2020-A/fasta/genome.fa");

    if !gtf.exists() || !fa.exists() {
        panic!("cannot run `test_iterator` integration test without the test data");
    }

    let mut gr = Grangers::from_gtf(gtf, true)?;
    let gene_id_s = gr.get_column_name("gene_id", false)?;

    let siter = gr.iter_sequences(fa, false, Some(&gene_id_s), options::OOBOption::Truncate)?;

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
    } = gr.get_sequences(fa, false, Some(&gene_id_s), options::OOBOption::Truncate)?;

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
