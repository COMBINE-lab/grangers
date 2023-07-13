use grangers::{options, Grangers};
use std::path::Path;
use std::pin::Pin;
use noodles::fasta::record::Record;
use anyhow;

#[test]
fn test_iterator() -> anyhow::Result<()> {

    let gtf = Path::new("data/refdata-gex-GRCh38-2020-A/genes/genes_small.gtf");
    let fa = Path::new("data/refdata-gex-GRCh38-2020-A/fasta/genome.fa");

    if !gtf.exists() || !fa.exists() {
        panic!("cannot run `test_iterator` integration test without the test data");
    }

    let mut gr = Grangers::from_gtf(gtf, true)?;
    let gene_id_s = gr.get_column_name("gene_id", false)?;

    let mut siter = gr.iter_sequences(fa, false, Some(&gene_id_s), options::OOBOption::Truncate)?;

    // Is this safe?
    let mut svec: Vec<Record> = Pin::into_inner(siter).collect();

    //while let Some(r) = siter.next() {
    //svec.push(r);
    //}

    println!("total records collected by iterator = {}", svec.len());

    let mut rvec: Vec<Record> = gr.get_sequences(fa, false, Some(&gene_id_s), options::OOBOption::Truncate)?
        .into_iter()
        .flatten()
        .collect();

    println!("total records collected by get_sequences = {}", rvec.len());

    svec.sort_unstable_by_key( |x| x.name().to_owned() );
    rvec.sort_unstable_by_key( |x| x.name().to_owned() );

    assert_eq!(svec, rvec);

    Ok(())
}
