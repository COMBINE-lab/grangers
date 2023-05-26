use std::env;
use std::path::PathBuf;
use std::time::{Duration, Instant};
use tracing_subscriber::prelude::*;
pub mod grangers;
// use crate::grangers::{FileFormat, FlankOptions};
use crate::grangers::{Grangers, IntervalType, MergeOptions};
use peak_alloc::PeakAlloc;
use polars::prelude::*;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() -> anyhow::Result<()> {
    let stdout_log = tracing_subscriber::fmt::layer().pretty();
    let subscriber = tracing_subscriber::Registry::default().with(stdout_log);
    tracing::subscriber::set_global_default(subscriber).unwrap();

    let args: Vec<String> = env::args().collect();
    let gtf_file = PathBuf::from(args.get(1).unwrap());

    // let _fasta_file = PathBuf::from(
    //     "/mnt/scratch3/alevin_fry_submission/refs/refdata-gex-GRCh38-2020-A/fasta/genome.fa",
    // );

    // println!("Start parsing GTF");
    // let start = Instant::now();

    // let duration: Duration = start.elapsed();
    // println!("Parsed GTF in {:?}", duration);

    let start = Instant::now();
    let mut gr = Grangers::from_gtf(gtf_file.as_path(), false)?;
    let duration: Duration = start.elapsed();
    println!("Built Grangers in {:?}", duration);
    println!("Grangers shape {:?}", gr.df().shape());

    let mo = MergeOptions::new(vec!["seqnames", "gene_id", "transcript_id"], false, 1)?;
    let start = Instant::now();
    gr.merge(&mo)?;
    let duration: Duration = start.elapsed();
    println!("Merged overlapping ranges in {:?}", duration);

    let start = Instant::now();
    gr.flank(10, None)?;
    let duration: Duration = start.elapsed();
    println!("Flanked ranges in {:?}", duration);

    let start = Instant::now();
    gr.introns("gene_id")?;
    let duration: Duration = start.elapsed();
    println!("Built intron's Grangers in {:?}", duration);

    let start = Instant::now();
    gr.extend(10, &grangers::ExtendOption::Both, false)?;
    let duration: Duration = start.elapsed();
    println!("Extended ranges in {:?}", duration);

    let start = Instant::now();
    gr.build_lapper(&mo.by)?;
    let duration: Duration = start.elapsed();
    println!("Built interval tree in {:?}", duration);

    let df = df!(
        "seqnames" => ["chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr2", "chr3"],
        "feature_type" => ["exon", "exon", "exon", "exon", "exon", "exon", "exon", "exon"],
        "start" => [1i64, 12, 1, 5, 22, 1, 5, 111],
        "end" => [10i64, 20, 10, 20, 30, 10, 30, 1111],
        "strand"=> ["+", "+", "+", "+", "+", "+", "-", "."],
        "gene_id" => [Some("g1"), Some("g1"), Some("g2"), Some("g2"), Some("g2"), Some("g3"), Some("g4"), None],
    )?;

    let mut gr = Grangers::new(df, None, None, None, IntervalType::Inclusive(1)).unwrap();
    println!("gr: {:?}", gr.df());
    let mo = MergeOptions {
        by: vec![
            "seqnames".to_string(),
            "gene_id".to_string(),
            "strand".to_string(),
        ],
        ignore_strand: false,
        slack: 1,
    };

    gr.drop_nulls(None)?;
    println!("drop_nulls' gr: \n{:?}", gr.df());

    let merged_gr = gr.merge(&mo)?;

    println!("merged gr: \n{:?}", merged_gr.df());

    let mut flanked_gr = gr.flank(10, None)?;
    println!("flanked gr: \n{:?}", flanked_gr.df());

    flanked_gr.extend(10, &grangers::ExtendOption::Both, false)?;
    println!("extended and flanked gr: \n{:?}", flanked_gr.df());

    // library(GenomicRanges)
    // gr <- GRanges(
    //     seqnames = Rle(rep("chr1", 9)),
    //     ranges = IRanges(c(101, 101, 101, 121,141, 201, 201, 201, 221), end = c(150, 150, 110, 130, 150, 250, 250, 210,250), names = head(letters, 9)),
    //     strand = Rle(strand(c(rep("+",5), rep("-", 4)))),
    //     score = 1:9,
    //     GC = seq(1, 0, length=9))
    // gr

    // pl.when().then().otherwise()
    // use schema to define data type
    // let mut schema = Schema::new();
    // Filling missing data

    // interval operation
    // merger overlapping intervals
    // find introns
    //

    // get stats
    let peak_mem = PEAK_ALLOC.peak_usage_as_gb();
    println!("Peak Memory usage was {} GB", peak_mem);
    Ok(())
}
