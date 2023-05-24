use crate::grangers::reader::*;
use anyhow;
use polars;
use polars::prelude::*;
use std::ops::Sub;
use std::path::PathBuf;
use std::time::{Duration, Instant};
// use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};
pub mod grangers;
use crate::grangers::grangers::{Grangers, MergeOptions};
use peak_alloc::PeakAlloc;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() -> anyhow::Result<()> {
    let gtf_file = PathBuf::from(
        "/mnt/scratch4/dongze/af_example_best_practices_book/af_xmpl_run/data/refdata-gex-GRCh38-2020-A/genes/genes.gtf.toy",
    );

    // let _fasta_file = PathBuf::from(
    //     "/mnt/scratch3/alevin_fry_submission/refs/refdata-gex-GRCh38-2020-A/fasta/genome.fa",
    // );

    // println!("Start parsing GTF");
    // let start = Instant::now();

    let gs = gtf::GStruct::from_gtf(gtf_file.as_path(), gtf::AttributeMode::Full)?;

    // let duration: Duration = start.elapsed();
    // println!("Parsed GTF in {:?}", duration);

    let start = Instant::now();

    let gr = Grangers::from_gstruct(gs)?;
    let df = gr.df().clone();
    let df: DataFrame = df
        .lazy()
        .select([all().exclude(["start"]), col("start").sub(lit(1000000))])
        .collect()?;
    let mut gr = Grangers::new(df, None, None, None)?;

    println!(
        "{:?}",
        gr.df()
            .select(["seqnames", "start", "gene_id", "transcript_id"])?
            .tail(Some(3))
    );
    let duration: Duration = start.elapsed();
    println!("build grangers in {:?}", duration);

    let mo = MergeOptions::new(vec!["seqnames", "gene_id", "transcript_id"], false, 0)?;
    let start = Instant::now();

    gr.build_lapper(&mo.by)?;

    let duration: Duration = start.elapsed();
    println!("build lapper in {:?}", duration);

    let start = Instant::now();
    gr = gr.flank(10, None)?;
    let duration: Duration = start.elapsed();
    println!("flank in {:?}", duration);
    println!("flanked gr {:?}", gr.df().head(Some(5)));

    // let duration = start.elapsed();
    // println!("Convert GStruct to Polars using {:?}", duration);

    // println!("{:?}",gr.df().head(Some(5)));
    // println!("{:?}",gr.strand()?.equal("+")?);
    let mut gs = grangers::reader::GStruct {
        seqid: vec![String::from("chr1"); 9],
        source: vec![String::from("HAVANA"); 9],
        feature_type: vec![
            String::from("gene"),
            String::from("transcript"),
            String::from("exon"),
            String::from("exon"),
            String::from("exon"),
            String::from("gene"),
            String::from("transcript"),
            String::from("exon"),
            String::from("exon"),
        ],
        start: vec![101, 101, 101, 121, 141, 201, 201, 201, 221],
        end: vec![150, 150, 110, 130, 150, 250, 250, 210, 250],
        score: vec![Some(10.0); 9],
        strand: vec![
            Some(String::from("+")),
            Some(String::from("+")),
            Some(String::from("+")),
            Some(String::from("+")),
            Some(String::from("+")),
            Some(String::from("-")),
            Some(String::from("-")),
            Some(String::from("-")),
            Some(String::from("-")),
        ],
        phase: vec![Some(String::from("0")); 9],
        attributes: grangers::reader::gtf::Attributes::new(AttributeMode::Full, FileType::GTF)?,
        misc: None,
        // comments: vec!["comment1".to_sring(), "comment2".to_string()],
        // directives: Some(vec!["directive1".to_string(), "directive2".to_string()]),
    };

    let gsr = &mut gs;
    gsr.attributes.file_type = FileType::GTF;
    gsr.attributes.essential.insert(
        String::from("gene_id"),
        vec![
            Some(String::from("g1")),
            Some(String::from("g1")),
            Some(String::from("g1")),
            Some(String::from("g1")),
            Some(String::from("g1")),
            Some(String::from("g2")),
            Some(String::from("g2")),
            Some(String::from("g2")),
            Some(String::from("g2")),
        ],
    );
    gsr.attributes.essential.insert(
        String::from("gene_name"),
        vec![
            Some(String::from("g1")),
            Some(String::from("g1")),
            Some(String::from("g1")),
            Some(String::from("g1")),
            Some(String::from("g1")),
            Some(String::from("g2")),
            Some(String::from("g2")),
            Some(String::from("g2")),
            Some(String::from("g2")),
        ],
    );
    gsr.attributes.essential.insert(
        String::from("transcript_id"),
        vec![
            None,
            Some(String::from("t1")),
            Some(String::from(String::from("t1"))),
            Some(String::from("t1")),
            None,
            Some(String::from("t2")),
            Some(String::from("t2")),
            Some(String::from("t2")),
            Some(String::from("t2")),
        ],
    );

    if let Some(extra) = &mut gsr.attributes.extra {
        extra.insert(
            String::from("gene_version"),
            vec![Some(String::from("1")); 9],
        );
    }

    let gr = Grangers::from_gstruct(gs)?;

    println!(
        "{:?}",
        gr.df()
            .select(["seqnames", "start", "end", "gene_id", "transcript_id"])?
    );

    let mut merged_gr = gr.merge(&MergeOptions::new(
        vec!["seqnames", "gene_id", "transcript_id"],
        false,
        0,
    )?)?;

    println!("merged_gr: \n{:?}", merged_gr.df().head(Some(5)));

    merged_gr.drop_nulls(None)?;
    println!("drop_nulls: \n{:?}", merged_gr.df().head(Some(5)));

    let flanked_gr = gr.flank(10, None)?;
    println!("flanked_gr: \n{:?}", flanked_gr.df().head(Some(5)));

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
    // let peak_mem = PEAK_ALLOC.peak_usage_as_gb();
    // println!("Peak Memory usage was {} GB", peak_mem);
    Ok(())
}
