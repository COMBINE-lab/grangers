use std::path::PathBuf;
use std::time::{Duration, Instant};
// use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};
pub mod grangers;
use crate::grangers::FileFormat;
use crate::grangers::{Grangers, IntervalType, MergeOptions};
use peak_alloc::PeakAlloc;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() -> anyhow::Result<()> {
    let gtf_file = PathBuf::from(
        "/mnt/scratch4/dongze/af_example_best_practices_book/af_xmpl_run/data/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
    );

    // let _fasta_file = PathBuf::from(
    //     "/mnt/scratch3/alevin_fry_submission/refs/refdata-gex-GRCh38-2020-A/fasta/genome.fa",
    // );

    // println!("Start parsing GTF");
    // let start = Instant::now();

    // let duration: Duration = start.elapsed();
    // println!("Parsed GTF in {:?}", duration);

    let start = Instant::now();
    let mut gr = Grangers::from_gtf(gtf_file.as_path(), true)?;
    let duration: Duration = start.elapsed();
    println!("build grangers in {:?}", duration);

    let mo = MergeOptions::new(vec!["seqnames", "gene_id", "transcript_id"], false, 1)?;
    let start = Instant::now();

    gr.build_lapper(&mo.by)?;

    let duration: Duration = start.elapsed();
    println!("build lapper in {:?}", duration);

    let start = Instant::now();

    let introns = gr.introns("gene_id")?;

    let duration: Duration = start.elapsed();
    println!("find introns in {:?}", duration);
    println!("introns: \n{:?}", introns.df().head(Some(5)));

    let start = Instant::now();
    gr = gr.flank(10, None)?;
    let duration: Duration = start.elapsed();
    println!("flank in {:?}", duration);
    println!("flanked gr {:?}", gr.df().head(Some(5)));

    let mut gs = grangers::GStruct {
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
        attributes: grangers::reader::gtf::Attributes::new(
            grangers::reader::AttributeMode::Full,
            FileFormat::GTF,
        )?,
        misc: None,
    };

    let gsr = &mut gs;
    gsr.attributes.file_type = FileFormat::GTF;
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
            Some(String::from("t1")),
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

    let mut gr = Grangers::from_gstruct(gs, IntervalType::Inclusive(1))?;

    println!(
        "{:?}",
        gr.df().select([
            "seqnames",
            "start",
            "end",
            "strand",
            "gene_id",
            "transcript_id"
        ])?
    );

    gr.drop_nulls(None)?;
    println!("drop_nulls: \n{:?}", gr.df().head(Some(5)));

    let mut merged_gr = gr.merge(&MergeOptions::new(
        vec!["seqnames", "gene_id", "transcript_id"],
        true,
        1,
    )?)?;

    println!("merged_gr: \n{:?}", merged_gr.df().head(Some(5)));

    merged_gr.drop_nulls(None)?;
    println!("drop_nulls: \n{:?}", merged_gr.df().head(Some(5)));

    let mut flanked_gr = gr.flank(10, None)?;
    println!("flanked_gr: \n{:?}", flanked_gr.df().head(Some(5)));

    flanked_gr.extend(10, &grangers::ExtendOption::Both, false)?;
    println!("added_gr: \n{:?}", flanked_gr.df().head(Some(5)));

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
