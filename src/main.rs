use std::env;
use std::path::PathBuf;
use std::time::{Duration, Instant};
use tracing::info;
use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};
pub mod grangers;
use crate::grangers::options::FieldColumns;
// use crate::grangers::{FileFormat, FlankOptions};
use crate::grangers::{options, Grangers, IntervalType};
use peak_alloc::PeakAlloc;
use polars::prelude::*;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() -> anyhow::Result<()> {
    // // Check the `RUST_LOG` variable for the logger level and
    // // respect the value found there. If this environment
    // // variable is not set then set the logging level to
    // // INFO.
    // tracing_subscriber::registry()
    //     .with(fmt::layer())
    //     .with(
    //         EnvFilter::builder()
    //             .with_default_directive(LevelFilter::INFO.into())
    //             .from_env_lossy(),
    //     )
    //     .init();

    // let args: Vec<String> = env::args().collect();
    // let gtf_file = PathBuf::from(args.get(1).unwrap());

    // // let _fasta_file = PathBuf::from(
    // //     "/mnt/scratch3/alevin_fry_submission/refs/refdata-gex-GRCh38-2020-A/fasta/genome.fa",
    // // );

    // // println!("Start parsing GTF");
    // // let start = Instant::now();

    // // let duration: Duration = start.elapsed();
    // // println!("Parsed GTF in {:?}", duration);

    // let start = Instant::now();
    // let mut gr = Grangers::from_gtf(gtf_file.as_path(), false)?;
    // let duration: Duration = start.elapsed();
    // info!("Built Grangers in {:?}", duration);
    // info!("Grangers shape {:?}", gr.df().shape());

    // let mo = options::MergeOptions::new(&["seqname", "gene_id", "transcript_id"], false, 1)?;
    // let start = Instant::now();
    // gr.merge(&mo)?;
    // let duration: Duration = start.elapsed();
    // info!("Merged overlapping ranges in {:?}", duration);

    // let start = Instant::now();
    // gr.flank(10, options::FlankOptions::default())?;
    // let duration: Duration = start.elapsed();
    // info!("Flanked ranges in {:?}", duration);

    // let start = Instant::now();
    // gr.introns(options::IntronsBy::Gene, None, true)?;
    // let duration: Duration = start.elapsed();
    // info!("Built intron's Grangers in {:?}", duration);

    // let start = Instant::now();
    // gr.extend(10, &options::ExtendOption::Both, false)?;
    // let duration: Duration = start.elapsed();
    // info!("Extended ranges in {:?}", duration);

    // let start = Instant::now();
    // gr.build_lapper(&mo.by)?;
    // let duration: Duration = start.elapsed();
    // info!("Built interval tree in {:?}", duration);

    // let df = df!(
    //     "seqname" => ["chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr2", "chr3"],
    //     "feature_type" => ["exon", "exon", "exon", "exon", "exon", "exon", "exon", "exon"],
    //     "start" => [1i64, 12, 1, 5, 22, 1, 5, 111],
    //     "end" => [10i64, 20, 10, 20, 30, 10, 30, 1111],
    //     "strand"=> ["+", "+", "+", "+", "+", "+", "-", "."],
    //     "gene_id" => [Some("g1"), Some("g1"), Some("g2"), Some("g2"), Some("g2"), Some("g3"), Some("g4"), None],
    // )?;

    // let mut gr = Grangers::new(df, None, None, None, IntervalType::Inclusive(1), FieldColumns::default()).unwrap();

    // // df.column("name")?.utf8()?.set(&df.column("name")?.is_null(), Some("."))?;
    // // println!("df: {:?}", df);

    // info!("gr: {:?}", gr.df());
    // let mo = options::MergeOptions {
    //     by: vec![
    //         "seqname".to_string(),
    //         "gene_id".to_string(),
    //         "strand".to_string(),
    //     ],
    //     ignore_strand: false,
    //     slack: 1,
    // };

    // gr.drop_nulls(None)?;
    // info!("drop_nulls' gr: \n{:?}", gr.df());

    // let merged_gr = gr.merge(&mo)?;

    // info!("merged gr: \n{:?}", merged_gr.df());

    // let mut flanked_gr = gr.flank(10, options::FlankOptions::default())?;
    // info!("flanked gr: \n{:?}", flanked_gr.df());

    // flanked_gr.extend(10, &options::ExtendOption::Both, false)?;
    // info!("extended and flanked gr: \n{:?}", flanked_gr.df());


    let mut df = df!(
        "seqname" => ["chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr2"],
        "start" => [1i64, 5, 1, 11, 22, 1, 5],
        "end" => [10i64, 10, 10, 20, 30, 10, 30],
        "strand"=> ["+", "+", "-", "-", "-", "+", "-"],
        "gene_id" => ["g1", "g1", "g2", "g2", "g2", "g3", "g4"],
    )
    .unwrap();

    let mut field_columns = FieldColumns::default();
    let interval_type = IntervalType::Inclusive(1);
    println!("{}", field_columns.is_valid(&df, true, false)?);
    field_columns.fix(&df, false)?;
    if interval_type.start_offset() != 0 {
        df.with_column(df.column(field_columns.start()).unwrap() - interval_type.start_offset())?;
    }

    if interval_type.end_offset() != 0 {
        df.with_column(df.column(field_columns.end()).unwrap() - interval_type.end_offset())?;
    }


        // instantiate a new Grangers struct
        let gr = Grangers {
            df,
            misc: None,
            seqinfo: None,
            lapper: None,
            interval_type,
            field_columns,
        };

        // validate
        gr.any_nulls(&gr.df().get_column_names(), true, true)?;






    // df.column("name")?.utf8()?.set(&df.column("name")?.is_null(), Some("."))?;

    // library(GenomicRanges)
    // gr <- GRanges(
    //     seqname = Rle(rep("chr1", 9)),
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
    info!("Peak Memory usage was {} GB", peak_mem);
    Ok(())
}
