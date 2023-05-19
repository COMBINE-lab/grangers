use crate::grangers::reader::*;
use std::ops::{Add, Sub, Mul};
use std::path::PathBuf;
use std::time::{Duration, Instant};
use polars::export::ahash::{HashMap, HashMapExt};
use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};
use polars::prelude::*;
use polars::lazy::prelude::*;
use polars::lazy::dsl::*;
pub mod grangers;

use peak_alloc::PeakAlloc;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() -> anyhow::Result<()> {
    // Check the `RUST_LOG` variable for the logger level and
    // respect the value found there. If this environment
    // variable is not set then set the logging level to
    // INFO.
    // tracing_subscriber::registry()
    //     .with(fmt::layer())
    //     .with(
    //         EnvFilter::builder()
    //             .with_default_directive(LevelFilter::INFO.into())
    //             .from_env_lossy(),
    //     )
    //     .init();

    // let gtf_file = PathBuf::from(
    //     "/mnt/scratch4/dongze/af_example_best_practices_book/af_xmpl_run/data/refdata-gex-GRCh38-2020-A/genes/genes.gtf.toy",
    // );

    // let _fasta_file = PathBuf::from(
    //     "/mnt/scratch3/alevin_fry_submission/refs/refdata-gex-GRCh38-2020-A/fasta/genome.fa",
    // );

    // println!("Start parsing GTF");
    // let start = Instant::now();

    // let gs = gtf::GStruct::from_gtf(gtf_file.as_path(), gtf::AttributeMode::Full)?;

    // let duration: Duration = start.elapsed();
    // println!("Parsed GTF in {:?}", duration);

    // let start = Instant::now();

    // let mut gr = gs.into_grangers()?;

    // let duration = start.elapsed();
    // println!("Convert GStruct to Polars using {:?}", duration);

    // println!("{:?}",gr.df().head(Some(5)));
    // println!("{:?}",gr.strand()?.equal("+")?);
    
    let mut gs = grangers::reader::GStruct {
            seqid: vec![String::from("chr1");9],
            source: vec![String::from("HAVANA");9],
            feature_type: vec![String::from("gene"), String::from("transcript"), String::from("exon"), String::from("exon"), String::from("exon"), String::from("gene"), String::from("transcript"), String::from("exon"), String::from("exon")],
            start: vec![101, 101, 101, 121,141, 201, 201, 201, 221],
            end: vec![150, 150, 110, 130, 150, 250, 250, 210,250],
            score: vec![Some(10.0);9],
            strand: vec![Some(String::from("+")), Some(String::from("+")), Some(String::from("+")), Some(String::from("+")), Some(String::from("+")), Some(String::from("-")), Some(String::from("-")), Some(String::from("-")), Some(String::from("-"))],
            phase: vec![Some(String::from("0"));9],
            attributes: grangers::reader::gtf::Attributes::new(AttributeMode::Full,FileType::GTF)?,
            comments: vec!["comment1".to_string(), "comment2".to_string()],
            directives: Some(vec!["directive1".to_string(), "directive2".to_string()]),
    };
    let gsr = &mut gs;
    gsr.attributes.file_type = FileType::GTF;
    gsr.attributes.essential.insert(String::from("gene_id"), vec![Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g2")),Some(String::from("g2")),Some(String::from("g2")),Some(String::from("g2"))]);
    gsr.attributes.essential.insert(String::from("gene_name"), vec![Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g2")),Some(String::from("g2")),Some(String::from("g2")),Some(String::from("g2"))]);
    gsr.attributes.essential.insert(String::from("transcript_id"), vec![None,Some(String::from("t1")),Some(String::from(String::from("t1"))),Some(String::from("t1")),None,Some(String::from("t2")),Some(String::from("t2")),Some(String::from("t2")),Some(String::from("t2"))]);
    
    if let Some(extra) = &mut gsr.attributes.extra {
        extra.insert(String::from("gene_version"), vec![Some(String::from("1"));9]);
    }

    let mut gr = gs.into_grangers()?;

    // println!("{:?}",gr.df().select(["start, end"]));

    // let start_flags = gr.flank(10, true, true, false)?;

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

    // start=TRUE
    // both=FALSE 
    // ignore.strand=FALSE
    let width: i64 = -10i64;
    let start=false;
    let both = false;
    let ignore_strand = false;

    // for a record
    // if the record is in the positive strand
        // - if start is true | if ignore_strand is true
            // start = start - width
            // - if both is true
                // end = start + width - 1
            // - else 
                // end = start - 1
        // else (start is false and ignore strand is false)
            // end = end + width
            // - if both is true
                // start = end - width + 1
            // - else 
                // end = start - 1
    // if ignore strand is true | if start is true and strand is positive  
    // then pin start
        // if col("start")
            // start = start - width
        


    // else pin end
        // 
    // then the start site is the anchor
        // - if rec is in positive strand
            // start = start - witdh

            let start_time: Instant = Instant::now();

    let mut df = gr.df.clone();
    let out = df.lazy()
        .with_column(
            when(ignore_strand).then(
                lit(true)
            ).otherwise(
                col("strand")
                    .eq(lit("-"))
                    .neq(lit(start))
            )
            .alias("start_flags")
        )
        .with_column(
                // when both is true
                when(both)
                    .then(
                        // when start_flag is true
                        when(col("start_flags").eq(lit(true)))
                        .then(col("start") - lit(width).abs())
                        // when start_flag is false
                        .otherwise(col("end") - lit(width).abs() + lit(1)) 
                    )
                // when both is false
                .otherwise(
                    // if width >= 0:
                    when(width >= 0).then(
                        // tstart = all_starts[idx] - abs(width) if sf else all_ends[idx] + 1

                        when(col("start_flags").eq(lit(true)))
                            .then(
                                col("start") - lit(width)
                            )
                            .otherwise(col("end") + lit(1))
                    ).otherwise(
                        // tstart = all_starts[idx] if sf else all_ends[idx] + abs(width) + 1
                        when(col("start_flags").eq(lit(true)))
                            .then(col("start"))
                            .otherwise(col("end") + lit(width) + lit(1)))
                    ).alias("start")
        )
        .with_column(
            // new_ends.append(tstart + (width * (2 if both else 1) - 1))
            col("start").add(
                (lit(width).abs().mul(when(lit(both)).then(lit(2)).otherwise(lit(1)))).sub(lit(1))
            ).alias("end")
        )
        .select([all().exclude(["start_flags"])])
        .collect()?;

        println!("{:?}", out.column("start"));
        println!("{:?}", out.column("end"));

    let duration: Duration = start_time.elapsed();
    println!("Parsed GTF in {:?}", duration);


    // let start = Instant::now();

    // let _chromsize = fasta::get_chromsize(fasta_file.as_path())?;
    // println!("Start counting chromsize");

    // let duration = start.elapsed();
    // println!("Count chromsize using {:?}", duration);
    
    // interval operation
    // merger overlapping intervals
    // find introns 
    // 

    // get stats
    // let peak_mem = PEAK_ALLOC.peak_usage_as_gb();
    // println!("Peak Memory usage was {} GB", peak_mem);
    Ok(())
}
