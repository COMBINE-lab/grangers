use crate::grangers::reader::*;
use std::path::PathBuf;
use std::time::{Duration, Instant};
use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};

pub mod grangers;

use peak_alloc::PeakAlloc;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() -> anyhow::Result<()> {
    // Check the `RUST_LOG` variable for the logger level and
    // respect the value found there. If this environment
    // variable is not set then set the logging level to
    // INFO.
    tracing_subscriber::registry()
        .with(fmt::layer())
        .with(
            EnvFilter::builder()
                .with_default_directive(LevelFilter::INFO.into())
                .from_env_lossy(),
        )
        .init();

    let gtf_file = PathBuf::from(
        "/mnt/scratch3/alevin_fry_submission/refs/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
    );

    let fasta_file = PathBuf::from(
        "/mnt/scratch3/alevin_fry_submission/refs/refdata-gex-GRCh38-2020-A/fasta/genome.fa",
    );

    println!("Start parsing GTF");
    let start = Instant::now();

    let gs = gtf::GStruct::from_gtf(gtf_file.as_path(), gtf::AttributeMode::Full)?;

    let duration: Duration = start.elapsed();
    println!("Parsed GTF in {:?}", duration);

    let start = Instant::now();

    let gr = gs.to_grangers()?;

    let duration = start.elapsed();
    println!("Convert GStruct to Polars using {:?}", duration);

    println!("{:?}", gr.df().head(Some(5)));

    let start = Instant::now();

    let chromsize = fasta::get_chromsize(fasta_file.as_path())?;
    println!("Start counting chromsize");

    let duration = start.elapsed();
    println!("Count chromsize using {:?}", duration);
    
    // interval operation
    // merger overlapping intervals
    // find introns 
    // 

    // get stats
    let peak_mem = PEAK_ALLOC.peak_usage_as_gb();
    println!("Peak Memory usage was {} GB", peak_mem);
    Ok(())
}
