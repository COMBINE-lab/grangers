use anyhow::Context;
use noodles::gff::record::attributes;
use std::collections::HashMap;
use std::path::PathBuf;
use std::str::FromStr;
// use noodles::gff::{Record, Line};
use grangers::*;
use noodles::gff::record::Phase as gff_Phase;
use noodles::gff::record::Strand as gff_Strand;
use noodles::{gff, gtf};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::{Duration, Instant};
use polars::prelude::*;
use polars::datatypes::CategoricalType;
use crate::grangers::grangers_utils::file_line_count;
use crate::grangers::reader::*;

pub mod grangers;

use peak_alloc::PeakAlloc;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;


fn main() -> anyhow::Result<()> {
    let gtf_file = PathBuf::from(
        "/mnt/scratch3/alevin_fry_submission/refs/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
    );

        let start = Instant::now();
        
        let gs = GStruct::from_gtf(gtf_file.as_path(),GStructMode::Essential)?;

        let duration: Duration = start.elapsed();
        println!("Parsed GTF in {:?}", duration);

        let start = Instant::now();
        
        let df = gs.to_df()?;

        println!("{:?}",df);
        
        let duration = start.elapsed();
        println!("Convert GStruct to Polars using {:?}", duration);

        let peak_mem = PEAK_ALLOC.peak_usage_as_gb();
        println!("Peak Memory usage was {} GB", peak_mem);
    Ok(())
}
