//! Grangers is a an (aspirationally) [GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)-like library
//! for use in [Rust](https://www.rust-lang.org/).  The goal of Grangers is to provide easy ingress
//! of genomic features into [Polars](https://pola.rs/) data frames, as well as a useful API for
//! processing and manipulation of those features.  While we believe Grangers can be useful and
//! helpful today, we are open to feedback, suggestions and ideas for improvement. If you'd like to
//! suggest some, please do so over on the [GitHub page](https://github.com/COMBINE-lab/grangers).

pub mod grangers_info;
pub mod grangers_utils;
pub mod options;
pub mod reader;
pub use grangers_info::{Grangers, GrangersRecordID, GrangersSequenceCollection};
