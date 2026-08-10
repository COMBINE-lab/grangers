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

// Both `polars` and `noodles` are part of grangers' public API surface --
// `Grangers::df` is a `polars::prelude::DataFrame`, the column accessors return
// `polars::prelude::Column`, and `GrangersSequenceCollection` is built from
// `noodles::fasta::Record`. Because a downstream crate that declares its own
// `polars`/`noodles` dependency can easily resolve a *different* version, and
// pre-1.0 crates change types across minor releases, doing so produces the
// notorious "expected `DataFrame`, found `DataFrame`" error.
//
// Re-exporting them here means consumers can write `grangers::polars::...` and be
// guaranteed to get the exact version grangers was compiled against, instead of
// hand-maintaining a matching version pin and feature list. Treat the enabled
// feature sets of these two crates as part of grangers' public API.
pub use noodles;
pub use polars;
