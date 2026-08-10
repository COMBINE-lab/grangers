//! Grangers is a an (aspirationally) [GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)-like library
//! for use in [Rust](https://www.rust-lang.org/).  The goal of Grangers is to provide easy ingress
//! of genomic features into [Polars](https://pola.rs/) data frames, as well as a useful API for
//! processing and manipulation of those features.  While we believe Grangers can be useful and
//! helpful today, we are open to feedback, suggestions and ideas for improvement. If you'd like to
//! suggest some, please do so over on the [GitHub page](https://github.com/COMBINE-lab/grangers).
//!
//! # Getting started
//!
//! [`Grangers`] is the central type. Build one from a GTF or GFF3 file, then use
//! the range algebra and sequence extraction on top of it:
//!
//! ```no_run
//! use grangers::Grangers;
//! use std::path::Path;
//!
//! # fn main() -> anyhow::Result<()> {
//! // `true` keeps only the essential attribute columns.
//! let gr = Grangers::from_gtf(Path::new("annotation.gtf"), true)?;
//! println!("{} records", gr.df().height());
//!
//! let exons = gr.exons(None, true)?;
//! let introns = gr.introns(None, None, None, true)?;
//! # Ok(())
//! # }
//! ```
//!
//! The main entry points are:
//!
//! - Reading: [`Grangers::from_gtf`], [`Grangers::from_gff`].
//! - Range algebra: `exons`, `introns`, `genes`, `transcripts`, `boundary`,
//!   `flank`, `extend`, `merge`, `gaps`, `setdiff`.
//! - Sequences: `get_sequences`, `get_transcript_sequences`, `iter_sequences`,
//!   and the `write_sequences` family.
//!
//! # Use the re-exported `polars` and `noodles`
//!
//! Both crates are part of this crate's public API surface: [`Grangers::df`]
//! returns a `polars` `DataFrame`, the column accessors return `polars` `Column`s,
//! and sequences are `noodles::fasta::Record`s. Both are pre-1.0, so a version
//! differing even in its minor component is a *different, incompatible type* —
//! which surfaces as the confusing ``expected `DataFrame`, found `DataFrame` ``.
//!
//! To make that impossible, grangers re-exports exactly what it was built
//! against. Prefer these to your own dependency entries:
//!
//! ```
//! use grangers::noodles;
//! use grangers::polars::prelude::*;
//! ```
//!
//! A corollary: the set of `polars` features grangers enables is effectively part
//! of its public API, and bumping `polars` or `noodles` is a breaking change for
//! downstream crates even when nothing else changes.
//!
//! # A note on strand
//!
//! GTF and GFF3 permit `.` (no strand) and `?` (unknown strand). Grangers models
//! strand as `+`/`-` only, so such records are coerced to `+`. The number coerced
//! is reported via a [`tracing`] warning when the file is parsed.

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
