# grangers

[![crates.io](https://img.shields.io/crates/v/grangers.svg)](https://crates.io/crates/grangers)

A Rust library for parsing and manipulating genomic annotations, in the spirit of
Bioconductor's [GenomicRanges]/[GenomicFeatures]. `grangers` reads GTF and GFF3
files into a [polars] `DataFrame`, gives you a range algebra over the result, and
can extract the corresponding sequences from a genome FASTA.

## Installing

```toml
[dependencies]
grangers = "0.6"
```

## A quick tour

```rust
use grangers::{options, Grangers};
use std::path::Path;

# fn main() -> anyhow::Result<()> {
// Read an annotation. `true` keeps only the essential attribute columns
// (gene_id, gene_name, transcript_id, exon_number); `false` keeps them all.
let gr = Grangers::from_gtf(Path::new("annotation.gtf"), true)?;

// The underlying polars DataFrame is public.
println!("{} records", gr.df().height());

// Range algebra: derive exons, then the introns between them.
let exons = gr.exons(None, true)?;
let introns = gr.introns(None, None, None, true)?;

// Pull the spliced transcript sequences out of a genome FASTA.
let seqs = gr.get_transcript_sequences(Path::new("genome.fa"), None, true)?;
# Ok(())
# }
```

`Grangers` also offers `genes`, `transcripts`, `boundary`, `flank`, `extend`,
`merge`, `gaps` and `setdiff`, plus filtering, sorting and interval-overlap
queries backed by [rust-lapper]. Sequence extraction is available in both
collect-everything (`get_sequences`), streaming (`iter_sequences`) and
write-straight-to-file (`write_sequences`) forms.

## Working with polars and noodles

`grangers` re-exports the exact versions of its two load-bearing dependencies:

```rust
use grangers::polars::prelude::*;
use grangers::noodles;
```

**Prefer these over your own `polars`/`noodles` dependency entries.** Both types
cross the `grangers` API boundary — `Grangers::df` is a `polars` `DataFrame`, and
sequences are `noodles::fasta::Record`s — and both crates are pre-1.0, so a
version that differs even in its minor component produces a distinct, incompatible
type. The symptom is the memorable ``expected `DataFrame`, found `DataFrame` ``.

Because of this, the polars feature set `grangers` enables is effectively part of
its public API, and a `polars` or `noodles` bump is a breaking change for
consumers even when nothing else changes.

## Note on strand

GTF and GFF3 both allow `.` (no strand) and `?` (unknown strand). `grangers`
models strand as `+`/`-` only, so such records are coerced to `+`; the number
coerced is reported through a `tracing` warning at parse time.

## Related projects

- [roers](https://github.com/COMBINE-lab/roers) — builds augmented (splici /
  spliceu) transcriptome references; the main consumer of this crate.
- [simpleaf](https://github.com/COMBINE-lab/simpleaf) — end-to-end single-cell
  processing built on [alevin-fry].

## License

BSD 3-Clause; see [LICENSE](LICENSE).

[GenomicRanges]: https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
[GenomicFeatures]: https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html
[polars]: https://pola.rs/
[rust-lapper]: https://crates.io/crates/rust-lapper
[alevin-fry]: https://github.com/COMBINE-lab/alevin-fry
