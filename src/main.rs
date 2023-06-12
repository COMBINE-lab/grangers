
use tracing::info;
use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};
pub mod grangers;
use peak_alloc::PeakAlloc;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() -> anyhow::Result<()> {
    // // Check the `RUST_LOG` variable for the logger level and
    // // respect the value found there. If this environment
    // // variable is not set then set the logging level to
    // // INFO.
    tracing_subscriber::registry()
        .with(fmt::layer())
        .with(
            EnvFilter::builder()
                .with_default_directive(LevelFilter::INFO.into())
                .from_env_lossy(),
        )
        .init();

    // let start = Instant::now();
    // gr.get_transcript_sequences(&fasta_file, None, true)?;
    // let duration: Duration = start.elapsed();
    // info!("extract transcript sequences in {:?}", duration);

    // gr.df = gr.df.head(Some(100000));
    // let start = Instant::now();
    // gr.get_sequences(&fasta_file, false, None, options::OOBOption::Skip)?;
    // let duration: Duration = start.elapsed();
    // info!("extract first 100,000 sequences in {:?}", duration);

    // let mo = options::MergeOptions::new(&["seqname", "gene_id", "transcript_id"], false, 1)?;
    // let start = Instant::now();
    // gr.merge(&["seqname", "gene_id", "transcript_id"], false, None)?;
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

    // 2. we quit if the required attributes are not valid:
    // - if transcript_id field dosn't exist
    // - if transcript_id field exists but is null for some exon features
    // - if gene_name and gene_id both are missing,
    // 3. we warn if one of gene_name and gene_id fields is missing, and copy the existing one as the missing one.
    // 4. If gene_name and/or gene_id fields contain null values, we imputet them:
    //      - If gene_name is missing, we impute it with gene_id
    //     - If gene_id is missing, we impute it with gene_name
    //     - If both are missing, we impute them with transcript_id

    // let df = df!(
    //     "seqname" => ["chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr2", "chr3"],
    //     "feature_type" => ["exon", "exon", "exon", "exon", "exon", "exon", "exon", "exon"],
    //     "start" => [1i64, 12, 1, 5, 22, 1, 5, 111],
    //     "end" => [10i64, 20, 10, 20, 30, 10, 30, 1111],
    //     "strand"=> ["+", "+", "+", "+", "+", "+", "-", "."],
    //     "gene_id" => [Some("g1"), Some("g1"), Some("g2"), Some("g2"), Some("g2"), Some("g3"), Some("g4"), None],
    // )?;

    // let mut gr = Grangers::new(df, None, None, None, IntervalType::Inclusive(1), FieldColumns::default(), false).unwrap();

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

    // let merged_gr = gr.merge(&["seqname", "gene_id", "strand"], false, None)?;

    // info!("merged gr: \n{:?}", merged_gr.df());

    // let mut flanked_gr = gr.flank(10, options::FlankOptions::default())?;
    // info!("flanked gr: \n{:?}", flanked_gr.df());

    // flanked_gr.extend(10, &options::ExtendOption::Both, false)?;
    // info!("extended and flanked gr: \n{:?}", flanked_gr.df());

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
