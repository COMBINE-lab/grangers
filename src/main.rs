use polars::lazy::dsl::concat_str;
use std::env;
use std::path::PathBuf;
use std::time::{Duration, Instant};
use tracing::info;
use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};
pub mod grangers;
use crate::grangers::{options, Grangers};
use peak_alloc::PeakAlloc;
use polars::prelude::*;
use tracing::warn;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

use clap::{builder::ArgPredicate, ArgGroup, Args, Subcommand};

/// The type of references we might create
/// to map against for quantification with
/// alevin-fry.
#[derive(Clone, Debug)]
pub enum ReferenceType {
    /// The spliced + intronic (splici) reference
    SplicedIntronic,
    /// The spliced + unspliced (splicu) reference
    SplicedUnspliced,
    Spliced,
}


fn ref_type_parser(s: &str) -> Result<ReferenceType, String> {
    match s {
        "spliced+intronic" | "splici" => Ok(ReferenceType::SplicedIntronic),
        "spliced+unspliced" | "spliceu" => Ok(ReferenceType::SplicedUnspliced),
        t => Err(format!("Do not recognize reference type {}", t)),
    }
}


#[derive(Debug, Subcommand)]
pub enum Commands {
    /// build the (expanded) reference index
    #[command(arg_required_else_help = true)]
    #[command(name = "make-spliced+intronic")]
    #[command(group(
        ArgGroup::new("reftype")
        .required(true)
        .args(["fasta", "ref_seq"])
    ))]
    MakeSplici {
        /// The path to a genome fasta file.
        #[arg(short, long, help_heading="Required Arguments", display_order = 1)]
        genome: PathBuf,


        /// The path to a GTF file.
        #[arg(short, long, help_heading="Required Arguments", display_order = 2)]
        genes: PathBuf,

        /// The path to a GTF file.
        #[arg(short, long, help_heading="Required Arguments", display_order = 2, default_value_t = 91 )]
        read_length: i64,


        /// reference genome to be used for the expanded reference construction
        #[arg(short, long, help_heading="Expanded Reference Options", display_order = 2, 
            requires_ifs([
                (ArgPredicate::IsPresent, "gtf") 
            ]),
            conflicts_with = "ref_seq")]
        fasta: Option<PathBuf>,

        /// reference GTF file to be used for the expanded reference construction
        #[arg(
            short,
            long,
            help_heading = "Expanded Reference Options",
            display_order = 3,
            requires = "fasta",
            conflicts_with = "ref_seq"
        )]
        gtf: Option<PathBuf>,

        /// the target read length the splici index will be built for
        #[arg(
            short,
            long,
            help_heading = "Expanded Reference Options",
            display_order = 4,
            requires = "fasta",
            conflicts_with = "ref_seq"
        )]
        rlen: Option<u32>,

        /// deduplicate identical sequences in pyroe when building an expanded reference  reference
        #[arg(
            long = "dedup",
            help_heading = "Expanded Reference Options",
            display_order = 5,
            requires = "fasta",
            conflicts_with = "ref_seq"
        )]
        dedup: bool,

        /// target sequences (provide target sequences directly; avoid expanded reference construction)
        #[arg(long, alias = "refseq", help_heading = "Direct Reference Options", display_order = 6,
              conflicts_with_all = ["dedup", "unspliced", "spliced", "rlen", "gtf", "fasta"])]
        ref_seq: Option<PathBuf>,

        /// path to FASTA file with extra spliced sequence to add to the index
        #[arg(
            long,
            help_heading = "Expanded Reference Options",
            display_order = 7,
            requires = "fasta",
            conflicts_with = "ref_seq"
        )]
        spliced: Option<PathBuf>,

        /// path to FASTA file with extra unspliced sequence to add to the index
        #[arg(
            long,
            help_heading = "Expanded Reference Options",
            display_order = 8,
            requires = "fasta",
            conflicts_with = "ref_seq"
        )]
        unspliced: Option<PathBuf>,

        /// use piscem instead of salmon for indexing and mapping
        #[arg(long, help_heading = "Piscem Index Options", display_order = 1)]
        use_piscem: bool,

        /// the value of m to be used to construct the piscem index (must be < k)
        #[arg(
            short = 'm',
            long = "minimizer-length",
            default_value_t = 19,
            requires = "use_piscem",
            help_heading = "Piscem Index Options",
            display_order = 2
        )]
        minimizer_length: u32,

        /// path to output directory (will be created if it doesn't exist)
        #[arg(short, long, display_order = 1)]
        output: PathBuf,

        /// overwrite existing files if the output directory is already populated
        #[arg(long, display_order = 6)]
        overwrite: bool,

        /// number of threads to use when running
        #[arg(short, long, default_value_t = 16, display_order = 2)]
        threads: u32,

        /// the value of k to be used to construct the index
        #[arg(
            short = 'k',
            long = "kmer-length",
            default_value_t = 31,
            display_order = 3
        )]
        kmer_length: u32,

        /// keep duplicated identical sequences when constructing the index
        #[arg(long, display_order = 4)]
        keep_duplicates: bool,

        /// if this flag is passed, build the sparse rather than dense index for mapping
        #[arg(
            short = 'p',
            long = "sparse",
            conflicts_with = "use_piscem",
            display_order = 5
        )]
        sparse: bool,
    },
    /// quantify a sample
    #[command(arg_required_else_help = true)]
    #[command(group(
            ArgGroup::new("filter")
            .required(true)
            .args(["expect_cells", "explicit_pl", "forced_cells", "knee", "unfiltered_pl"])
            ))]
    #[command(group(
            ArgGroup::new("input-type")
            .required(true)
            .args(["index", "map_dir"])
            ))]
    Quant {
        /// chemistry
        #[arg(short, long)]
        chemistry: String,

        /// output directory
        #[arg(short, long)]
        output: PathBuf,

        /// number of threads to use when running
        #[arg(short, long, default_value_t = 16)]
        threads: u32,

        /// path to index
        #[arg(
            short = 'i',
            long = "index",
            help_heading = "Mapping Options",
            requires_ifs([
                (ArgPredicate::IsPresent, "reads1"),
                (ArgPredicate::IsPresent, "reads2")
            ])
        )]
        index: Option<PathBuf>,

        /// comma-separated list of paths to read 1 files
        #[arg(
            short = '1',
            long = "reads1",
            help_heading = "Mapping Options",
            value_delimiter = ',',
            requires = "index",
            conflicts_with = "map_dir"
        )]
        reads1: Option<Vec<PathBuf>>,

        /// comma-separated list of paths to read 2 files
        #[arg(
            short = '2',
            long = "reads2",
            help_heading = "Mapping Options",
            value_delimiter = ',',
            requires = "index",
            conflicts_with = "map_dir"
        )]
        reads2: Option<Vec<PathBuf>>,

        /// use selective-alignment for mapping (instead of pseudoalignment with structural
        /// constraints).
        #[arg(short = 's', long, help_heading = "Mapping Options")]
        use_selective_alignment: bool,

        /// use piscem for mapping (requires that index points to the piscem index)
        #[arg(long, requires = "index", help_heading = "Mapping Options")]
        use_piscem: bool,

        /// path to a mapped output directory containing a RAD file to skip mapping
        #[arg(long = "map-dir", conflicts_with_all = ["index", "reads1", "reads2"], help_heading = "Mapping Options")]
        map_dir: Option<PathBuf>,

        /// use knee filtering mode
        #[arg(short, long, help_heading = "Permit List Generation Options")]
        knee: bool,

        /// use unfiltered permit list
        #[arg(short, long, help_heading = "Permit List Generation Options")]
        unfiltered_pl: Option<Option<PathBuf>>,

        /// use forced number of cells
        #[arg(short, long, help_heading = "Permit List Generation Options")]
        forced_cells: Option<usize>,

        /// use a filtered, explicit permit list
        #[arg(short = 'x', long, help_heading = "Permit List Generation Options")]
        explicit_pl: Option<PathBuf>,

        /// use expected number of cells
        #[arg(short, long, help_heading = "Permit List Generation Options")]
        expect_cells: Option<usize>,

        /// The expected direction/orientation of alignments in the chemistry being processed. If
        /// not provided, will default to `fw` for 10xv2/10xv3, otherwise `both`.
        #[arg(short = 'd', long, help_heading="Permit List Generation Options", value_parser = clap::builder::PossibleValuesParser::new(["fw", "rc", "both"]))]
        expected_ori: Option<String>,

        /// minimum read count threshold for a cell to be retained/processed; only used with --unfiltered-pl
        #[arg(
            long,
            help_heading = "Permit List Generation Options",
            default_value_t = 10
        )]
        min_reads: usize,

        /// transcript to gene map
        #[arg(short = 'm', long, help_heading = "UMI Resolution Options")]
        t2g_map: Option<PathBuf>,

        /// resolution mode
        #[arg(short, long, help_heading = "UMI Resolution Options", value_parser = clap::builder::PossibleValuesParser::new(["cr-like", "cr-like-em", "parsimony", "parsimony-em", "parsimony-gene", "parsimony-gene-em"]))]
        resolution: String,
    },
}



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

    let ref_typ = ReferenceType::SplicedUnspliced;

    let args: Vec<String> = env::args().collect();
    let gtf_file = PathBuf::from(args.get(1).unwrap());
    let fasta_file = PathBuf::from(args.get(2).unwrap());
    let out_dir = PathBuf::from(args.get(3).unwrap());

    // create the folder if it doesn't exist
    std::fs::create_dir_all(&out_dir)?;
    let out_fa = out_dir.join("splici_fl86.fa");

    // 1. we read the gtf file as grangers. This will make sure that the eight fields are there.
    let start = Instant::now();
    let gr = Grangers::from_gtf(gtf_file.as_path(), true)?;
    let duration: Duration = start.elapsed();
    info!("Built Grangers in {:?}", duration);
    info!("Grangers shape {:?}", gr.df().shape());

    // we get the exon df and validate it
    // this will make sure that each exon has a valid transcript ID, and the exon numbers are valid
    let mut exon_gr = gr.exons(None, true)?;

    let mut fc = exon_gr.field_columns().clone();
    let df = exon_gr.df_mut();

    // we then make sure that the gene_id and gene_name fields are not both missing
    if fc.gene_id().is_none() && fc.gene_name().is_none() {
        anyhow::bail!(
            "The input GTF file must have either gene_id or gene_name field. Cannot proceed"
        );
    } else if fc.gene_id().is_none() {
        warn!("The input GTF file do not have a gene_id field. We will use gene_name as gene_id");
        // we get gene name and rename it to gene_id
        let mut gene_id = df.column(fc.gene_name().unwrap())?.clone();
        gene_id.rename("gene_id");
        fc.update("gene_id", "gene_id")?;
        // push to the df
        df.with_column(gene_id)?;
    } else if fc.gene_name().is_none() {
        warn!(
            "The input GTF file do not have a gene_name field. We will use gene_id as gene_name."
        );
        // we get gene id and rename it to gene_name
        let mut gene_name = df.column(fc.gene_id().unwrap())?.clone();
        gene_name.rename("gene_name");
        fc.update("gene_name", "gene_name")?;
        // push to the df
        df.with_column(gene_name)?;
    }

    exon_gr.field_columns = fc;
    let gene_id_s = exon_gr.get_column_name("gene_id", false)?.to_string();
    let gene_id = gene_id_s.as_str();
    let gene_name_s = exon_gr.get_column_name("gene_name", false)?.to_string();
    let gene_name = gene_name_s.as_str();
    let transcript_id_s = exon_gr.get_column_name("transcript_id", false)?.to_string();
    let transcript_id = transcript_id_s.as_str();

    // Next, we fill the missing gene_id and gene_name fields
    if exon_gr.any_nulls(&[gene_id, gene_name], false, false)? {
        warn!("Found missing gene_id and/or gene_name; Imputing. If both missing, will impute using transcript_id; Otherwise, will impute using the existing one.");
        exon_gr.df = exon_gr
            .df
            .lazy()
            .with_columns([
                when(col(gene_id).is_null())
                    .then(
                        when(col(gene_name).is_null())
                            .then(col(transcript_id))
                            .otherwise(col(gene_name)),
                    )
                    .otherwise(col(gene_id))
                    .alias(gene_id),
                when(col(gene_name).is_null())
                    .then(
                        when(col(gene_id).is_null())
                            .then(col(transcript_id))
                            .otherwise(col(gene_id)),
                    )
                    .otherwise(col(gene_name))
                    .alias(gene_name),
            ])
            .collect()?;
    }

    // Next, we get the gene_name to id mapping
    let mut gene_id_to_name =
        exon_gr
            .df()
            .select([gene_id, gene_name])?
            .unique(None, UniqueKeepStrategy::Any, None)?;

    info!(
        "Extracting the sequence of {} transcripts",
        exon_gr.df().column("transcript_id")?.n_unique()?
    );

    // Next, we write the transcript seuqences
    exon_gr.write_transcript_sequences(&fasta_file, &out_fa, None, true, false)?;

    match ref_typ {
        ReferenceType::Spliced => {
            // do nothing
        },
        ReferenceType::SplicedIntronic => {
            // Then, we get the introns
            let mut intron_gr = exon_gr.introns(None, None, None, true)?;
        
            intron_gr.extend(86, &options::ExtendOption::Both, false)?;
        
            // Then, we merge the overlapping introns
            intron_gr = intron_gr.merge(
                &[intron_gr.get_column_name("gene_id", false)?],
                false,
                None,
                None,
            )?;
        
            intron_gr.add_order(Some(&["gene_id"]), "intron_number", Some(1), true)?;
            intron_gr.df = intron_gr
                .df
                .lazy()
                .with_column(concat_str([col("gene_id"), col("intron_number")], "-I").alias("intron_id"))
                .collect()?;
        
            intron_gr.write_sequences(
                &fasta_file,
                &out_fa,
                false,
                Some("intron_id"),
                options::OOBOption::Truncate,
                true,
            )?;
        },
        ReferenceType::SplicedUnspliced => {
            // Then, we get the introns
            let mut intron_gr = exon_gr.genes(None, true)?;
            intron_gr.write_sequences(
                &fasta_file,
                &out_fa,
                false,
                Some("gene_id"),
                options::OOBOption::Truncate,
                true,
            )?;
        },
    };

    let mut file = std::fs::File::create(out_dir.join("gene_id_to_name.tsv"))?;
    CsvWriter::new(&mut file)
        .has_header(false)
        .with_delimiter(b'\t')
        .finish(&mut gene_id_to_name)?;



    // next, we write transcripts and unspliced/itrons
    // 2. we quit if the required attributes are not valid:
    // - if transcript_id field dosn't exist
    // - if transcript_id field exists but is null for some exon features
    // - if gene_name and gene_id both are missing,
    // 3. we warn if one of gene_name and gene_id fields is missing, and copy the existing one as the missing one.
    // 4. If gene_name and/or gene_id fields contain null values, we imputet them:
    //      - If gene_name is missing, we impute it with gene_id
    //      - If gene_id is missing, we impute it with gene_name

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
