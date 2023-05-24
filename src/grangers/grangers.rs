use super::reader::fasta::SeqInfo;
use crate::grangers::reader;
use anyhow::{bail, Context};
use polars::lazy::dsl::col;
use polars::{lazy::prelude::*, prelude::*, series::Series};
use rust_lapper::{Interval, Lapper};
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::{
    collections::HashSet,
    ops::{Add, Mul, Sub},
};
use tracing::{info, warn};

const ESSENTIAL_FIELDS: [&str; 4] = ["seqnames", "start", "end", "strand"];

// GTF files are 1-based with closed intervals.
/// The Grangers struct contains the following fields:
/// - df: the underlying polars dataframe
/// - comments: the comments in the GTF file
/// - chromsize: the chromosome size (set as none before calling `add_chromsize`)
/// - directives: the directives in the GFF file (set as none for GTF files)
#[derive(Clone)]
pub struct Grangers {
    // file_type: Option<reader::FileType>,
    df: DataFrame,
    misc: Option<HashMap<String, Vec<String>>>,
    // comments: Option<Vec<String>>,
    seqinfo: Option<SeqInfo>,
    // directives: Option<Vec<String>>,
    // attribute_names: Option<Vec<String>>,
    lapper: Option<Lapper<u64, Vec<String>>>,
}

// IO
impl Grangers {
    /// add one column into the dataframe
    // TODO: use this in the unstranded case
    pub fn add_column<T: SeriesTrait>(&mut self, series: Series) -> anyhow::Result<()> {
        if let Err(_) = self.column("strand") {
            warn!("The dataframe does not contain the required column `seqnames`; Cannot proceed");
            self.df.with_column(series)?;
        }
        Ok(())
    }
    fn validate(df: &mut DataFrame) -> anyhow::Result<()> {
        // check essential fields
        for ef in ESSENTIAL_FIELDS {
            if let Err(_) = df.column(ef) {
                bail!(
                    "The dataframe does not contain the required column {}; Cannot proceed",
                    ef
                );
            }
        }

        // check strand
        if let Err(_) = df.column("strand") {
            warn!("The dataframe does not contain the required column seqnames; Cannot proceed");
            df.with_column(Series::new_null("strand", df.height()))?;
        }

        Ok(())
    }
    /// Instantiate a new Grangers struct according to
    /// - a range dataframe that contain the ranges of the genomic features
    /// - an optional SeqInfo struct contains the reference information (usually chromosome)
    /// - an optional HashMap<String, Vec<String>>> contains additional information.
    /// The datafram should contain the ranges of the genomic features, and the SeqInfo struct should contain the chromosome sizes. The  
    pub fn new(
        mut df: DataFrame,
        seqinfo: Option<SeqInfo>,
        misc: Option<HashMap<String, Vec<String>>>,
        lapper: Option<Lapper<u64, Vec<String>>>,
    ) -> anyhow::Result<Grangers> {
        // validate the dataframe
        Grangers::validate(&mut df)?;

        Ok(Grangers {
            df,
            misc,
            seqinfo,
            lapper,
        })
    }

    pub fn from_gstruct(gstruct: reader::GStruct) -> anyhow::Result<Grangers> {
        // create dataframe!
        // we want to make some columns categorical because of this https://docs.rs/polars/latest/polars/docs/performance/index.html
        // fields
        let mut df_vec = vec![
            Series::new("seqnames", gstruct.seqid),
            // .cast(&DataType::Categorical(None)).unwrap(),
            Series::new("source", gstruct.source),
            // .cast(&DataType::Categorical(None)).unwrap(),
            Series::new("feature_type", gstruct.feature_type),
            // .cast(&DataType::Categorical(None)).unwrap(),
            Series::new("start", gstruct.start),
            Series::new("end", gstruct.end),
            Series::new("score", gstruct.score),
            Series::new("strand", gstruct.strand),
            // .cast(&DataType::Categorical(None)).unwrap()
            Series::new("phase", gstruct.phase),
            // .cast(&DataType::Categorical(None)).unwrap(),
        ];

        //for essential attributes
        for (k, v) in gstruct.attributes.essential {
            let s = Series::new(k.as_str(), v);
            // .cast(&DataType::Categorical(None)).unwrap();
            df_vec.push(s);
        }

        // for extra attributes
        if let Some(attributes) = gstruct.attributes.extra {
            for (k, v) in attributes {
                df_vec.push(Series::new(k.as_str(), v))
            }
        }
        let df = DataFrame::new(df_vec)?;
        let gr = Grangers::new(df, None, gstruct.misc, None)?;
        Ok(gr)
    }

    // Build the Grangers struct from a GTF file.\
    // Attributes except gene_id, gene_name and transcript_id
    // will be ignored if only_essential is set to true.\
    // For a human GRh38 GTF file, the extra attributes
    // can take more than 2GB of memory.
    pub fn from_gtf(file_path: &std::path::Path, only_essential: bool) -> anyhow::Result<Grangers> {
        let am = reader::AttributeMode::from(!only_essential);
        let gstruct = reader::GStruct::from_gtf(file_path, am)?;
        let gr = Grangers::from_gstruct(gstruct)?;
        Ok(gr)
    }

    // Build the Grangers struct from a GFF file.\
    // Attributes except ID, gene_id, gene_name and transcript_id
    // will be ignored if only_essential is set to true.\
    // For a human GRh38 GTF file, the extra attributes
    // can take more than 2GB of memory.
    pub fn from_gff(file_path: &std::path::Path, only_essential: bool) -> anyhow::Result<Grangers> {
        let am = reader::AttributeMode::from(!only_essential);
        let gstruct = reader::GStruct::from_gff(file_path, am)?;
        let gr = Grangers::from_gstruct(gstruct)?;
        Ok(gr)
    }

    /// add seqinfo to the Grangers struct according to a fasta file
    pub fn add_seqinfo<T: AsRef<Path>>(&mut self, genome_file: T) -> anyhow::Result<()> {
        self.seqinfo = Some(SeqInfo::from_fasta(genome_file)?);
        Ok(())
    }

    pub fn write(&self, file_path: &std::path::Path) -> anyhow::Result<()> {
        fs::create_dir_all(file_path.parent().with_context(|| {
            format!(
                "Could not get the parent directory of the given output file path {:?}",
                file_path.as_os_str()
            )
        })?)?;
        // match self.file_type {
        //     reader::FileType::GTF => {
        //         unimplemented!()
        //     }
        //     reader::FileType::GFF => {
        //         unimplemented!()
        //     }
        // }
        unimplemented!()

        // Ok(())
    }
}

// get struct fields
impl Grangers {
    /// get the reference of the underlying dataframe
    pub fn df(&self) -> &DataFrame {
        &self.df
    }

    /// get the reference of the seqinfo
    pub fn seqinfo(&self) -> Option<&SeqInfo> {
        self.seqinfo.as_ref()
    }

    /// get the mutable reference of the underlying dataframe
    pub fn df_mut(&mut self) -> &mut DataFrame {
        &mut self.df
    }

    /// get the mutable reference to the seqinfo
    pub fn seqinfo_mut(&mut self) -> Option<&mut SeqInfo> {
        self.seqinfo.as_mut()
    }
}

// get record fields
impl Grangers {
    fn column(&self, col_name: &str) -> anyhow::Result<&Series> {
        self.df
            .column(col_name)
            .with_context(|| format!("Could not get the column {} from the dataframe.", col_name))
    }

    /// get the reference to the seqnames column
    pub fn seqid(&self) -> anyhow::Result<&Series> {
        self.column("seqnames")
    }

    /// get the reference to the source column
    pub fn start(&self) -> anyhow::Result<&Series> {
        self.column("start")
    }

    /// get the reference to the end column
    pub fn end(&self) -> anyhow::Result<&Series> {
        self.column("end")
    }

    /// get the reference to the strand column
    pub fn strand(&self) -> anyhow::Result<&Series> {
        self.column("strand")
    }

    /// get the reference to the score column
    pub fn score(&self) -> anyhow::Result<&Series> {
        self.column("score")
    }

    /// get the reference to the phase (frame) column
    pub fn phase(&self) -> anyhow::Result<&Series> {
        self.column("phase")
    }

    /// get the reference to the type column
    pub fn feature_type(&self) -> anyhow::Result<&Series> {
        self.column("feature_type")
    }

    /// get the reference to the lapper (interval tree) of the Grangers struct
    pub fn lapper(&self) -> &Option<Lapper<u64, Vec<String>>> {
        &self.lapper
    }

    /// get the start, end, and strand columns as a dataframe
    pub fn range(&self) -> anyhow::Result<DataFrame> {
        let range = self.df.select(["start", "end", "strand"])?;
        Ok(range)
    }

    /// check if the GRanges object has a column of a given name
    fn is_column(&self, col_name: &str) -> bool {
        let df = self.df();
        let mut valid = false;
        if df.get_column_names().contains(&col_name) {
            valid = true;
        }
        valid
    }
}

// implement GenomicFeatures for Grangers
impl Grangers {
    /// get the intronic sequences of each gene by the "gene_id" attribute.
    /// This function calls the gaps function internally.
    pub fn intron_by_gene(&self) -> anyhow::Result<Grangers> {
        Ok(self.clone())
    }

    /// get the intronic sequences of each gene by the "gene_id" attribute.
    /// This function calls the gaps function internally, so you need to provide .
    pub fn intron_by_transcript(&self) -> anyhow::Result<Grangers> {
        Ok(self.clone())
    }

    /// The function consumes a Grangers object, flank the genomic features of it by a given width, and returns a new Grangers.
    /// The logic for flanking - (from the BiocPy/GenomicRanges)

    /// - If `start` is `True` for a given range, the flanking occurs at the start,
    /// otherwise the end.
    /// - The `widths` of the flanks are given by the `width` parameter.
    /// The widths can be negative, in which case the flanking region is
    /// reversed so that it represents a prefix or suffix of the range.

    /// Example:
    ///     gr.flank(3, True), where x indicates a range in gr and
    ///     - indicates the resulting flanking region:\

    ///         ---xxxxxxx\

    ///     If start were FALSE, the range in gr becomes\

    ///         xxxxxxx---\

    ///     For negative width, i.e. gr.flank(x, -3, FALSE),
    ///         where * indicates the overlap between x and the result:\

    ///         xxxx***\

    ///     If both is True, then, for all ranges in x,
    ///         the flanking regions are extended into
    ///         (or out of, if width is negative) the range,
    ///         so that the result straddles the given endpoint
    ///         and has twice the width given by width.

    ///     This is illustrated below for gr.flank(3, both=TRUE):\

    ///         ---***xxxx\

    /// Args:
    ///     gr: A `Grangers` object.
    ///     width (int): width to flank by.
    ///     flank_option (FlankOption, optional):
    ///     start (bool, optional): only flank starts?. Defaults to True.
    ///     both (bool, optional): both starts and ends?. Defaults to False.
    ///     ignoreStrand (bool, optional): ignore strand?. Defaults to False.

    /// Returns:
    ///     Grangers: a new `Grangers` object with the flanked ranges.
    pub fn flank(&self, width: i64, option: Option<FlankOption>) -> anyhow::Result<Grangers> {
        let start;
        let both;
        let ignore_strand;
        if let Some(option) = option {
            start = option.start;
            both = option.both;
            ignore_strand = option.ignore_strand;
        } else {
            start = true;
            both = false;
            ignore_strand = false;
        }
        let df = self
            .df()
            .clone()
            .lazy()
            .with_column(
                when(ignore_strand)
                    .then(lit(true))
                    .otherwise(col("strand").eq(lit("-")).neq(lit(start)))
                    .alias("start_flags"),
            )
            .with_column(
                // when both is true
                when(both)
                    .then(
                        // when start_flag is true
                        when(col("start_flags").eq(lit(true)))
                            .then(col("start") - lit(width).abs())
                            // when start_flag is false
                            .otherwise(col("end") - lit(width).abs() + lit(1)),
                    )
                    // when both is false
                    .otherwise(
                        // if width >= 0:
                        when(width >= 0)
                            .then(
                                // tstart = all_starts[idx] - abs(width) if sf else all_ends[idx] + 1
                                when(col("start_flags").eq(lit(true)))
                                    .then(col("start") - lit(width))
                                    .otherwise(col("end") + lit(1)),
                            )
                            .otherwise(
                                // tstart = all_starts[idx] if sf else all_ends[idx] + abs(width) + 1
                                when(col("start_flags").eq(lit(true)))
                                    .then(col("start"))
                                    .otherwise(col("end") + lit(width) + lit(1)),
                            ),
                    )
                    .alias("start"),
            )
            .select([
                // everything except end and start_flags
                all().exclude(["end", "start_flags"]),
                // new_ends.append(tstart + (width * (2 if both else 1) - 1))
                col("start")
                    .add(
                        (lit(width)
                            .abs()
                            .mul(when(lit(both)).then(lit(2)).otherwise(lit(1))))
                        .sub(lit(1)),
                    )
                    .alias("end"),
            ])
            .collect()?;
        Ok(Grangers {
            df,
            seqinfo: self.seqinfo.clone(),
            misc: self.misc.clone(),
            lapper: self.lapper.clone(),
        })
    }

    /// Find the set difference of genomic intervals with `other`.
    /// The `on` and `boundary_on` arguments are used as anchors to find the corresponding intervals in `self` and `boundary` respectively.
    /// For the boundary dataframe, all values in the `boundary_on` column should be unique, because it is used for defining the boundary of each feature.
    // TODO: implement this after figuring out the interval inclusive/exclusive
    // The idea will be to add two more rows to the gaps result, one from 1 to the smallest start and another one from the largest end to the end of the chromosome
    pub fn setdiff(
        &self,
        _boundary: Grangers,
        _on: &str,
        _boundary_on: &str,
    ) -> anyhow::Result<Grangers> {
        unimplemented!();
    }

    /// this function turns the seqinfo of the Grangers object into a boundary Grangers object.
    pub fn seqinfo_to_bounary(&self) -> anyhow::Result<Grangers> {
        unimplemented!();
    }

    /// Find the gaps between features in each group identified by the provided `by` vector.
    /// As this function will call the `merge` function first, so it takes a `MergeOptions` as the parameter.
    pub fn gaps(&self, options: &MergeOptions) -> anyhow::Result<Grangers> {
        // check if the `by` vector contains valid column namess
        for col_name in options.by.iter() {
            if self.is_column(col_name.as_ref()) {
                bail!(
                    " `by` contains non-existing column - {}. Cannot proceed",
                    col_name
                );
            }
        }

        // merge returns a sorted and merged Grangers object
        let gr = self.merge(&options)?;
        let df = Grangers::apply(self.df(), &options.by, 0, options.ignore_strand, apply_gaps)?;
        Ok(Grangers {
            df,
            seqinfo: gr.seqinfo,
            misc: gr.misc,
            lapper: gr.lapper,
        })
    }

    /// drop rows inplace that contain missing values in the given columns.
    /// if `cols` is None, all columns will be checked.
    pub fn drop_nulls(&mut self, cols: Option<&Vec<String>>) -> anyhow::Result<()> {
        // check the validity of the column names
        let check_cols = match cols {
            Some(cols) => cols.to_owned(),
            None => self
                .df()
                .get_column_names()
                .into_iter()
                .map(|n| n.to_string())
                .collect(),
        };
        for col in check_cols.iter() {
            if !self.is_column(col) {
                bail!("Column {} does not exist in the Grangers object.", col);
            }
        }

        *self.df_mut() = self.df().drop_nulls(Some(&check_cols))?;
        Ok(())
    }

    /// merge the features by the given columns via the `by` argument to generate a new Grangers object.
    /// The `by` columns cannot have any missing value. If yours' do, run something like `gr.drop_nulls(&vec!["by_col1".to_string(), "by_col2".to_string()])` first.
    /// ### Argument:
    /// Merge Option: a struct containing the following fields:
    /// - `by`: a vector of string representing which group(s) to merge by. Each string should be a valid column name.
    /// - `ignore_strand`: whether to ignore the strand information when merging.
    /// - `slack`: the maximum distance between two features to be merged.
    pub fn merge(&self, options: &MergeOptions) -> anyhow::Result<Grangers> {
        // check if the `by` vector contains valid column namess
        for col_name in options.by.iter() {
            if !self.is_column(col_name.as_ref()) {
                bail!(
                    " `by` contains non-existing column - {}. Cannot proceed",
                    col_name
                );
            }
        }
        let df = Grangers::apply(
            self.df(),
            &options.by,
            options.slack,
            options.ignore_strand,
            apply_merge,
        )?;

        Ok(Grangers {
            df,
            seqinfo: self.seqinfo.clone(),
            misc: self.misc.clone(),
            lapper: self.lapper.clone(),
        })
    }

    fn apply<F>(
        df: &DataFrame,
        by: &Vec<String>,
        slack: i64,
        ignore_strand: bool,
        apply_fn: F,
    ) -> anyhow::Result<DataFrame>
    where
        F: Fn(Series, i64) -> Result<Option<polars::prelude::Series>, PolarsError>
            + Copy
            + std::marker::Send
            + std::marker::Sync
            + 'static,
    {
        // we take the selected columns and add two more columns: start and end
        let mut selected = by.clone();
        selected.append(vec![String::from("start"), String::from("end")].as_mut());
        // make the thing as polars friendly
        let selection = selected
            .iter()
            .map(|n| col(n.as_str()))
            .collect::<Vec<Expr>>();

        // the lazy API of polars takes the ownership of a dataframe
        let mut df = df.select(&selected)?;

        // let's see polars' way of checking missing values saying df.isna().sum()
        if df
            .null_count()
            .sum()
            .get_row(0)?
            .0
            .into_iter()
            .any(|c| c != AnyValue::UInt32(0))
        {
            warn!("Found null value(s) in the selected columns. As null will be used for grouping, we recommend dropping all null values by calling gr.drops_nulls() beforehand.")
        }

        // we will do the following
        // 1. sort the dataframe by the `by` columns + start and end columns
        // 2. group by the `by` columns
        // 3. build the new start and end columns by applying the `apply_fn` function
        // 4. explode the new start and end columns (makes the columns tall instead of wide)
        // 5. drop the intermediate columns
        df = df
            .lazy()
            .sort_by_exprs(selection, vec![false; selected.len()], false)
            .groupby(by.iter().map(|s| s.as_str()).collect::<Vec<&str>>())
            .agg([
                all().exclude(["start", "end"]).first(),
                // process two columns at once
                // Notice the df is sorted
                as_struct(&[col("start"), col("end")])
                    .apply(
                        move |s| apply_fn(s, slack),
                        GetOutput::from_type(DataType::List((DataType::Int64).into())),
                    )
                    .alias("pos"),
            ])
            .explode(["pos"])
            // with_columns returns all columns and adds extra
            // as we can't drop a non-existing column, we need to add a dummy column
            .with_columns([
                col("pos").arr().get(lit(0)).alias("start"),
                col("pos").arr().get(lit(1)).alias("end"),
                lit(NULL).cast(DataType::Utf8).alias(if ignore_strand {
                    "strand"
                } else {
                    "ignore_strand"
                }),
            ])
            .with_column(lit(NULL).cast(DataType::Utf8).alias("ignore_strand"))
            .select([all().exclude(["pos", "ignore_strand"])])
            // rearrange the columns
            .select([
                col("seqnames"),
                col("start"),
                col("end"),
                col("strand"),
                all().exclude(["seqnames", "start", "end", "strand"]),
            ])
            .collect()?;

        Ok(df)
    }

    /// merge the features by the given columns via the `by` argument to generate a new Grangers object.
    /// *** Argument:
    /// - `by`: a vector of string representing which group(s) to merge by. Each string should be a valid column name.
    /// - `ignore_strand`: whether to ignore the strand information when merging.
    /// - `slack`: the maximum distance between two features to be merged.
    // TODO: add slack to this function
    pub fn _lapper_merge(&mut self, by: Vec<String>, _slack: i64) -> anyhow::Result<Grangers> {
        self.build_lapper(&by)?;

        let mut lapper = if let Some(lapper) = self.lapper.clone() {
            lapper
        } else {
            bail!("Could not find the lapper that was just built; Please report this bug!")
        };

        lapper.merge_overlaps();

        if !lapper.overlaps_merged {
            info!("Did not find overlapping features. Nothing to merge.")
        }

        let lapper_vecs = LapperVecs::new(&lapper, &by);
        let df = lapper_vecs.to_df()?;
        Grangers::new(df, self.seqinfo.clone(), self.misc.clone(), Some(lapper))
    }

    /// Instantiate a Lapper struct for interval search and overlap detection
    /// The Lapper struct works for unsigned integer only, so the start and end values will be converted to unsigned integer.
    /// The range is [start, end), i.e. inclusive start and exclusive end. Note that for GTF/GFF, the start is 1-based and the end is closed; for BED, the start is 0-based and open but the end is closed.
    // TODO: figure out the interval ex/inclusive and slack issue and then improve it
    pub fn build_lapper(&mut self, meta_cols: &Vec<String>) -> anyhow::Result<()> {
        // rust_lapper
        // only unsigned
        // [start, stop) Inclusive start, exclusive of stop
        type Iv = Interval<u64, Vec<String>>;

        let df = self.df();
        let mut meta_vec = vec![Vec::<String>::with_capacity(meta_cols.len()); df.height()];
        let mut iters = df
            .columns(meta_cols)?
            .iter()
            .map(|s| s.iter())
            .collect::<Vec<_>>();

        // build meta_vec
        // iterate over the iterator of each selected column
        for row in 0..df.height() {
            for iter in &mut iters {
                let value = iter
                    .next()
                    .expect("should have as many iterations as rows")
                    // .cast(&DataType::Utf8)? // this will cause a bug saying cannot cast non numeric to numeric
                    .to_string();
                // process value
                if let Some(v) = meta_vec.get_mut(row) {
                    v.push(value);
                } else {
                    bail!("meta_vec length is not equal to the dataframe height! Please Report this bug!")
                }
            }
        }

        // build a lapper from the dataframe
        // check if negative start or end
        let min_start = df
            .column("start")?
            .cast(&DataType::Int64)?
            .i64()?
            .min()
            .with_context(|| {
                format!("Could not get the minimum start value from the dataframe.")
            })?;

        // let mut pos_vec = vec![Vec::<String>::with_capacity(meta_cols.len()); df.height()];
        let mut iters = df
            .columns(["start", "end"])?
            .iter()
            .map(|s| s.iter())
            .collect::<Vec<_>>();

        let mut lapper_tree_vec = Vec::with_capacity(df.height());
        // build pos_vec
        for (_rid, meta) in meta_vec.into_iter().enumerate() {
            let s: i64 = iters[0]
                .next()
                .expect("should have as many iterations as rows")
                .cast(&DataType::Int64)?
                .try_extract()?;
            let e: i64 = iters[1]
                .next()
                .expect("should have as many iterations as rows")
                .cast(&DataType::Int64)?
                .try_extract()?;

            // lappers take unsigned integer
            let (start, end) = if !min_start.is_negative() {
                (s as u64, e as u64)
            } else {
                ((s + min_start.abs()) as u64, (e + min_start.abs()) as u64)
            };

            lapper_tree_vec.push(Iv {
                start,
                stop: end,
                val: meta,
            });
        }

        // rust-lapper
        let start = std::time::Instant::now();
        self.lapper = Some(Lapper::new(lapper_tree_vec));
        let duration: std::time::Duration = start.elapsed();
        info!("build rust-lapper in {:?}", duration);
        Ok(())
    }
}

struct LapperVecs {
    start: Vec<i64>,
    end: Vec<i64>,
    val: Vec<Vec<String>>,
    val_names: Vec<String>,
}

impl LapperVecs {
    fn new(lapper: &Lapper<u64, Vec<String>>, meta_cols: &Vec<String>) -> LapperVecs {
        // initialize the vectors
        let mut start: Vec<i64> = Vec::with_capacity(lapper.len());
        let mut end: Vec<i64> = Vec::with_capacity(lapper.len());
        let mut val: Vec<Vec<String>> =
            vec![Vec::<String>::with_capacity(lapper.len()); meta_cols.len()];

        // build the vectors
        for iv in lapper.iter() {
            start.push(iv.start as i64);
            end.push(iv.stop as i64);
            for (i, v) in iv.val.iter().enumerate() {
                val.get_mut(i)
                    .expect("lapper vec val is shorter than designed")
                    .push(v.to_string());
            }
        }

        LapperVecs {
            start,
            end,
            val,
            val_names: meta_cols.clone(),
        }
    }

    fn to_df(self) -> anyhow::Result<DataFrame> {
        let mut df_vec = vec![
            Series::new("start", self.start),
            Series::new("end", self.end),
        ];

        for (name, value) in self.val_names.into_iter().zip(self.val.into_iter()) {
            df_vec.push(Series::new(name.as_str(), value));
        }

        let df = DataFrame::new(df_vec)?;
        Ok(df)
    }
}

#[derive(Copy, Clone)]
pub struct FlankOption {
    start: bool,
    both: bool,
    ignore_strand: bool,
}

impl FlankOption {
    pub fn default() -> FlankOption {
        FlankOption {
            start: true,
            both: false,
            ignore_strand: false,
        }
    }

    pub fn new(start: bool, both: bool, ignore_strand: bool) -> FlankOption {
        FlankOption {
            start,
            both,
            ignore_strand,
        }
    }
}

/// Options used for the merge function
/// - by: a vector of string representing which column(s) to merge by. Each string should be a valid column name.
/// - slack: the maximum distance between two features to be merged.
/// - output_count: whether to output the count of ranges of the merged features.
pub struct MergeOptions {
    pub by: Vec<String>,
    pub slack: i64,
    pub ignore_strand: bool,
}

impl MergeOptions {
    pub fn default() -> MergeOptions {
        MergeOptions {
            by: vec![String::from("seqnames"), String::from("strand")],
            slack: 0,
            ignore_strand: false,
        }
    }

    pub fn new<T: ToString>(
        by: Vec<T>,
        ignore_strand: bool,
        slack: i64,
    ) -> anyhow::Result<MergeOptions> {
        // avoid duplicated columns
        let mut by_hash: HashSet<String> = by.into_iter().map(|n| n.to_string()).collect();

        if by_hash.take(&String::from("start")).is_some()
            | by_hash.take(&String::from("end")).is_some()
        {
            bail!("The provided `by` vector cannot contain the start or end column")
        };

        if ignore_strand {
            if by_hash.take(&String::from("strand")).is_some() {
                warn!("Remove `strand` from the provided `by` vector as the ignored_strand flag is set. Please instantiate it manually if you intend to not include it.")
            }
        } else {
            by_hash.insert(String::from("strand"));
        }

        // add chromosome name and strand if needed
        if by_hash.insert(String::from("seqnames")) {
            warn!("Added `seqnames` to the `by` vector as it is required.")
        };

        Ok(MergeOptions {
            by: by_hash.into_iter().collect(),
            slack,
            ignore_strand,
        })
    }
}

fn apply_merge(s: Series, slack: i64) -> Result<Option<polars::prelude::Series>, PolarsError> {
    // get the two columns from the struct
    let ca: StructChunked = s.struct_()?.clone();

    // get the start and end series
    let start_series = &ca.fields()[0];
    let end_series = &ca.fields()[1];

    // downcast the `Series` to their known type and turn them into iterators
    let mut start_iter = start_series.i64()?.into_iter();
    let mut end_iter = end_series.i64()?.into_iter();

    // initialize variables for finding groups
    // we sorted the group, so the most left feature is the first one

    let (mut window_start, mut window_end) = if let (Some(Some(start)), Some(Some(end))) =
        (start_iter.next(), end_iter.next())
    {
        (start, end)
    } else {
        // this should not happen as we dropped all null values
        // rust will always use anyhow result by default
        return Result::<Option<polars::prelude::Series>, polars::prelude::PolarsError>::Ok(Some(
            Series::new_empty("pos", &DataType::List((DataType::Int64).into())),
        ));
    };
    // initialize variables for new features
    let mut out_list: Vec<Series> = Vec::with_capacity(start_series.len());
    let mut curr_count = 1;

    // iter each feature
    // we sorted the group, so the most left feature is the first one
    for (id, (start, end)) in start_iter.zip(end_iter).enumerate() {
        let (curr_start, curr_end) = if let (Some(start), Some(end)) = (start, end) {
            (start, end)
        } else {
            // rust will always use anyhow result by default
            return Result::<Option<polars::prelude::Series>, polars::prelude::PolarsError>::Ok(
                Some(Series::new_empty(
                    "pos",
                    &DataType::List((DataType::Int64).into()),
                )),
            );
        };

        // we know the df is sorted and the window starts from the leftmost feature
        // we want to check four cases:
        // 1. the feature is within the window
        // 2. the feature overlaps the window on the right end
        // 3. the window is within the feature

        if ((curr_start - slack) <= window_end) | // case 1 and 2
            ((curr_start <= window_start) & (curr_end >= window_end))
        // case 3
        {
            // extend the group
            // update group start and end
            curr_count += 1;
            // start is sorted so we only need to check end
            if curr_end > window_end {
                window_end = curr_end;
            }
        } else {
            out_list.push(Series::new(
                id.to_string().as_str(),
                [window_start, window_end],
            ));
            // out_list.push(Series::new(id.to_string().as_str(), [window_start, window_end, curr_count]));

            window_start = curr_start;
            window_end = curr_end;
            curr_count = 1;
        }
    }

    // Dont forget the last group
    if curr_count > 1 {
        out_list.push(Series::new("one more", [window_start, window_end]));
        // out_list.push(Series::new("one more", [window_start, window_end, curr_count]));
    }

    let ls = Series::new("pos", out_list);
    Result::<Option<polars::prelude::Series>, PolarsError>::Ok(Some(ls))
}

/// The apply function used for `gaps()`. The prerequesite is that the df is sorted. This is done in `apply()`.
/// The idea is that for the merged features, after dropped the first item, the `start` column, after substracting by one, can be used as the end column of the gaps.
/// The `end` column of the merged features, after dropped the last item and add 1, can be used as the start column of the gaps.
// TODO: The implementation is now assuming the intervals are inclusive. This should be changed to be more flexible.
fn apply_gaps(s: Series, _slack: i64) -> Result<Option<polars::prelude::Series>, PolarsError> {
    // get the two columns from the struct
    let ca: StructChunked = s.struct_()?.clone();

    // get the start and end series
    let start_series = &ca.fields()[0];
    let end_series = &ca.fields()[1];

    // downcast the `Series` to their known type and turn them into iterators
    let mut start_iter = start_series.i64()?.into_iter();
    let mut end_iter = end_series.i64()?.into_iter();

    // initialize variables for new features
    let mut out_list: Vec<Series> = Vec::with_capacity(start_series.len());

    // we drop the first item of the start column, this will be used as the end column of the gaps (still need to substract by 1)
    start_iter.next();

    // we iterate over the rest items in start
    for (id, next_feat_start) in start_iter.enumerate() {
        // we take the start and end out of the Option
        // the start of a gap is the end of the previous feature + 1
        // the end of a gap is the start of the current feature - 1
        let (gap_start, gap_end) = if let (Some(next_feat_start), Some(Some(prev_feat_end))) =
            (next_feat_start, end_iter.next())
        {
            // (gap start, gap end)
            (prev_feat_end + 1, next_feat_start - 1)
        } else {
            // rust will always use anyhow result by default
            return Result::<Option<polars::prelude::Series>, polars::prelude::PolarsError>::Ok(
                Some(Series::new_empty(
                    "pos",
                    &DataType::List((DataType::Int64).into()),
                )),
            );
        };

        out_list.push(Series::new(id.to_string().as_str(), [gap_start, gap_end]));
    }

    let ls = Series::new("pos", out_list);
    Result::<Option<polars::prelude::Series>, PolarsError>::Ok(Some(ls))
}

#[cfg(test)]
mod tests {
    // use polars::prelude::*;

    use super::*;
    use crate::gtf::{AttributeMode, Attributes, FileType, GStruct};
    #[test]
    fn test_graners() {
        // let df = df!(
        //     "seqid" => ["chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1"],
        //     "source" => ["HAVANA", "HAVANA", "HAVANA", "HAVANA", "HAVANA", "HAVANA", "HAVANA", "HAVANA", "HAVANA"],
        //     "feature_type" => ["gene", "transcript", "exon", "exon", "exon", "gene", "transcript", "exon", "exon"],
        //     "start" => [1, 1, 1, 21, 41, 101, 101, 101, 121],
        //     "end" => [50, 50, 10, 30, 50, 150, 150, 110,150],
        //     "score" => [100;9],
        //     "strand" => ["+", "+", "+", "+", "+", "-", "-", "-", "-"],
        //     "phase" => [0;9],
        //     "gene_id" =>["g1","g1","g1","g1","g1","g2","g2","g2","g2"],
        //     "gene_name" => ["g1","g1","g1","g1","g1","g2","g2","g2","g2"],
        //     "transcript_id" => [None,Some("t1"),Some("t1"),Some("t1"),None,Some("t2"),Some("t2"),Some("t2"),Some("t2")],
        // ).unwrap();
        // let comments = vec!["comment1".to_string(), "comment2".to_string()];
        // let directives = Some(vec!["directive1".to_string(), "directive2".to_string()]);
        // let file_type = reader::FileType::GTF;

        let mut gs = GStruct {
            seqid: vec![String::from("chr1"); 9],
            source: vec![String::from("HAVANA"); 9],
            feature_type: vec![
                String::from("gene"),
                String::from("transcript"),
                String::from("exon"),
                String::from("exon"),
                String::from("exon"),
                String::from("gene"),
                String::from("transcript"),
                String::from("exon"),
                String::from("exon"),
            ],
            start: vec![1, 1, 1, 21, 41, 101, 101, 101, 121],
            end: vec![50, 50, 10, 30, 50, 150, 150, 110, 150],
            score: vec![Some(10.0); 9],
            strand: vec![
                Some(String::from("+")),
                Some(String::from("+")),
                Some(String::from("+")),
                Some(String::from("+")),
                Some(String::from("+")),
                Some(String::from("-")),
                Some(String::from("-")),
                Some(String::from("-")),
                Some(String::from("-")),
            ],
            phase: vec![Some(String::from("0")); 9],
            attributes: Attributes::new(AttributeMode::Full, FileType::GTF).unwrap(),
            misc: Some(HashMap::new()),
        };
        let gsr = &mut gs;
        gsr.attributes.file_type = FileType::GTF;
        gsr.attributes.essential.insert(
            String::from("gene_id"),
            vec![
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g2")),
                Some(String::from("g2")),
                Some(String::from("g2")),
                Some(String::from("g2")),
            ],
        );
        gsr.attributes.essential.insert(
            String::from("gene_name"),
            vec![
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g2")),
                Some(String::from("g2")),
                Some(String::from("g2")),
                Some(String::from("g2")),
            ],
        );
        gsr.attributes.essential.insert(
            String::from("transcript_id"),
            vec![
                None,
                Some(String::from("t1")),
                Some(String::from(String::from("t1"))),
                Some(String::from("t1")),
                None,
                Some(String::from("t2")),
                Some(String::from("t2")),
                Some(String::from("t2")),
                Some(String::from("t2")),
            ],
        );

        if let Some(extra) = &mut gsr.attributes.extra {
            extra.insert(
                String::from("gene_version"),
                vec![Some(String::from("1")); 9],
            );
        }
        let _gr = Grangers::from_gstruct(gs).unwrap();

        // test builder
        // assert_eq!(gr.df(), &df);
        // assert_eq!(gr.comments(), &comments);
        // assert_eq!(gr.directives(), directives.as_ref());
        // assert!(gr.file_type.is_gtf());
    }

    #[test]
    fn test_flank() {
        let gr = get_toy_gr().unwrap();

        // test flank with default parameters
        let fo = FlankOption {
            start: true,
            both: false,
            ignore_strand: false,
        };
        let gr1 = gr.flank(10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![91, 91, 91, 111, 131, 251, 251, 211, 251];
        let end: Vec<i64> = vec![100, 100, 100, 120, 140, 260, 260, 220, 260];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        let gr1 = gr.flank(-10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![101, 101, 101, 121, 141, 241, 241, 201, 241];
        let end: Vec<i64> = vec![110, 110, 110, 130, 150, 250, 250, 210, 250];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        // test flank with default parameters and both=true
        let fo = FlankOption {
            start: true,
            both: true,
            ignore_strand: false,
        };
        let gr1 = gr.flank(10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![91, 91, 91, 111, 131, 241, 241, 201, 241];
        let end: Vec<i64> = vec![110, 110, 110, 130, 150, 260, 260, 220, 260];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        let gr1 = gr.flank(-10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![91, 91, 91, 111, 131, 241, 241, 201, 241];
        let end: Vec<i64> = vec![110, 110, 110, 130, 150, 260, 260, 220, 260];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        // test flank with start = false
        let fo = FlankOption {
            start: false,
            both: false,
            ignore_strand: false,
        };
        let gr1 = gr.flank(10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![151, 151, 111, 131, 151, 191, 191, 191, 211];
        let end: Vec<i64> = vec![160, 160, 120, 140, 160, 200, 200, 200, 220];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        let gr1 = gr.flank(-10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![141, 141, 101, 121, 141, 201, 201, 201, 221];
        let end: Vec<i64> = vec![150, 150, 110, 130, 150, 210, 210, 210, 230];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        // test flank with start = false and both=true
        let fo = FlankOption {
            start: false,
            both: true,
            ignore_strand: false,
        };
        let gr1 = gr.flank(10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![141, 141, 101, 121, 141, 191, 191, 191, 211];
        let end: Vec<i64> = vec![160, 160, 120, 140, 160, 210, 210, 210, 230];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        let gr1 = gr.flank(-10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![141, 141, 101, 121, 141, 191, 191, 191, 211];
        let end: Vec<i64> = vec![160, 160, 120, 140, 160, 210, 210, 210, 230];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        // test flank with ignore_strand: true
        let fo = FlankOption {
            start: true,
            both: false,
            ignore_strand: true,
        };
        let gr1 = gr.flank(10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![91, 91, 91, 111, 131, 191, 191, 191, 211];
        let end: Vec<i64> = vec![100, 100, 100, 120, 140, 200, 200, 200, 220];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        let gr1 = gr.flank(-10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![101, 101, 101, 121, 141, 201, 201, 201, 221];
        let end: Vec<i64> = vec![110, 110, 110, 130, 150, 210, 210, 210, 230];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        // test flank with ignore_strand: true and both=true
        let fo = FlankOption {
            start: true,
            both: true,
            ignore_strand: true,
        };
        let gr1 = gr.flank(10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![91, 91, 91, 111, 131, 191, 191, 191, 211];
        let end: Vec<i64> = vec![110, 110, 110, 130, 150, 210, 210, 210, 230];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        let gr1 = gr.flank(-10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![91, 91, 91, 111, 131, 191, 191, 191, 211];
        let end: Vec<i64> = vec![110, 110, 110, 130, 150, 210, 210, 210, 230];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );
    }

    fn get_toy_gr() -> anyhow::Result<Grangers> {
        let mut gs = GStruct {
            seqid: vec![String::from("chr1"); 9],
            source: vec![String::from("HAVANA"); 9],
            feature_type: vec![
                String::from("gene"),
                String::from("transcript"),
                String::from("exon"),
                String::from("exon"),
                String::from("exon"),
                String::from("gene"),
                String::from("transcript"),
                String::from("exon"),
                String::from("exon"),
            ],
            start: vec![101, 101, 101, 121, 141, 201, 201, 201, 221],
            end: vec![150, 150, 110, 130, 150, 250, 250, 210, 250],
            score: vec![Some(10.0); 9],
            strand: vec![
                Some(String::from("+")),
                Some(String::from("+")),
                Some(String::from("+")),
                Some(String::from("+")),
                Some(String::from("+")),
                Some(String::from("-")),
                Some(String::from("-")),
                Some(String::from("-")),
                Some(String::from("-")),
            ],
            phase: vec![Some(String::from("0")); 9],
            attributes: Attributes::new(AttributeMode::Full, FileType::GTF)?,
            misc: Some(HashMap::new()),
        };
        let gsr = &mut gs;
        gsr.attributes.file_type = FileType::GTF;
        gsr.attributes.essential.insert(
            String::from("gene_id"),
            vec![
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g2")),
                Some(String::from("g2")),
                Some(String::from("g2")),
                Some(String::from("g2")),
            ],
        );
        gsr.attributes.essential.insert(
            String::from("gene_name"),
            vec![
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g2")),
                Some(String::from("g2")),
                Some(String::from("g2")),
                Some(String::from("g2")),
            ],
        );
        gsr.attributes.essential.insert(
            String::from("transcript_id"),
            vec![
                None,
                Some(String::from("t1")),
                Some(String::from(String::from("t1"))),
                Some(String::from("t1")),
                None,
                Some(String::from("t2")),
                Some(String::from("t2")),
                Some(String::from("t2")),
                Some(String::from("t2")),
            ],
        );

        if let Some(extra) = &mut gsr.attributes.extra {
            extra.insert(
                String::from("gene_version"),
                vec![Some(String::from("1")); 9],
            );
        }

        let mut gr = Grangers::from_gstruct(gs)?;
        Ok(gr)
    }
}
