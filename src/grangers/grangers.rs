use crate::grangers::reader;
use crate::grangers::reader::fasta::SeqInfo;
use anyhow::{bail, Context};
use polars::lazy::dsl::col;
use polars::{lazy::prelude::*, prelude::*, series::Series};
use rust_lapper::{Interval, Lapper};
use std::collections::HashMap;
use std::default::Default;
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
/// - misc: the additional information
/// - seqinfo: the reference information
/// - lapper: the lapper interval tree
///
/// **Notice** that Granges uses 1-based closed intervals for the ranges.
/// If your ranges are not like this, when instantiating new Grangers,
/// you should use the `interval_type` parameter to help the builder
/// to convert the ranges to 1-based closed intervals.
#[derive(Clone)]
pub struct Grangers {
    /// The underlying dataframe
    df: DataFrame,
    /// The additional information
    misc: Option<HashMap<String, Vec<String>>>,
    /// The reference information
    seqinfo: Option<SeqInfo>,
    /// The lapper interval tree
    lapper: Option<Lapper<u64, Vec<String>>>,
    /// The interval type
    interval_type: IntervalType,
}

// IO
impl Grangers {
    /// add or replace a column in the dataframe
    // TODO: use this in the unstranded case
    pub fn add_column<T: SeriesTrait>(&mut self, series: Series) -> anyhow::Result<()> {
        self.df.with_column(series)?;

        Ok(())
    }
    fn validate(df: &mut DataFrame) -> anyhow::Result<()> {
        // check essential fields
        for ef in ESSENTIAL_FIELDS {
            if df.column(ef).is_err() {
                bail!(
                    "The dataframe does not contain the required column {}; Cannot proceed",
                    ef
                );
            }
        }

        // check strand
        if df.column("strand").is_err() {
            warn!("The dataframe does not contain the strand column; adding a null strand column.");
            df.with_column(Series::new_null("strand", df.height()))?;
        }

        Ok(())
    }
    /// Instantiate a new Grangers struct according to
    /// - a range dataframe that contain the ranges of the genomic features
    /// - an optional SeqInfo struct contains the reference information (usually chromosome)
    /// - an optional HashMap<String, Vec<String>>> contains additional information.
    /// - an optional Lapper interval tree
    /// - an optional IntervalType representing the coordinate system and interval type of the ranges. If passing None, the default
    /// The datafram should contain the ranges of the genomic features, and the SeqInfo struct should contain the chromosome sizes. The  
    pub fn new(
        mut df: DataFrame,
        seqinfo: Option<SeqInfo>,
        misc: Option<HashMap<String, Vec<String>>>,
        lapper: Option<Lapper<u64, Vec<String>>>,
        interval_type: IntervalType,
    ) -> anyhow::Result<Grangers> {
        // validate the dataframe
        Grangers::validate(&mut df)?;

        if interval_type.start_offset() != 0 {
            df.with_column(df.column("start").unwrap() - interval_type.start_offset())?;
        }

        if interval_type.end_offset() != 0 {
            df.with_column(df.column("end").unwrap() - interval_type.end_offset())?;
        }

        Ok(Grangers {
            df,
            misc,
            seqinfo,
            lapper,
            interval_type,
        })
    }

    pub fn from_gstruct(
        gstruct: reader::GStruct,
        interval_type: IntervalType,
    ) -> anyhow::Result<Grangers> {
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
        let gr = Grangers::new(df, None, gstruct.misc, None, interval_type)?;
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
        let gr = Grangers::from_gstruct(gstruct, IntervalType::Inclusive(1))?;
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
        let gr = Grangers::from_gstruct(gstruct, IntervalType::Inclusive(1))?;
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
        //     reader::FileFormat::GTF => {
        //         unimplemented!()
        //     }
        //     reader::FileFormat::GFF => {
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

    /// get the interval type
    pub fn interval_type(&self) -> &IntervalType {
        &self.interval_type
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

    /// sort the dataframe
    pub fn sort_df_by<T>(&mut self, by: &[&str], descending: Vec<bool>) -> anyhow::Result<()> {
        self.df.sort(by, descending)?;
        Ok(())
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
    /// get the intronic sequences of each gene or transcript according to the given column name.
    /// The `by` parameter should specify the column name of the gene or transcript ID.
    pub fn introns<T: ToString>(&self, by: T) -> anyhow::Result<Grangers> {
        let by_string = by.to_string();
        let by = by_string.as_str();
        if !self.is_column(by) {
            anyhow::bail!("The column {} does not exist in the Grangers struct.", by);
        }

        // we requires that the "feature_type" column exists, which
        // defines the type of each genomic feature, i.e., gene, transcript, exon, UTR etc
        if !self.is_column("feature_type") {
            anyhow::bail!("The column feature_type does not exist in the Grangers struct. This is required for identifying exon records. Cannot proceed");
        }

        let mo = MergeOptions::new(vec![by], false, 1)?;
        let mut gr = self.clone();

        // polars way to subset
        let anyvalue_exon = AnyValue::Utf8("exon");

        gr.df = gr.df.filter(
            &gr.df
                .column("feature_type")?
                .iter()
                .map(|f| f.eq(&anyvalue_exon))
                .collect(),
        )?;

        if gr.column(by)?.null_count() > 0 {
            anyhow::bail!("The column {} contains null values in the by features. Cannot use it to find introns. Please fill up null values and try again.", by);
        }

        gr.gaps(&mo)
    }

    /// extend each genomic feature by a given length from the start, end, or both sides.
    pub fn extend(
        &mut self,
        length: i32,
        extend_option: &ExtendOption,
        ignore_strand: bool,
    ) -> anyhow::Result<()> {
        // if contains null value in strand, we cannot do strand-specific extension
        if (!ignore_strand) & (extend_option != &ExtendOption::Both)
            && self.column("strand")?.is_null().any()
                | self
                    .column("strand")?
                    .unique()?
                    .is_in(&Series::new("missing strand", ["*"]))?
                    .any()
        {
            bail!("The strand column contains null values. Please remove them first or set ignore_strand to true.")
        }

        if let ExtendOption::Both = extend_option {
            self.df
                .with_column(self.df.column("start")?.clone() - length)?;
            self.df
                .with_column(self.df.column("end")?.clone() + length)?;
            return Ok(());
        }

        if ignore_strand {
            match extend_option {
                ExtendOption::Start => {
                    self.df
                        .with_column(self.df.column("start")?.clone() - length)?;
                    return Ok(());
                }
                ExtendOption::End => {
                    self.df
                        .with_column(self.df.column("end")?.clone() + length)?;
                    return Ok(());
                }
                _ => {}
            }
        } else {
            let mut df = self.df().select(["start", "end", "strand"])?;
            df = df
                .lazy()
                .with_columns([
                    // we first consider the start site
                    // when the strand is + and extend from start, or the strand is - and extend from end, we extend the start site
                    when(
                        col("strand")
                            .eq(lit("+"))
                            .eq(lit(extend_option == &ExtendOption::Start))
                            .or(col("strand")
                                .eq(lit("-"))
                                .eq(lit(extend_option == &ExtendOption::End))),
                    )
                    .then(col("start").sub(lit(length)))
                    .otherwise(col("start"))
                    .alias("start"),
                    // then the end site
                    // when the strand is - and extend from start, or the strand is + and extend from end, we extend the end site
                    when(
                        col("strand")
                            .eq(lit("-"))
                            .eq(lit(extend_option == &ExtendOption::Start))
                            .or(col("strand")
                                .eq(lit("+"))
                                .eq(lit(extend_option == &ExtendOption::End))),
                    )
                    .then(col("end").add(lit(length)))
                    .otherwise(col("end"))
                    .alias("end"),
                ])
                .collect()?;

            // Then we update df
            self.df.with_column(df.column("start")?.clone())?;
            self.df.with_column(df.column("end")?.clone())?;
        }
        Ok(())
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
    ///     flank_option (FlankOptions, optional):
    ///     start (bool, optional): only flank starts?. Defaults to True.
    ///     both (bool, optional): both starts and ends?. Defaults to False.
    ///     ignoreStrand (bool, optional): ignore strand?. Defaults to False.

    /// Returns:
    ///     Grangers: a new `Grangers` object with the flanked ranges.
    // flank doesn't not related to interval thing
    pub fn flank(&self, width: i64, options: Option<FlankOptions>) -> anyhow::Result<Grangers> {
        let start;
        let both;
        let ignore_strand;
        if let Some(options) = options {
            start = options.start;
            both = options.both;
            ignore_strand = options.ignore_strand;
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
            .select(
                self.df()
                    .get_column_names()
                    .iter()
                    .map(|x| col(x))
                    .collect::<Vec<Expr>>(),
            )
            .collect()?;
        Ok(Grangers {
            df,
            seqinfo: self.seqinfo.clone(),
            misc: self.misc.clone(),
            lapper: self.lapper.clone(),
            interval_type: self.interval_type,
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
            if !self.is_column(col_name.as_ref()) {
                bail!(
                    " `by` contains non-existing column - {}. Cannot proceed",
                    col_name
                );
            }
        }

        // merge returns a sorted and merged Grangers object
        let mut gr = self.merge(options)?;
        gr.df = Grangers::apply(gr.df(), &options.by, 0, options.ignore_strand, apply_gaps)?;
        Ok(gr)
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
    /// - `slack`: the maximum distance between two features to be merged. Slack should be a non-negative integer. For example we have three intervals, [1,2], [3,4] and [4,5]
    ///     - slack = 1 means [1,2] and [3,4] will be merged as [1,4]. This is the default (and desired) behavior unless you do not want to merge adjacent intervals.
    ///     - slack = 0 means [1,2] and [2,3] will stay separated although they are adjacent.
    ///     - slack=2 means[1,2] and [4,5] will be merged though there is a base [3,3] in between that separate them. The slack here is equivalent to the `min.gapwidth` argument in the GenomicRanges::reduce() function in R.
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
        let mut by = options.by.to_owned();
        if !options.ignore_strand & !options.by.contains(&"strand".to_string()) {
            warn!("ignore_strand is set to false. Added `strand` to the `by` vector");
            by.push(String::from("strand"));
        }

        let df = Grangers::apply(
            self.df(),
            &by,
            options.slack,
            options.ignore_strand,
            apply_merge,
        )?;

        Ok(Grangers {
            df,
            seqinfo: self.seqinfo.clone(),
            misc: self.misc.clone(),
            lapper: self.lapper.clone(),
            interval_type: self.interval_type,
        })
    }

    fn apply<F>(
        df: &DataFrame,
        by: &[String],
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
        let mut selected = by.to_owned();
        selected.append(vec![String::from("start"), String::from("end")].as_mut());
        // make the thing as polars friendly
        // let selection = selected
        //     .iter()
        //     .map(|n| col(n.as_str()))
        //     .collect::<Vec<Expr>>();

        // we want to sort the dataframe by first essential columns, then the rest columns in `by`
        let mut sorted_by_exprs_essential = vec![col("seqnames"), col("start"), col("end")];
        let mut sorted_by_desc_essential = vec![false, false, true];
        if !ignore_strand {
            sorted_by_exprs_essential.push(col("strand"));
            sorted_by_desc_essential.push(false);
        }

        let mut sorted_by_exprs: Vec<Expr> = by
            .iter()
            .filter(|&n| !sorted_by_exprs_essential.contains(&col(n.as_str())))
            .map(|n| col(n.as_str()))
            .collect();

        let mut sorted_by_desc = vec![false; sorted_by_exprs.len()];
        sorted_by_exprs.extend(sorted_by_exprs_essential);
        sorted_by_desc.extend(sorted_by_desc_essential);

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
            .sort_by_exprs(&sorted_by_exprs, &sorted_by_desc, false)
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
                    // .arr().sort(SortOptions {descending: if ignore_strand {false} else {col("strand").first() == lit("-")}, nulls_last: false, multithreaded:  true})
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
            .drop_nulls(Some(vec![cols(["start", "end"])]))
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
            // groupby is multithreaded, so the order do not preserve
            .sort_by_exprs(sorted_by_exprs, sorted_by_desc, false)
            .collect()?;

        Ok(df)
    }

    /// merge the features by the given columns via the `by` argument to generate a new Grangers object.
    /// *** Argument:
    /// - `by`: a vector of string representing which group(s) to merge by. Each string should be a valid column name.
    /// - `ignore_strand`: whether to ignore the strand information when merging.
    /// - `slack`: the maximum distance between two features to be merged.
    // TODO: add slack to this function
    pub fn lapper_merge(&mut self, by: Vec<String>, _slack: i64) -> anyhow::Result<Grangers> {
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
        let df = lapper_vecs.into_df()?;
        Ok(Grangers {
            df,
            seqinfo: self.seqinfo.clone(),
            misc: self.misc.clone(),
            lapper: Some(lapper),
            interval_type: self.interval_type,
        })
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
                "Could not get the minimum start value from the dataframe.".to_string()
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

    fn into_df(self) -> anyhow::Result<DataFrame> {
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
pub struct FlankOptions {
    start: bool,
    both: bool,
    ignore_strand: bool,
}

impl Default for FlankOptions {
    fn default() -> FlankOptions {
        FlankOptions {
            start: true,
            both: false,
            ignore_strand: false,
        }
    }
}

impl FlankOptions {
    pub fn new(start: bool, both: bool, ignore_strand: bool) -> FlankOptions {
        FlankOptions {
            start,
            both,
            ignore_strand,
        }
    }
}

/// Options used for the merge function
/// - by: a vector of string representing which column(s) to merge by. Each string should be a valid column name.
/// - slack: the minimum gap between two features to be merged. It should be a positive integer.
/// - output_count: whether to output the count of ranges of the merged features.
pub struct MergeOptions {
    pub by: Vec<String>,
    pub slack: i64,
    pub ignore_strand: bool,
}

impl Default for MergeOptions {
    fn default() -> MergeOptions {
        MergeOptions {
            by: vec![String::from("seqnames"), String::from("strand")],
            slack: 1,
            ignore_strand: false,
        }
    }
}

impl MergeOptions {
    pub fn new<T: ToString>(
        by: Vec<T>,
        ignore_strand: bool,
        slack: i64,
    ) -> anyhow::Result<MergeOptions> {
        // avoid duplicated columns
        let mut by_hash: HashSet<String> = by.into_iter().map(|n| n.to_string()).collect();

        if slack < 1 {
            warn!("It usually doen't make sense to set a non-positive slack.")
        }

        if by_hash.take(&String::from("start")).is_some()
            | by_hash.take(&String::from("end")).is_some()
        {
            bail!("The provided `by` vector cannot contain the start or end column")
        };

        if ignore_strand {
            if by_hash.take(&String::from("strand")).is_some() {
                warn!("Remove `strand` from the provided `by` vector as the ignored_strand flag is set.")
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

    let (mut window_start, mut window_end) =
        if let (Some(Some(start)), Some(Some(end))) = (start_iter.next(), end_iter.next()) {
            (start, end)
        } else {
            // this should not happen as we dropped all null values
            // rust will always use anyhow result by default
            return Result::<Option<polars::prelude::Series>, PolarsError>::Ok(None);

            // return Result::<Option<polars::prelude::Series>, polars::prelude::PolarsError>::Ok(Some(
            //     Series::new_empty("pos", &DataType::List((DataType::Int64).into())),
            // ));
        };
    // initialize variables for new features
    let mut out_list: Vec<Series> = Vec::with_capacity(start_series.len());

    // iter each feature
    // we sorted the group, so the most left feature is the first one
    for (id, (start, end)) in start_iter.zip(end_iter).enumerate() {
        let (curr_start, curr_end) = if let (Some(start), Some(end)) = (start, end) {
            (start, end)
        } else {
            // rust will always use anyhow result by default
            return Result::<Option<polars::prelude::Series>, polars::prelude::PolarsError>::Err(
                polars::prelude::PolarsError::ComputeError(
                    "Found missing value in the start or end column. This should not happen."
                        .into(),
                ),
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
            // start is sorted so we only need to check end
            if curr_end > window_end {
                window_end = curr_end;
            }
        } else {
            // if window_start >= window_end {
            out_list.push(Series::new(
                id.to_string().as_str(),
                [window_start, window_end],
            ));
            // }
            window_start = curr_start;
            window_end = curr_end;
        }
    }

    // Dont forget the last group
    out_list.push(Series::new("one more", [window_start, window_end]));

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

    // if we have only one feature, we return an empty list
    if start_series.len() == 1 {
        // return an empty list
        return Result::<Option<polars::prelude::Series>, PolarsError>::Ok(None);

        // return Result::<Option<polars::prelude::Series>, PolarsError>::Ok(Some(
        //     Series::new_empty("pos", &DataType::List((DataType::Int64).into())),
        // ));
    }

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
            return Result::<Option<polars::prelude::Series>, polars::prelude::PolarsError>::Err(
                polars::prelude::PolarsError::ComputeError(
                    "Found missing value in the start or end column. This should not happen."
                        .into(),
                ),
            );
        };

        out_list.push(Series::new(id.to_string().as_str(), [gap_start, gap_end]));
    }

    let ls = Series::new("pos", out_list).cast(&DataType::List((DataType::Int64).into()))?;
    Result::<Option<polars::prelude::Series>, PolarsError>::Ok(Some(ls))
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum ExtendOption {
    /// Extend the feature to the start
    Start,
    /// Extend the feature to the end
    End,
    /// Extend the feature to both sides
    Both,
}

#[derive(Clone, Copy)]
/// This enum is used for specifying the interval type of the input ranges.\
/// Each variant takes a `i64` to represent the coordinate system.\
/// For example, `Inclusive(1)` means 1-based [start,end]. This is also the format used in Grangers.\
pub enum IntervalType {
    /// Inclusive interval, e.g. [start, end]. GTF and GFF use this.
    Inclusive(i64),
    /// Exclusive interval, e.g. (start, end). VCF uses this.
    Exclusive(i64),
    /// Left inclusive, right exclusive, e.g. [start, end).
    LeftInclusive(i64),
    /// Left exclusive, right inclusive, e.g. (start, end]. BED uses this.
    RightInclusive(i64),
}

impl Default for IntervalType {
    fn default() -> Self {
        IntervalType::Inclusive(1)
    }
}

impl IntervalType {
    pub fn from<T: ToString>(file_type: T) -> Self {
        match file_type.to_string().to_lowercase().as_str() {
            "gtf" | "gff" | "bam" | "sam" => IntervalType::Inclusive(1),
            "bed" => IntervalType::RightInclusive(0),
            _ => panic!("The file type is not supported"),
        }
    }
    pub fn start_offset(&self) -> i64 {
        // 1 - c is for coordinate
        // the other 1 is for exclusive
        match self {
            IntervalType::Inclusive(c) => 1 - c,
            IntervalType::LeftInclusive(c) => 1 - c,
            IntervalType::RightInclusive(c) => 1 + 1 - c,
            IntervalType::Exclusive(c) => 1 + 1 - c,
        }
    }
    pub fn end_offset(&self) -> i64 {
        // 1 - c is for coordinate
        // -1 is for exclusive
        match self {
            IntervalType::Inclusive(c) => 1 - c,
            IntervalType::LeftInclusive(c) => -1 + 1 - c,
            IntervalType::RightInclusive(c) => 1 - c,
            IntervalType::Exclusive(c) => -1 + 1 - c,
        }
    }
}

#[cfg(test)]
mod tests {
    // use polars::prelude::*;
    use super::*;
    use crate::grangers::reader::gtf::{AttributeMode, Attributes, GStruct};

    use crate::grangers::grangers_utils::FileFormat;
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
        // let file_type = reader::FileFormat::GTF;

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
            attributes: Attributes::new(AttributeMode::Full, FileFormat::GTF).unwrap(),
            misc: Some(HashMap::new()),
        };
        let gsr = &mut gs;
        gsr.attributes.file_type = FileFormat::GTF;
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
        let _gr = Grangers::from_gstruct(gs, IntervalType::Inclusive(1)).unwrap();

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
        let fo = FlankOptions {
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
        let fo = FlankOptions {
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
        let fo = FlankOptions {
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
        let fo = FlankOptions {
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
        let fo = FlankOptions {
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
        let fo = FlankOptions {
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
            attributes: Attributes::new(AttributeMode::Full, FileFormat::GTF)?,
            misc: Some(HashMap::new()),
        };
        let gsr = &mut gs;
        gsr.attributes.file_type = FileFormat::GTF;
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

        let gr = Grangers::from_gstruct(gs, IntervalType::Inclusive(1))?;
        Ok(gr)
    }

    #[test]
    fn test_merge() {
        let df = df!(
            "seqnames" => ["chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr2"],
            "start" => [1i64, 5, 1, 11, 22, 1, 5],
            "end" => [10i64, 10, 10, 20, 30, 10, 30],
            "strand"=> ["+", "+", "-", "-", "-", "+", "-"],
            "gene_id" => ["g1", "g1", "g2", "g2", "g2", "g3", "g4"],
        )
        .unwrap();

        let gr = Grangers::new(df, None, None, None, IntervalType::Inclusive(1)).unwrap();
        // println!("gr: {:?}", gr.df());
        let mo = MergeOptions::new(vec!["seqnames", "gene_id"], false, 1).unwrap();

        // default setting
        let merged_gr: Grangers = gr.merge(&mo).unwrap();
        // println!("merged_gr: {:?}", merged_gr.df());
        let start: Vec<i64> = vec![1i64, 1, 22, 1, 5];
        let end: Vec<i64> = vec![10i64, 20, 30, 10, 30];
        assert_eq!(
            merged_gr
                .start()
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
            merged_gr
                .end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        // test ignore strand
        let mo = MergeOptions::new(vec!["seqnames"], true, 1).unwrap();

        let merged_gr = gr.merge(&mo).unwrap();
        // println!("merged_gr: {:?}", merged_gr.df());
        let start: Vec<i64> = vec![1i64, 22, 1];
        let end: Vec<i64> = vec![20i64, 30, 30];
        assert_eq!(
            merged_gr
                .start()
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
            merged_gr
                .end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        // slack=0
        // test ignore strand
        let mo = MergeOptions::new(vec!["seqnames", "gene_id"], false, 0).unwrap();

        let merged_gr: Grangers = gr.merge(&mo).unwrap();
        // println!("merged_gr: {:?}", merged_gr.df());
        let start: Vec<i64> = vec![1i64, 1, 11, 22, 1, 5];
        let end: Vec<i64> = vec![10i64, 10, 20, 30, 10, 30];
        assert_eq!(
            merged_gr
                .start()
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
            merged_gr
                .end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        // slack=2

        // test ignore strand
        let mo = MergeOptions::new(vec!["seqnames", "gene_id"], false, 2).unwrap();

        let merged_gr: Grangers = gr.merge(&mo).unwrap();
        // println!("merged_gr: {:?}", merged_gr.df());
        let start: Vec<i64> = vec![1i64, 1, 1, 5];
        let end: Vec<i64> = vec![10i64, 30, 10, 30];
        assert_eq!(
            merged_gr
                .start()
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
            merged_gr
                .end()
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

    #[test]
    fn test_gaps() {
        let df = df!(
            "seqnames" => ["chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr2"],
            "start" => [1i64, 12, 1, 5, 22, 1, 5],
            "end" => [10i64, 20, 10, 20, 30, 10, 30],
            "strand"=> ["+", "+", "+", "+", "+", "+", "-"],
            "gene_id" => ["g1", "g1", "g2", "g2", "g2", "g3", "g4"],
        )
        .unwrap();

        let gr = Grangers::new(df, None, None, None, IntervalType::Inclusive(1)).unwrap();
        // println!("gr: {:?}", gr.df());
        let mo = MergeOptions::new(vec!["seqnames", "gene_id"], false, 1).unwrap();

        // default setting
        let gr1: Grangers = gr.gaps(&mo).unwrap();
        // println!("gr1: {:?}", gr1.df());
        let start: Vec<i64> = vec![11i64, 21];
        let end: Vec<i64> = vec![11i64, 21];
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

    #[test]
    fn test_extend() {
        let say = false;

        let df = df!(
            "seqnames" => ["chr1", "chr1"],
            "start" => [1i64, 50],
            "end" => [10i64, 60],
            "strand"=> ["+", "-"],
            "gene_id" => ["g1", "g2"],
        )
        .unwrap();

        let gr = Grangers::new(df, None, None, None, IntervalType::Inclusive(1)).unwrap();

        if say {
            println!("gr: {:?}", gr.df());
        }

        // extend from start stranded
        let mut gr1 = gr.clone();
        gr1.extend(5, &ExtendOption::Start, false).unwrap();

        if say {
            println!("gr1: {:?}", gr1.df());
        }
        let start: Vec<i64> = vec![-4i64, 50];
        let end: Vec<i64> = vec![10i64, 65];
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

        // extend from start unstranded
        let mut gr1 = gr.clone();
        gr1.extend(5, &ExtendOption::Start, true).unwrap();

        if say {
            println!("gr1: {:?}", gr.df());
        }
        let start: Vec<i64> = vec![-4i64, 45];
        let end: Vec<i64> = vec![10i64, 60];
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

        // extend from end stranded
        let mut gr1 = gr.clone();
        gr1.extend(5, &ExtendOption::End, false).unwrap();

        if say {
            println!("gr1: {:?}", gr.df());
        }
        let start: Vec<i64> = vec![1i64, 45];
        let end: Vec<i64> = vec![15i64, 60];
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

        // extend from start unstranded
        let mut gr1 = gr.clone();
        gr1.extend(5, &ExtendOption::End, true).unwrap();

        if say {
            println!("gr1: {:?}", gr.df());
        }
        let start: Vec<i64> = vec![1i64, 50];
        let end: Vec<i64> = vec![15i64, 65];
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

        // extend from both
        let mut gr1 = gr.clone();
        gr1.extend(5, &ExtendOption::Both, true).unwrap();

        if say {
            println!("gr1: {:?}", gr.df());
        }
        let start: Vec<i64> = vec![-4i64, 45];
        let end: Vec<i64> = vec![15i64, 65];
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

    #[test]
    fn test_intorons() {
        let say = true;

        let df = df!(
            "seqnames" => ["chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1"],
            "feature_type" => ["gene", "transcript", "exon", "exon", "transcript", "exon", "exon"],
            "start" => [1i64, 1, 1, 71, 71, 71, 101],
            "end" => [200i64, 80, 20, 80, 150, 80, 150],
            "strand"=> ["+", "+", "+", "+", "+", "+", "+"],
            "gene_id" => ["g1", "g1", "g1", "g1", "g1", "g1", "g1"],
            "transcript_id" => [None, Some("t1"), Some("t1"), Some("t1"), Some("t2"), Some("t2"), Some("t2")],
        )
        .unwrap();

        let gr = Grangers::new(df, None, None, None, IntervalType::Inclusive(1)).unwrap();
        if say {
            println!("gr: {:?}", gr.df());
        }

        // extend from both
        let gr1 = gr.introns("gene_id").unwrap();
        if say {
            println!("gr1: {:?}", gr1.df());
        }
        let start: Vec<i64> = vec![21i64, 81];
        let end: Vec<i64> = vec![70i64, 100];
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

        // extend from both
        let gr1 = gr.introns("transcript_id").unwrap();
        if say {
            println!("gr1: {:?}", gr1.df());
        }
        let start: Vec<i64> = vec![21i64, 81];
        let end: Vec<i64> = vec![70i64, 100];
        let tid = vec![String::from("t1"), String::from("t2")];
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

        assert_eq!(
            gr1.column("transcript_id")
                .unwrap()
                .utf8()
                .unwrap()
                .into_iter()
                .map(|x| x.unwrap().to_string())
                .collect::<Vec<String>>(),
            tid
        );
    }
}
