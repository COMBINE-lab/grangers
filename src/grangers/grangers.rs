use crate::grangers::grangers_utils::*;
use crate::grangers::options::*;
use crate::grangers::reader;
use crate::grangers::reader::fasta::SeqInfo;
use anyhow::{bail, Context};
use noodles::fasta;
pub use noodles::fasta::record::Sequence;
use polars::lazy::dsl::col;
use polars::{lazy::prelude::*, prelude::*, series::Series};
use rust_lapper::{Interval, Lapper};
use std::collections::HashMap;
use std::fs;
use std::io::BufReader;
use std::ops::{Add, Mul, Sub};
use std::path::Path;
use tracing::{info, warn};

// const ESSENTIAL_FIELDS: [&str; 4] = ["seqname", "start", "end", "strand"];

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
    /// The name of the columns that are used to identify the genomic features
    field_columns: FieldColumns,
}

// IO
impl Grangers {
    /// add or replace a column in the dataframe
    // TODO: use this in the unstranded case
    pub fn add_column<T: SeriesTrait>(&mut self, series: Series) -> anyhow::Result<()> {
        self.df.with_column(series)?;
        Ok(())
    }

    /// check if the essential fields contain null values
    /// Each value in fields should either be a column name or a field of the FieldColumns struct
    fn any_nulls<T: AsRef<str>>(&self, fields: &[T], complain: bool) -> anyhow::Result<bool> {
        let mut valid = true;
        let df = self.df();
        let fc = self.field_columns();

        for col in fields {
            // get the name
            let col = if let Some(col) = fc.get_colname(col, false) {
                col
            } else if df.column(col.as_ref()).is_ok() {
                col.as_ref()
            } else {
                bail!(
                    "The column {} does not exist in the Grangers struct. Cannot check null values",
                    col.as_ref()
                )
            };

            if df.column(col)?.null_count() > 0 {
                valid = false;
                if complain {
                    warn!("The column {} contains null values. This will cause problems for most Grangers functions.", col);
                }
            }
        }
        if (!valid) & complain {
            warn!("You can drop null values by using the Grangers::drop_nulls() method.")
        }

        Ok(valid)
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
        field_columns: FieldColumns,
    ) -> anyhow::Result<Grangers> {
        let field_columns = if let Some(validated) = field_columns.validate(&df, true)? {
            validated
        } else {
            field_columns
        };

        // if the interval type is not inclusive, we need to convert it to inclusive
        if interval_type.start_offset() != 0 {
            df.with_column(df.column("start").unwrap() - interval_type.start_offset())?;
        }

        if interval_type.end_offset() != 0 {
            df.with_column(df.column("end").unwrap() - interval_type.end_offset())?;
        }

        // instantiate a new Grangers struct
        let gr = Grangers {
            df,
            misc,
            seqinfo,
            lapper,
            interval_type,
            field_columns,
        };

        // validate
        gr.any_nulls(&gr.df().get_column_names(), true)?;

        Ok(gr)
    }

    pub fn from_gstruct(
        gstruct: reader::GStruct,
        interval_type: IntervalType,
    ) -> anyhow::Result<Grangers> {
        // create dataframe!
        // we want to make some columns categorical because of this https://docs.rs/polars/latest/polars/docs/performance/index.html
        // fields
        let mut df_vec = vec![
            Series::new("seqname", gstruct.seqid),
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
        let gr = Grangers::new(
            df,
            None,
            gstruct.misc,
            None,
            interval_type,
            FieldColumns::default(),
        )?;
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
    /// get the reference of the field_columns
    pub fn field_columns(&self) -> &FieldColumns {
        &self.field_columns
    }

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

    /// get the reference to the seqname column
    pub fn seqid(&self) -> anyhow::Result<&Series> {
        self.column("seqname")
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
    /// get the intronic sequences of each gene or transcript or other custom groups.
    /// The `by` parameter should specify group used for identifying introns. The exons within each group will be merged before computing the introns.\
    /// The `exon_name` parameter is used for idenfitying exon records from the `feature_type` field column. If set as None, "exon" will be used.
    /// - This function requires that there is a valid feature_type field column for identifying exon records,
    /// and all exon records have a valid value for all fields defined in the grangers.field_columns.
    /// If succeed, it will return a Grangers struct with intronic ranges.
    pub fn introns(
        &self,
        by: IntronsBy,
        exon_name: Option<&str>,
        multithreaded: bool,
    ) -> anyhow::Result<Grangers> {
        // get exon records only
        // if this call succeeds, we can make sure that the exon records are all valid
        let exon_gr = self.exons(exon_name, multithreaded)?;
        let fc = exon_gr.field_columns();

        // parse `by`
        // we need to make sure `by` points to a valid column
        let by =  if self.is_column(by.as_ref()) {
            by.as_ref()
        }  else if let Some(by) = fc.get_colname(by.as_ref(), false) {
            by
        }else {
            bail!("The column {} neither exists in the dataframe nor represents a field of FieldColumns. Cannot use it to group exons", by.as_ref())
        };

        let mo = MergeOptions::new(&[by], false, 1)?;
        exon_gr.gaps(&mo)
    }

    /// filter exon records and deduplicate if needed according to the by parameter.\
    /// This function takes an optional `feature_name` value to identify exon records. If set as None, "exon" will be used. It is used for identifying exon records. This value should match exons' `feature_type` in the dataframe. \
    /// This function will not work if there are invalid exon records. The criteria are that all records marked as the defined `exon_name` in the defined `feature_type` column must have:
    /// - a valid seqname
    /// - a valid start
    /// - a valid end
    /// - a valid strand ("+" or "-")
    /// - a valid transcript_id
    /// - a valid gene_id
    /// - a valid exon_id
    pub fn exons(
        &self,
        feature_name: Option<&str>,
        multithreaded: bool,
    ) -> anyhow::Result<Grangers> {
        // validate the field_columns
        let mut fc = if let Some(fc) = self.field_columns().validate(self.df(), false)? {
            fc
        } else {
            self.field_columns().to_owned()
        };

        // parse feature_name
        let feature_name = if let Some(en) = feature_name {
            en.to_string()
        } else {
            "exon".to_string()
        };

        // feature_type can have null values, they will be ignored
        let feature_type = if let Some(ft) = &fc.feature_type {
            ft.as_str()
        } else {
            bail!("The {:?} column does not exist in the Grangers struct. This is required for identifying exon records. Cannot proceed", fc.feature_type);
        };

        if self.column(feature_type)?.null_count() > 0 {
            warn!("Found rows with a null `{}` value. These rows will be ignored when selecting exon records.", feature_type)
        }

        // polars way to subset
        let mut exon_df = self.df().filter(
            &self
                .df()
                .column(feature_type)?
                .iter()
                .map(|f| f.eq(&AnyValue::Utf8(feature_name.as_str())))
                .collect(),
        )?;

        // then we need to check if the exon records are valid
        check_col(&exon_df, Some(fc.seqname()))?;
        let start = check_col(&exon_df, Some(fc.start()))?;
        let end = check_col(&exon_df, Some(fc.end()))?;
        let strand = check_col(&exon_df, Some(fc.strand()))?;
        let transcript_id = check_col(&exon_df, fc.transcript_id())?;
        let _gene_id = check_col(&exon_df, fc.gene_id())?;
        let _exon_id = check_col(&exon_df, fc.exon_id())?;

        // make sure that stand is valid:
        // - the exons of each transcript are on the same strand
        // - the strand column does not contain other values than "+" or "-"
        let tx_strand = exon_df.select([transcript_id.as_str(), strand.as_str()])?
            .lazy()
            .groupby([transcript_id.as_str()])
            .agg([
                col(strand.as_str())
                    .unique()
                    .count()
                    .gt(lit(1)
                ).alias("is_solo")])
            .collect()?;
        if tx_strand.column("is_solo")?.bool()?.any() {
            bail!("The strand column contains null values. Cannot proceed.")
        }

        if !exon_df
            .column(strand.as_str())?
            .unique()?
            .is_in(&Series::new("valid strands", ["+", "-"]))?
            .all()
        {
            bail!("The strand column contains values that are not \"+\" or \"-\". Cannot proceed.")
        }

        // make sure start and end are positive
        if let Some(start_min) = exon_df.column(start.as_str())?.i64()?.min() {
            if start_min < 1 {
                bail!(
                    "The {} column contains values less than 1. Cannot proceed.",
                    start
                )
            }
        } else {
            bail!(
                "The Cannot get min value in the {} column. Cannot proceed.",
                start
            )
        }

        if let Some(end_min) = exon_df.column(end.as_str())?.i64()?.min() {
            if end_min < 1 {
                bail!(
                    "The {} column contains values less than 1. Cannot proceed.",
                    end.as_str()
                )
            }
        } else {
            bail!(
                "The Cannot get min value in the {} column. Cannot proceed.",
                end.as_str()
            )
        }
        
        // if there is an exon number field, it should be valid
        // if there is no exon number field, we will add one
        let exon_number = if let Some(exon_number) = fc.exon_number.clone() {
            if exon_df.column(exon_number.as_str())?.null_count() > 0 {
                bail!("The exon number column {} contains null values. Cannot proceed. You can delete this column so that the exons() method will compute the number according to the start and strand fields.", exon_number);
            }
            exon_number
        } else {
            // update exon number in fc
            fc.exon_number = Some("exon_number".to_string());
            // TODO: add an exon_number column if it doesn't exist
            exon_df = exon_df
                .lazy()
                .with_columns([
                    all(),
                    when(col(fc.strand()).first().eq(lit("+")))
                        .then(
                            col(fc.start())
                                .arg_sort(SortOptions {
                                    descending: false,
                                    nulls_last: false,
                                    multithreaded,
                                }).add(lit(1))
                        )
                        .otherwise(
                            col(fc.start())
                                .arg_sort(SortOptions {
                                    descending: true,
                                    nulls_last: false,
                                    multithreaded,
                                }).add(lit(1))
                        )
                        .over([transcript_id.as_str()])
                        .cast(DataType::Utf8)
                        .alias(fc.exon_number().unwrap()),
                    
                ])
                .collect()?;
            "exon_number".to_string()
        };

        // sort at the end according to exon_number
        // exon number is stored as a string, so we need to cast it to int
        exon_df = exon_df
            .lazy()
            .select([all()
                .sort_by(
                    [col(exon_number.as_str()).cast(DataType::UInt64)],
                    [false],
                )
                .over([col(transcript_id.as_str())])])
            .collect()?;

        // well done!
        let gr = Grangers::new(
            exon_df,
            self.seqinfo.clone(),
            self.misc.clone(),
            None,
            self.interval_type,
            fc,
        )?;

        Ok(gr)
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
                    .is_in(&Series::new("missing strand", ["."]))?
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
    pub fn flank(&self, width: i64, options: FlankOptions) -> anyhow::Result<Grangers> {
        let df = self
            .df()
            .clone()
            .lazy()
            .with_column(
                when(options.ignore_strand)
                    .then(lit(true))
                    .otherwise(col("strand").eq(lit("-")).neq(lit(options.start)))
                    .alias("start_flags"),
            )
            .with_column(
                // when both is true
                when(options.both)
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
                            .mul(when(lit(options.both)).then(lit(2)).otherwise(lit(1))))
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
            field_columns: self.field_columns.clone(),
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
        gr.df = gr.apply(&options.by, 0, options.ignore_strand, apply_gaps)?;
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

        let df = self.apply(
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
            field_columns: self.field_columns.clone(),
        })
    }

    fn apply<F>(
        &self,
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
        let df = self.df();
        let fc = self.field_columns();
        let seqname = check_col(df, Some(fc.seqname()))?;
        let start = check_col(df, Some(fc.start()))?;
        let end = check_col(df, Some(fc.end()))?;
        let strand = check_col(df, Some(fc.strand()))?;

        // we take the selected columns and add two more columns: start and end
        let mut selected = by.iter().map(|s| s.as_str()).collect::<Vec<&str>>();
        selected.append(vec![start.as_str(), end.as_str()].as_mut());

        // we want to sort the dataframe by first by columns, then the essential columns
        let mut sorted_by_exprs_essential = vec![col(seqname.as_str()), col(start.as_str()), col(end.as_str())];
        // we sort start in ascending order and end in descending order so that in each group, 
        let mut sorted_by_desc_essential = vec![false, false, true];
        if !ignore_strand {
            sorted_by_exprs_essential.push(col(strand.as_str()));
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
            warn!("Found null value(s) in the selected columns -- {:?}. As null will be used for grouping, we recommend dropping all null values by calling gr.drops_nulls() beforehand.", selected)
        }

        // we make sure that start and end have no null values
        if df.column("start")?.null_count() > 0 {
            bail!("The start column contains null values. Cannot proceed.")
        };


        // we will do the following
        // 1. sort the dataframe by the `by` columns + start and end columns
        // 2. group by the `by` columns
        // 3. build the new start and end columns by applying the `apply_fn` function
        // 4. explode the new start and end columns (makes the columns tall instead of wide)
        // 5. drop the intermediate columns
        df = df
            .lazy()
            // TODO: This can be replaced by select([all().sort(essentials).over(groups)]). Not sure if it is faster
            .sort_by_exprs(&sorted_by_exprs, &sorted_by_desc, false)
            .groupby(by.iter().map(|s| s.as_str()).collect::<Vec<&str>>())
            .agg([
                all().exclude([start.as_str(), end.as_str()]).first(),
                // process two columns at once
                // Notice the df is sorted
                as_struct(&[col(start.as_str()), col(end.as_str())])
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
                col("pos").arr().get(lit(0)).alias(start.as_str()),
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
                col("seqname"),
                col("start"),
                col("end"),
                col("strand"),
                all().exclude(["seqname", "start", "end", "strand"]),
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
            field_columns: self.field_columns.clone(),
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

// implement get sequence functions for Grangers
impl Grangers {
    /// Extract the transcript sequences in the Grangers object from the provided reference file.
    /// This function works only if the features are well defined:
    /// - Exon features cannot have a null "transcript_id".
    /// - The exons of a transcript should have the same "seqname" and "strand".
    /// - The exons of a transcript should not overlap with each other.
    /// - Each value in the "seqname" column should represent a sequence in the reference file.

    pub fn get_transcript_sequences<T: AsRef<Path>>(
        &self,
        fasta_path: T,
        exon_name: Option<&str>,
        multithreaded: bool,
    ) -> anyhow::Result<Vec<Sequence>> {
        // get exon_gr
        // exons() ensures that all exon records are valid,
        // and they have a valid exon number
        let exon_gr = self.exons(exon_name, multithreaded)?;
        let fc = exon_gr.field_columns();

        // all these fields are valid after exon()
        // TODO: check if ignoring extra fields helps memory usage
        // let selection = [
        //     fc.seqname(),
        //     fc.start(),
        //     fc.end(),
        //     fc.strand(),
        //     fc.gene_id().unwrap(),
        //     fc.transcript_id().unwrap(),
        //     fc.exon_id().unwrap(),
        //     fc.exon_number().unwrap(),
        // ];

        // Now, we read the fasta file and process each reference sequence at a time
        let reader = std::fs::File::open(fasta_path).map(BufReader::new)?;
        let mut reader = noodles::fasta::Reader::new(reader);
        // let mut seq_vec: Vec<Option<Sequence>> = vec![None; exon_gr.df().height()];
        let mut transcript_seq_vec: Vec<Sequence> = Vec::with_capacity(
            self.df()
                .column(fc.transcript_id().unwrap())?
                .unique()?
                .len(),
        );

        let mut chr_df: DataFrame;

        // we iterate the fasta reader. For each fasta reacord (usually chromosome), we do
        // 1. subset the dataframe by the seqname (chromosome name)
        // 2. for each gene, we get the sequence of all its exons
        // 3. for each transcript, we join the transcripts' exon sequences to get the sequence of the transcript
        for result in reader.records() {
            let record = result?;

            // filter the dataframe by the chromosome name
            // it is possible that the chromosome name
            let anyvalue_name =
                AnyValue::Utf8(record.name().strip_suffix(' ').unwrap_or(record.name()));

            // TODO: we can directly filter the gene name here.
            chr_df = exon_gr.df().filter(
                &self
                    .df()
                    .column(fc.seqname())?
                    .iter()
                    .map(|f| f.eq(&anyvalue_name))
                    .collect(),
            )?;

            // check if exons are in the range of the reference sequence
            if let Some(end_max) = chr_df.column("end")?.i64()?.max() {
                if end_max > record.sequence().len() as i64 {
                    bail!("Found exons that exceed the length of the reference sequence. Cannot proceed")
                }
            } else {
                bail!("Could not get the maximum end value of the exons. Cannot proceed")
            }

            // we get the sequence of a chromosome at a time
            let seq_vec = Grangers::get_sequences_fasta_record(&chr_df, &record, &OOBOption::Skip)?;
            if seq_vec
                .iter()
                .map(|f| f.is_none())
                .fold(0, |acc, x| acc + x as usize)
                > 0
            {
                bail!("Found invalid exons that exceed the chromosome length; Cannot proceed")
            }

            // we assemble the transcript sequences
            // exons() will sort the exons by the exon number
            // get_sequences_fasta_record() will take care of the strands
            // so here we just need to join the exon sequences
            // let mut transcript_seq_vec = vec![None; chr_df.height()];
            let mut tx_id_iter = chr_df
                .column("transcript_id")?
                .utf8()?
                .into_iter()
                .peekable();
            let mut curr_tx = if let Some(id) = tx_id_iter
                .peek()
                .with_context(|| "Could not get the first transcript id")?
            {
                id.to_string()
            } else {
                bail!("Could not get the first transcript id")
            };
            // This is the vector that stores the exon sequences of the current transcript
            // each element is a base, represented by its u8 value

            let mut exon_u8_vec: Vec<u8> = Vec::new();

            for (tx_id, seq) in tx_id_iter.zip(seq_vec.into_iter()) {
                if let (Some(tx_id), Some(seq)) = (tx_id, seq) {
                    // first we want to check if the transcript id is the same as the previous one
                    if tx_id == curr_tx {
                        // if it is the same, we extend the exon_vec with the current sequence
                        exon_u8_vec.extend(seq.as_ref().iter());
                    } else {
                        // // if it is not the same, we create a Sequence and push it to seq_vec
                        transcript_seq_vec.push(Sequence::from_iter(exon_u8_vec.clone()));
                        exon_u8_vec.clear();
                        // update the current transcript id
                        curr_tx = tx_id.to_string();
                    }
                } else {
                    bail!("Found null transcript id or empty exon sequence. This should not happen, please report this bug.")
                }
            }
        }
        Ok(transcript_seq_vec)
    }

    /// Extract the sequence of the features in the Grangers object from the provided reference file.
    /// Currently only fasta file is supported. This function four field columns: seqname, start, end, and strand.
    /// Arguments:
    /// - `genome_path`: the path to the reference genome file.
    /// - `file_format`: the format of the reference genome file. Currently only fasta is supported.
    /// - `oob_option`: the option for out-of-boundary positions. It can be either `Truncate` or `Skip`. If `Truncate`, the out-of-boundary positions will be truncated to the start or end of the sequence. If `Skip`, a None will be returned for features with OOB positions
    /// The function outputs the extracted sequence as a vector of `Option<Sequence>`. If the feature has OOB positions and the oob_option is set as `Skip`, the corresponding element in the vector will be None. The order of the vector follows the row order of the dataframe of the Grangers object.
    pub fn get_sequences_fasta<T: AsRef<Path>>(
        &self,
        fasta_path: T,
        ignore_strand: bool,
        oob_option: &OOBOption,
    ) -> anyhow::Result<Vec<Option<Sequence>>>
// anyhow::Result<Vec<fasta::record::Sequence>>
    {
        let fc = self.field_columns();
        // we need to map the sequence back to the original row order of the dataframe
        // So, we need to have a minimum copy of the dataset, which contains only the essential fields,
        // and add one more column representing the row order of the original dataframe
        let mut df = self
            .df
            .select([fc.seqname(), fc.start(), fc.end(), fc.strand()])?;
        df.with_column(Series::new(
            "col_order",
            (0..df.height() as u32).collect::<Vec<u32>>(),
        ))?;

        // if ignore strand, set the strand to +
        if ignore_strand {
            df.with_column(Series::new(fc.strand(), vec!["+"; df.height()]))?;
        }
        let mut chr_df: DataFrame;

        let reader = std::fs::File::open(fasta_path).map(BufReader::new)?;
        let mut reader = noodles::fasta::Reader::new(reader);
        let mut seq_vec: Vec<Option<Sequence>> = vec![None; df.height()];
        // we iterate the fasta reader. For each fasta reacord (usually chromosome), we do
        // 1. subset the dataframe by the chromosome name
        // 2. get the sequence of the features in the dataframe on that fasta record
        // 3. insert the sequence into the sequence vector according to the row order
        for result in reader.records() {
            let record = result?;

            // filter the dataframe by the chromosome name
            let anyvalue_name =
                AnyValue::Utf8(record.name().strip_suffix(' ').unwrap_or(record.name()));

            chr_df = df.filter(
                &self
                    .df()
                    .column(fc.seqname())?
                    .iter()
                    .map(|f| f.eq(&anyvalue_name))
                    .collect(),
            )?;

            // we get the sequence of a chromosome at a time
            let chr_seq_vec = Grangers::get_sequences_fasta_record(&chr_df, &record, oob_option)?;

            //
            for (idx, seq) in chr_df
                .column("col_order")?
                .u32()?
                .into_iter()
                .zip(chr_seq_vec.into_iter())
            {
                seq_vec[idx.unwrap() as usize] = seq;
            }
        }

        Ok(seq_vec)
    }

    /// Get the sequences of the intervals from one fasta record.
    /// The provided dataframe should be filtered by the reference name.
    ///  
    /// This function uses the fasta module from noodles.
    /// The fasta sequence struct is 1-based and inclusive, same as Grangers.
    /// To get [1,2,3,4,5] in rust, use 1..=5
    /// If the position less than 1 or exceeds the length of the sequence, it will return None.
    fn get_sequences_fasta_record(
        df: &DataFrame,
        record: &fasta::record::Record,
        oob_option: &OOBOption,
    ) -> anyhow::Result<Vec<Option<Sequence>>> {
        // initialize seq vector
        if df.column("seqname")?.unique()?.len() > 1 {
            bail!("The dataframe contains more than one reference name. Please filter the dataframe by the reference name first.")
        }

        let mut seq_vec = Vec::with_capacity(df.height());
        let ses = df.columns(["start", "end", "strand"])?;
        for ((start, end), strand) in ses[0]
            .i64()?
            .into_iter()
            .zip(ses[1].i64()?.into_iter())
            .zip(ses[2].utf8()?.into_iter())
        {
            if let (Some(start), Some(end)) = (start, end) {
                let (start, end) = if oob_option == &OOBOption::Truncate {
                    (
                        noodles::core::Position::try_from(std::cmp::max(1, start as usize))?,
                        noodles::core::Position::try_from(std::cmp::min(
                            record.sequence().len(),
                            end as usize,
                        ))?,
                    )
                } else {
                    (
                        noodles::core::Position::try_from(start as usize)?,
                        noodles::core::Position::try_from(end as usize)?,
                    )
                };
                let seq = record.sequence().slice(start..=end);

                if strand == Some("-") {
                    if let Some(seq) = seq {
                        seq_vec.push(Some(seq.complement().rev().collect::<Result<_, _>>()?))
                    }
                } else {
                    seq_vec.push(seq);
                };
            }
        }
        Ok(seq_vec)
    }
}

pub fn argsort1based<T: Ord>(data: &[T], descending: bool) -> Vec<usize> {
    let mut indices = (1..=data.len()).collect::<Vec<_>>();
    indices.sort_by_key(|&i| &data[i - 1]);
    if descending {
        indices.reverse();
    }
    indices
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
    // we sorted the group in the apply(), so the most left feature (and the widest one if there are many) is on the top
    let (mut window_start, mut window_end) =
        if let (Some(Some(start)), Some(Some(end))) = (start_iter.next(), end_iter.next()) {
            (start, end)
        } else {
            // this should not happen as we dropped all null values
            // rust will always use anyhow result by default
            return Result::<Option<polars::prelude::Series>, PolarsError>::Err(PolarsError::ComputeError("Found missing value in the start or end column. Cannot proceed.".into()));

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

#[cfg(test)]
mod tests {
    // use polars::prelude::*;
    use super::*;
    use crate::grangers::grangers_utils::*;
    use crate::grangers::options::*;
    use crate::grangers::reader::gtf::{AttributeMode, Attributes, GStruct};

    use crate::grangers::grangers_utils::FileFormat;

    const SAY: bool = true;
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

        if SAY {
            println!("gr: {:?}", gr.df());
        }

        // test flank with default parameters
        let fo = FlankOptions {
            start: true,
            both: false,
            ignore_strand: false,
        };
        let gr1 = gr.flank(10, fo).unwrap();
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

        let gr1 = gr.flank(-10, fo).unwrap();
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
        let gr1 = gr.flank(10, fo).unwrap();
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

        let gr1 = gr.flank(-10, fo).unwrap();
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
        let gr1 = gr.flank(10, fo).unwrap();
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

        let gr1 = gr.flank(-10, fo).unwrap();
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
        let gr1 = gr.flank(10, fo).unwrap();
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

        let gr1 = gr.flank(-10, fo).unwrap();
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
        let gr1 = gr.flank(10, fo).unwrap();
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

        let gr1 = gr.flank(-10, fo).unwrap();
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
        let gr1 = gr.flank(10, fo).unwrap();
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

        let gr1 = gr.flank(-10, fo).unwrap();
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
            "seqname" => ["chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr2"],
            "start" => [1i64, 5, 1, 11, 22, 1, 5],
            "end" => [10i64, 10, 10, 20, 30, 10, 30],
            "strand"=> ["+", "+", "-", "-", "-", "+", "-"],
            "gene_id" => ["g1", "g1", "g2", "g2", "g2", "g3", "g4"],
        )
        .unwrap();

        let gr = Grangers::new(
            df,
            None,
            None,
            None,
            IntervalType::Inclusive(1),
            FieldColumns::default(),
        )
        .unwrap();

        if SAY {
            println!("gr: {:?}", gr.df());
        }
        let mo = MergeOptions::new(&vec!["seqname", "gene_id"], false, 1).unwrap();

        // default setting
        let gr1: Grangers = gr.merge(&mo).unwrap();

        if SAY {
            println!("gr1: {:?}", gr1.df());
        }
        let start: Vec<i64> = vec![1i64, 1, 22, 1, 5];
        let end: Vec<i64> = vec![10i64, 20, 30, 10, 30];
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

        // test ignore strand
        let mo = MergeOptions::new(&vec!["seqname"], true, 1).unwrap();

        let gr1 = gr.merge(&mo).unwrap();

        if SAY {
            println!("gr1: {:?}", gr1.df());
        }
        let start: Vec<i64> = vec![1i64, 22, 1];
        let end: Vec<i64> = vec![20i64, 30, 30];
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

        // slack=0
        // test ignore strand
        let mo = MergeOptions::new(&vec!["seqname", "gene_id"], false, 0).unwrap();

        let gr1: Grangers = gr.merge(&mo).unwrap();

        if SAY {
            println!("gr1: {:?}", gr1.df());
        }
        let start: Vec<i64> = vec![1i64, 1, 11, 22, 1, 5];
        let end: Vec<i64> = vec![10i64, 10, 20, 30, 10, 30];
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

        // slack=2

        // test ignore strand
        let mo = MergeOptions::new(&vec!["seqname", "gene_id"], false, 2).unwrap();

        let gr1: Grangers = gr.merge(&mo).unwrap();

        if SAY {
            println!("gr1: {:?}", gr1.df());
        }
        let start: Vec<i64> = vec![1i64, 1, 1, 5];
        let end: Vec<i64> = vec![10i64, 30, 10, 30];
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
    fn test_gaps() {
        let df = df!(
            "seqname" => ["chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr2"],
            "start" => [1i64, 12, 1, 5, 22, 1, 5],
            "end" => [10i64, 20, 10, 20, 30, 10, 30],
            "strand"=> ["+", "+", "+", "+", "+", "+", "-"],
            "gene_id" => ["g1", "g1", "g2", "g2", "g2", "g3", "g4"],
        )
        .unwrap();

        let gr = Grangers::new(
            df,
            None,
            None,
            None,
            IntervalType::Inclusive(1),
            FieldColumns::default(),
        )
        .unwrap();
        if SAY {
            println!("gr: {:?}", gr.df());
        }

        let mo = MergeOptions::new(&vec!["seqname", "gene_id"], false, 1).unwrap();

        // default setting
        let gr1: Grangers = gr.gaps(&mo).unwrap();

        if SAY {
            println!("gr1: {:?}", gr1.df());
        }

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
        let df = df!(
            "seqname" => ["chr1", "chr1"],
            "start" => [1i64, 50],
            "end" => [10i64, 60],
            "strand"=> ["+", "-"],
            "gene_id" => ["g1", "g2"],
        )
        .unwrap();

        let gr = Grangers::new(
            df,
            None,
            None,
            None,
            IntervalType::Inclusive(1),
            FieldColumns::default(),
        )
        .unwrap();

        if SAY {
            println!("gr: {:?}", gr.df());
        }

        // extend from start stranded
        let mut gr1 = gr.clone();
        gr1.extend(5, &ExtendOption::Start, false).unwrap();

        if SAY {
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

        if SAY {
            println!("gr1: {:?}", gr1.df());
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

        if SAY {
            println!("gr1: {:?}", gr1.df());
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

        if SAY {
            println!("gr1: {:?}", gr1.df());
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

        if SAY {
            println!("gr1: {:?}", gr1.df());
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
    fn test_exons() {
        let df = df!(
            "seqname" => ["chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1"],
            "feature_type" => ["gene", "transcript", "exon", "exon", "transcript", "exon", "exon"],
            "start" => [1i64, 1, 1, 71, 71, 71, 101],
            "end" => [200i64, 80, 20, 80, 150, 80, 150],
            "strand"=> ["+", "+", "+", "+", "-", "-", "-"],
            "gene_id" => ["g1", "g1", "g1", "g1", "g1", "g1", "g1"],
            "transcript_id" => [None, Some("t1"), Some("t1"), Some("t1"), Some("t2"), Some("t2"), Some("t2")],
            "exon_id" => [None, None, Some("e1"), Some("e2"), None, Some("e3"), Some("e4")],
        )
        .unwrap();

        let gr = Grangers::new(
            df,
            None,
            None,
            None,
            IntervalType::Inclusive(1),
            FieldColumns::default(),
        )
        .unwrap();
        if SAY {
            println!("gr: {:?}", gr.df());
        }

        // extend from both
        let gr1 = gr.exons(None, true).unwrap();
        if SAY {
            println!("gr1: {:?}", gr1.df());
        }
        let start: Vec<i64> = vec![1i64, 71, 101, 71];
        let end: Vec<i64> = vec![20i64, 80, 150, 80];
        let exon_number = vec![1i64,2,1,2];
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
            gr1.column("exon_number")
                .unwrap()
                .cast(&DataType::Int64)
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            exon_number
        );

        let df = df!(
            "seqname" => ["chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1"],
            "feature_type" => ["gene", "transcript", "exon", "exon", "transcript", "exon", "exon"],
            "start" => [1i64, 1, 1, 71, 71, 71, 101],
            "end" => [200i64, 80, 20, 80, 150, 80, 150],
            "strand"=> ["+", "+", "+", "+", "-", "-", "-"],
            "gene_id" => ["g1", "g1", "g1", "g1", "g1", "g1", "g1"],
            "transcript_id" => [None, Some("t1"), Some("t1"), Some("t1"), Some("t2"), Some("t2"), Some("t2")],
            "exon_id" => [None, None, Some("e1"), Some("e2"), None, Some("e3"), Some("e4")],
            "exon_number" => [None, None, Some("1"), Some("2"), None, Some("1"), Some("2")],
        )
        .unwrap();

        let gr = Grangers::new(
            df,
            None,
            None,
            None,
            IntervalType::Inclusive(1),
            FieldColumns::default(),
        )
        .unwrap();
        if SAY {
            println!("gr: {:?}", gr.df());
        }

        // extend from both
        let gr1 = gr.exons(None, true).unwrap();
        if SAY {
            println!("gr1: {:?}", gr1.df());
        }
        let start: Vec<i64> = vec![1i64, 71, 71, 101];
        let end: Vec<i64> = vec![20i64, 80, 80, 150];
        let exon_number = vec![1i64,2,1,2];
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
            gr1.column("exon_number")
                .unwrap()
                .cast(&DataType::Int64)
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            exon_number
        );

    }
    #[test]
    fn test_intorons() {
        let df = df!(
            "seqname" => ["chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1"],
            "feature_type" => ["gene", "transcript", "exon", "exon", "transcript", "exon", "exon"],
            "start" => [1i64, 1, 1, 71, 71, 71, 101],
            "end" => [200i64, 80, 20, 80, 150, 80, 150],
            "strand"=> ["+", "+", "+", "+", "+", "+", "+"],
            "gene_id" => ["g1", "g1", "g1", "g1", "g1", "g1", "g1"],
            "transcript_id" => [None, Some("t1"), Some("t1"), Some("t1"), Some("t2"), Some("t2"), Some("t2")],
            "exon_id" => [None, None, Some("e1"), Some("e2"), None, Some("e3"), Some("e4")],
        )
        .unwrap();

        let gr = Grangers::new(
            df,
            None,
            None,
            None,
            IntervalType::Inclusive(1),
            FieldColumns::default(),
        )
        .unwrap();
        if SAY {
            println!("gr: {:?}", gr.df());
        }

        // extend from both
        let gr1 = gr.introns(IntronsBy::Gene, None, true).unwrap();
        if SAY {
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
        let gr1 = gr.introns(IntronsBy::Transcript, None, true).unwrap();
        if SAY {
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
