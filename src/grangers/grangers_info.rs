// TODO:
// 1. update write and get sequence functions to use the same implementation
use crate::grangers::grangers_utils::*;
use crate::grangers::options::*;
use crate::grangers::reader;
use crate::grangers::reader::fasta::SeqInfo;
use anyhow::{bail, Context};
use noodles::fasta;
pub use noodles::fasta::record::{Definition, Record, Sequence};
use polars::{lazy::prelude::*, prelude::*, series::Series};
use rust_lapper::{Interval, Lapper};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::BufReader;
use std::ops::{Add, Mul, Sub};
use std::path::Path;
use std::result::Result::Ok;
use tracing::debug;
use tracing::{info, warn};

// GTF files are 1-based with closed intervals.
/// The Grangers struct contains the following fields:
/// - df: the underlying polars dataframe
/// - misc: the additional information
/// - seqinfo: the reference information
/// - lapper: the lapper interval tree
/// - interval type: the interval type (inclusive/exclusive, 0-based/1-based)
/// - field_columns: the name of the columns that are used to identify the genomic features
///
/// **Notice** that Granges uses 1-based closed intervals for the ranges.
/// If your ranges are not like this, when instantiating new Grangers,
/// you should use the `interval_type` parameter to help the builder
/// to convert the ranges to 1-based closed intervals.
#[derive(Clone)]
pub struct Grangers {
    /// The underlying dataframe
    pub df: DataFrame,
    /// The additional information
    pub misc: Option<HashMap<String, Vec<String>>>,
    /// The reference information
    pub seqinfo: Option<SeqInfo>,
    /// The lapper interval tree
    pub lappers: Option<HashMap<[String;2], Lapper<u64, (usize, Vec<String>)>>>,
    lappers_ignore_strand: Option<bool>,
    /// The interval type
    pub interval_type: IntervalType,
    /// The name of the columns that are used to identify the genomic features
    pub field_columns: FieldColumns,
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
    pub fn any_nulls<T: AsRef<str>>(
        &self,
        fields: &[T],
        is_warn: bool,
        is_bail: bool,
    ) -> anyhow::Result<bool> {
        self.field_columns().is_valid(self.df(), true, true)?;
        let mut any_nulls = false;
        let df = self.df();

        for col in fields {
            if df
                .column(&self.get_column_name(col.as_ref(), false)?)?
                .null_count()
                > 0
            {
                any_nulls = true;
                if is_warn {
                    warn!("The dataframe contains null values in the given fields -- {:?}. This will cause problems for most Grangers functions.", col.as_ref());
                }
            }
        }

        df.drop_nulls(Some(&["start", "end"]))?;
        if (any_nulls) & is_bail {
            let fields_str = fields.iter().map(|s| s.as_ref()).collect::<Vec<_>>();
            bail!("The dataframe contains null values in the given fields -- {:?}. You can drop null values by calling `df.drop_nulls(Some(&{:?}))`", fields_str,fields_str);
        }

        if (any_nulls) & is_warn {
            let fields_str = fields.iter().map(|s| s.as_ref()).collect::<Vec<_>>();
            warn!(
                "You can drop null values by calling `df.drop_nulls(Some(&{:?}))`",
                fields_str
            )
        }

        Ok(any_nulls)
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
        // lappers: Option<HashMap<[String;2], Lapper<u64, (usize, Vec<String>)>>>,
        interval_type: IntervalType,
        mut field_columns: FieldColumns,
        verbose: bool,
    ) -> anyhow::Result<Grangers> {
        // we reject if the field_column is not valid
        if !field_columns.is_valid(&df, verbose, false)? {
            field_columns.fix(&df, false)?;
        }
        // if the interval type is not inclusive, we need to convert it to inclusive
        if interval_type.start_offset() != 0 {
            df.with_column(
                df.column(field_columns.start()).unwrap() - interval_type.start_offset(),
            )?;
        }

        if interval_type.end_offset() != 0 {
            df.with_column(df.column(field_columns.end()).unwrap() - interval_type.end_offset())?;
        }

        // instantiate a new Grangers struct
        let gr = Grangers {
            df,
            misc,
            seqinfo,
            lappers: None,
            lappers_ignore_strand: None,
            interval_type,
            field_columns,
        };

        // validate
        gr.any_nulls(&gr.field_columns().essential_fields(), verbose, true)?;

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
            Series::new("source", gstruct.source),
            Series::new("feature_type", gstruct.feature_type),
            Series::new("start", gstruct.start),
            Series::new("end", gstruct.end),
            Series::new("score", gstruct.score),
            Series::new("strand", gstruct.strand),
            Series::new("phase", gstruct.phase),
        ];

        //for essential attributes
        for (k, v) in gstruct.attributes.essential {
            if !v.is_empty() {
                df_vec.push(Series::new(k.as_str(), v));
            };
        }

        // for extra attributes
        if let Some(attributes) = gstruct.attributes.extra {
            for (k, v) in attributes {
                let s = if v.is_empty() {
                    Series::new_null(k.as_str(), gstruct.attributes.tally)
                } else {
                    Series::new(k.as_str(), v)
                };
                df_vec.push(s);
            }
        }
        let df = DataFrame::new(df_vec)?;
        let gr = Grangers::new(
            df,
            None,
            gstruct.misc,
            interval_type,
            FieldColumns::default(),
            true,
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

    pub fn get_gtf_df<T: AsRef<Path>>(&self, _file_path: T) -> anyhow::Result<DataFrame> {
        // get a copy of the dataframe
        let df = self.df();
        let mut fc = self.field_columns.clone();
        let mut attr_cols: HashSet<&str> = df.get_column_names().into_iter().collect();
        let mut existing_field_cols = Vec::new();
        let mut missing_field_cols = Vec::new();
        // let attr_cols = Vec::new();

        // we devide the columns into three groups
        // 1. existing fields: the fields that are already in the dataframe
        // 2. missing fields: the fields that are not in the dataframe
        // 3. attr_cols: attribute columns
        for field in GXFFIELDS {
            if let Ok(col) = self.get_column_name(field, false) {
                fc.update(field, col.as_str())?;
                attr_cols.remove(col.as_str());
                existing_field_cols.push(col);
            } else {
                fc.update(field, field)?;
                missing_field_cols.push(field);
            }
        }

        // then, the left elements in the df_cols are extra fields
        // we will concat the name with its value to make it as a valid GTF attribute column, and finally concat all attributes into a single column
        // key1 "value1"; key2 "value2"; key3 "value3"

        // for existing fields, we select them
        let mut expr_vec = existing_field_cols
            .iter()
            .map(|col_name| col(col_name.as_str()))
            .collect::<Vec<_>>();

        // for missing fields, we add each of them as a new column
        expr_vec.extend(
            missing_field_cols
                .iter()
                .map(|&col_name| lit(".").alias(col_name)),
        );

        // for attribute columns, we concat the name with its value
        expr_vec.extend(attr_cols.iter().map(|&col_name| {
            (when(col(col_name).is_not_null())
                .then(lit(col_name) + lit(" \"") + col(col_name).cast(DataType::Utf8) + lit("\";"))
                .otherwise(lit("")))
            .alias(col_name)
        }));

        // then, we prepare the final datafram for polar csv writer
        let out_df = self.df().clone()
            .lazy()
            .select(
                expr_vec
            )
            .select([
                col(fc.field("seqname").expect("Could not get the seqname field. Please report this issue in our GitHub repo.")),
                col(fc.field("source").expect("Could not get the source field. Please report this issue in our GitHub repo.")),
                col(fc.field("feature_type").expect("Could not get the feature_type field. Please report this issue in our GitHub repo.")),
                col(fc.field("start").expect("Could not get the start field. Please report this issue in our GitHub repo.")),
                col(fc.field("end").expect("Could not get the end field. Please report this issue in our GitHub repo.")),
                col(fc.field("score").expect("Could not get the score field. Please report this issue in our GitHub repo.")),
                col(fc.field("strand").expect("Could not get the strand field. Please report this issue in our GitHub repo.")),
                col(fc.field("phase").expect("Could not get the phase field. Please report this issue in our GitHub repo.")),
                concat_str(attr_cols.iter().map(|&c| col(c)).collect::<Vec<_>>(), "").alias("attributes"),
            ])
            .fill_nan(lit("."))
            .fill_null(lit("."))
            .collect()?;

        Ok(out_df)
    }

    pub fn write_gtf<T: AsRef<Path>>(&self, file_path: T) -> anyhow::Result<()> {
        let file_path = file_path.as_ref();

        // create the folder if it doesn't exist
        fs::create_dir_all(file_path.parent().with_context(|| {
            format!(
                "Could not get the parent directory of the given output file path {:?}",
                file_path.as_os_str()
            )
        })?)?;

        let mut out_df = self.get_gtf_df(file_path)?;

        let mut file = std::fs::File::create(file_path)?;
        CsvWriter::new(&mut file)
            .has_header(false)
            .with_delimiter(b'\t')
            .with_null_value(".".to_string())
            .finish(&mut out_df)?;

        Ok(())
    }
}

// get struct fields
impl Grangers {
    /// get the reference of the field_columns
    pub fn field_columns(&self) -> &FieldColumns {
        &self.field_columns
    }

    /// get the mutable reference of the field_columns
    pub fn field_columns_mut(&mut self) -> &FieldColumns {
        &mut self.field_columns
    }

    pub fn field_columns_checked(
        &self,
        is_warn: bool,
        is_bail: bool,
    ) -> anyhow::Result<&FieldColumns> {
        self.field_columns().is_valid(self.df(), is_warn, is_bail)?;
        Ok(self.field_columns())
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

    pub fn filter<T: AsRef<str>>(&self, by: T, values: &[T]) -> anyhow::Result<Grangers> {
        let column = self.get_column_name(by.as_ref(), false)?;
        let df = self
            .df()
            .filter(&self.df().column(&column)?.is_in(&Series::new(
                "values",
                values.iter().map(|s| s.as_ref()).collect::<Vec<&str>>(),
            ))?)?;
        Grangers::new(
            df,
            self.seqinfo().cloned(),
            self.misc.clone(),
            IntervalType::default(),
            self.field_columns().clone(),
            true,
        )
    }
}

impl Grangers {
    /// get the column name of the given name. The difference between this function and get_column_name is that this function returns a &str (not owned) to the column name, while the other one returns an owned string of the column name. If it is a field of the FieldColumns struct, return the corresponding column name; If it is a column name, return itself\
    pub fn get_column_name_str<T: AsRef<str>>(
        &self,
        name: T,
        bail_null: bool,
    ) -> anyhow::Result<&str> {
        let name = name.as_ref();
        // if it is a column name, return itself
        let fc = self.field_columns_checked(false, true)?;

        let name = if let Some(col) = fc.field(name) {
            col
        } else if self.df().get_column_names().contains(&name) {
            self.df().column(name)?.name()
        } else {
            bail!("{} is neither a column in the dataframe nor a field of FieldColumns. Cannot proceed", name)
        };

        if bail_null && self.df().column(name)?.null_count() > 0 {
            bail!("The column {} contains null values. Cannot proceed.", name)
        }

        Ok(name)
    }

    /// get the column name of the given name. If it is a field of the FieldColumns struct, return the corresponding column name; If it is a column name, return itself\
    pub fn get_column_name<T: AsRef<str>>(
        &self,
        name: T,
        bail_null: bool,
    ) -> anyhow::Result<String> {
        let name = name.as_ref();
        // if it is a column name, return itself
        let fc = self.field_columns_checked(false, true)?;

        let name = if let Some(col) = fc.field(name) {
            col
        } else if self.df().get_column_names().contains(&name) {
            self.df().column(name)?.name()
        } else {
            bail!("{} is neither a column in the dataframe nor a field of FieldColumns. Cannot proceed", name)
        };

        if bail_null && self.df().column(name)?.null_count() > 0 {
            bail!("The column {} contains null values. Cannot proceed.", name)
        }

        Ok(name.to_string())
    }

    pub fn column<T: AsRef<str>>(&self, name: T) -> anyhow::Result<&Series> {
        let direct = self.df().column(name.as_ref());

        // first check if it is a direct column
        // if not, try to get the column from the field_columns
        let col = if direct.is_ok() {
            direct?
        } else {
            // make sure that FieldColumns is valid
            let col = self.get_column_name_str(name.as_ref(), false)?;
            self.df().column(col)?
        };
        Ok(col)
    }

    pub fn columns<T: AsRef<str>>(&self, names: &[T]) -> anyhow::Result<Vec<&Series>> {
        let mut cols = Vec::new();
        for name in names {
            cols.push(self.column(name)?);
        }
        Ok(cols)
    }

    /// get the reference to the seqname column
    pub fn seqname(&self) -> anyhow::Result<&Series> {
        self.column(self.field_columns.seqname())
    }

    /// get the reference to the source column
    pub fn start(&self) -> anyhow::Result<&Series> {
        self.column(self.field_columns.start())
    }

    /// get the reference to the end column
    pub fn end(&self) -> anyhow::Result<&Series> {
        self.column(self.field_columns.end())
    }

    /// get the reference to the strand column
    pub fn strand(&self) -> anyhow::Result<&Series> {
        self.column(self.field_columns.strand())
    }

    /// get the reference to the score column
    pub fn score(&self) -> anyhow::Result<&Series> {
        self.column(
            self.field_columns
                .score()
                .with_context(|| "Could not get the score column from the dataframe.")?,
        )
    }

    /// get the reference to the phase (frame) column
    pub fn phase(&self) -> anyhow::Result<&Series> {
        self.column(
            self.field_columns
                .phase()
                .with_context(|| "Could not get the score column from the dataframe.")?,
        )
    }

    /// get the reference to the type column
    pub fn feature_type(&self) -> anyhow::Result<&Series> {
        self.column(
            self.field_columns
                .feature_type()
                .with_context(|| "Could not get the score column from the dataframe.")?,
        )
    }

    /// get the reference to the lappers of the Grangers struct
    // TODO
    // pub fn lappers(&self) -> &Option<Lapper<u64, Vec<String>>> {
    //     &self.lappers
    // }

    /// get the start, end, and strand columns as a dataframe
    pub fn range(&self) -> anyhow::Result<DataFrame> {
        let range = self.df.select([
            self.field_columns().start(),
            self.field_columns().end(),
            self.field_columns().strand(),
        ])?;
        Ok(range)
    }

    /// check if the GRanges object has a column of a given name
    pub fn is_column<T: AsRef<str>>(&self, name: T) -> bool {
        self.column(name).is_ok()
    }
}

// validate Grangers
impl Grangers {
    ///validate the Grangers struct
    /// - return error if the field_columns is invalid
    /// - return error if the dataframe contains null values in the essential fields
    /// - return warning if the dataframe contains null values in the additional fields
    pub fn validate(&self, is_warn: bool, is_bail: bool) -> anyhow::Result<bool> {
        if self.df().height() == 0 {
            if is_bail {
                bail!("The dataframe is empty. Cannot proceed.")
            } else {
                return Ok(false);
            }
        }
        // field columns will check if the fields exist in the dataframe
        let valid_fc = self.field_columns.is_valid(self.df(), false, true)?;

        if is_warn & !valid_fc {
            warn!("The field_columns is not valid. You can use Grangers::fix_field_columns() to fix it.");
            return Ok(false);
        }

        // Then we check null values in the essential fields
        let essential_nulls =
            self.any_nulls(&self.field_columns().essential_fields(), false, is_bail)?;

        if is_warn & essential_nulls {
            warn!("The dataframe contains null values in the essential fields - seqname, start, end and strand. You can use Grangers::drop_nulls() to drop them.");
            return Ok(false);
        }

        // then additional fields
        if is_warn {
            self.any_nulls(&self.df().get_column_names(), is_warn, false)?;
        }
        Ok(true)
    }

    /// call try fix the invalid fields in self.field_columns
    pub fn fix_field_columns(&mut self, is_warn: bool) -> anyhow::Result<()> {
        let mut field_columns = self.field_columns().clone();
        field_columns.fix(self.df(), is_warn)?;
        self.field_columns = field_columns;
        Ok(())
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
        by: Option<&str>,
        exon_feature: Option<&str>,
        keep_columns: Option<&[&str]>,
        multithreaded: bool,
    ) -> anyhow::Result<Grangers> {
        // get exon records only
        // if this call succeeds, we can make sure that the exon records are all valid
        let exon_gr = self.exons(exon_feature, multithreaded)?;
        let gene_id = self.field_columns().gene_id();
        let gene_name = self.field_columns().gene_name();
        let mut kc = Vec::new();
        if let Some(gene_id) = gene_id {
            kc.push(gene_id);
        }
        if let Some(gene_name) = gene_name {
            kc.push(gene_name);
        }

        if let Some(keep_columns) = keep_columns {
            kc.extend_from_slice(keep_columns);
        }

        let by_str = if let Some(by) = by {
            exon_gr.get_column_name_str(by, true)?
        } else {
            exon_gr.get_column_name_str("transcript_id", true)?
        };

        exon_gr.gaps(&[by_str], false, None, Some(&kc))
    }

    /// get the range of each gene. The range of each gene will be the union of the ranges of all exons of the gene.\
    /// Therefore, this function calls exons() internally to get the exon ranges.\
    /// To make this function work, the grangers must have well-defined exon records.\
    pub fn genes(
        &self,
        exon_feature: Option<&str>,
        multithreaded: bool,
    ) -> anyhow::Result<Grangers> {
        self.validate(false, true)?;
        let gene_id = self.get_column_name_str("gene_id", true)?;
        self.boundary(gene_id, exon_feature, multithreaded)
    }

    /// get the range of each transcript. The range of each gene will be the union of the ranges of all exons of the gene.\
    /// Therefore, this function calls exons() internally to get the exon ranges.\
    /// To make this function work, the grangers must have well-defined exon records.\
    pub fn transcripts(
        &self,
        exon_feature: Option<&str>,
        multithreaded: bool,
    ) -> anyhow::Result<Grangers> {
        let transcript_id = self.get_column_name_str("transcript_id", true)?;
        self.boundary(transcript_id, exon_feature, multithreaded)
    }

    /// get the range of each group in the given field column. The field column can be either a field of the FieldColumns struct, or a column in the Grangers's dataframe. \
    /// To make sense, one should provide the ID column of genes or transcripts. The range of each group will be the union of the ranges of all exons of the gene.\
    /// Therefore, this function calls exons() internally to get the exon ranges.\
    /// To make this function work, the grangers must have well-defined exon records.\
    pub fn boundary(
        &self,
        by: &str,
        exon_feature: Option<&str>,
        multithreaded: bool,
    ) -> anyhow::Result<Grangers> {
        self.validate(false, true)?;
        let mut exon_gr = self.exons(exon_feature, multithreaded)?;
        let fc = self.field_columns();
        let seqname = fc.seqname();
        let start = fc.start();
        let end = fc.end();
        let strand = fc.strand();
        let by = self.get_column_name_str(by, true)?;

        // check if genes are well defined: all features of a gene should have a valid seqname, start, end, and strand
        let any_invalid = exon_gr
            .df()
            .select([seqname, strand, by])?
            .lazy()
            .groupby([by])
            .agg([
                col(seqname)
                    .unique()
                    .count()
                    .neq(lit(1))
                    .alias("seqname_any"),
                col(strand).unique().count().neq(lit(1)).alias("strand_any"),
            ])
            .select([col("seqname_any").any(), col("strand_any").any()])
            .collect()?
            .get_row(0)?
            .0
            .into_iter()
            .any(|c| c != AnyValue::Boolean(false));

        if any_invalid {
            bail!("The genes are not well defined. All features of a gene should be defined in the same seqname and strand. Cannot proceed.")
        };

        exon_gr.df = exon_gr
            .df
            .lazy()
            .groupby([seqname, by, strand])
            .agg([col(start).min(), col(end).max()])
            .collect()?;

        exon_gr.fix_field_columns(false)?;

        Ok(exon_gr)
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
        exon_feature: Option<&str>,
        multithreaded: bool,
    ) -> anyhow::Result<Grangers> {
        // validate itself
        // essential fields have no null values
        // all fields correspond to a column in the dataframe
        self.validate(false, true)?;

        let exon_feature = if let Some(exon_feature) = exon_feature {
            exon_feature
        } else {
            "exon"
        };

        // feature_type can have null values, they will be ignored
        let feature_type = self.get_column_name_str("feature_type", false)?;

        if self.column(feature_type)?.null_count() > 0 {
            warn!("Found rows with a null `{}` value. These rows will be ignored when selecting exon records.", feature_type)
        }

        // polars way to subset
        let mut exon_gr = self.filter(feature_type, &[exon_feature])?;

        // We know fields are valid, then we need to check nulls
        let mut fc = self.field_columns().clone();
        let seqname_s = self.get_column_name("seqname", true)?;
        let seqname = seqname_s.as_str();
        let start_s = self.get_column_name("start", true)?;
        let start = start_s.as_str();
        let end_s = self.get_column_name("end", true)?;
        let end = end_s.as_str();
        let strand_s = self.get_column_name("strand", true)?;
        let strand = strand_s.as_str();
        let transcript_id_s = self.get_column_name("transcript_id", false)?;
        let transcript_id = transcript_id_s.as_str();

        // make sure that transcript_id is not null
        if exon_gr.column(transcript_id)?.null_count() > 0 {
            bail!("Found exon features with a null transcript_id; Cannot proceed")
        }

        // make sure that strand is valid
        if !exon_gr
            .column(strand)?
            .unique()?
            .is_in(&Series::new("valid strands", ["+", "-"]))?
            .all()
        {
            bail!("Found exons that do not have a valid strand (+ or -). Cannot proceed.")
        }

        // make sure that stand is valid:
        // - the exons of each transcript are on the same strand
        // - the strand column does not contain other values than "+" or "-"
        let tx_strand = exon_gr
            .df()
            .select([seqname, transcript_id, strand])?
            .lazy()
            .groupby([seqname, transcript_id])
            .agg([col(strand).unique().count().gt(lit(1)).alias("is_solo")])
            .collect()?;
        if tx_strand.column("is_solo")?.bool()?.any() {
            bail!(
                "Found transcripts with exons from multiple chromosomes or strands; Cannot proceed"
            )
        }

        // make sure start and end are positive
        if let Some(start_min) = exon_gr.column(start)?.i64()?.min() {
            if start_min < 1 {
                bail!("Found exons with non-positive start position. Cannot proceed.")
            }
        } else {
            bail!(
                "Cannot get min value in the {} column. Cannot proceed.",
                start
            )
        }

        if let Some(end_min) = exon_gr.column(end)?.i64()?.min() {
            if end_min < 1 {
                bail!("Found exons with non-positive start position. Cannot proceed.")
            }
        } else {
            bail!(
                "Cannot get min value in the {} column. Cannot proceed.",
                end
            )
        }

        // if there is an exon number field, it should contain no null values
        // otherwise we will compute the exon number from exon start position
        if fc.exon_number.is_some() && exon_gr.column(fc.exon_number().unwrap())?.null_count() > 0 {
            warn!("The {} column contains null values. Will compute the exon number from exon start position .", fc.exon_number().unwrap());
            fc.exon_number = None;
        }

        // if there is no exon number field, we will compute the exon number from exon start position
        let exon_number = if let Some(exon_number) = fc.exon_number() {
            exon_number.to_string()
        } else {
            // update exon number in fc
            fc.exon_number = Some("exon_number".to_string());

            exon_gr.add_order(
                Some(&[transcript_id]),
                "exon_number",
                Some(1),
                multithreaded,
            )?;

            "exon_number".to_string()
        };

        // sort at the end according to exon_number
        // exon number is stored as a string, so we need to cast it to int
        exon_gr.df = exon_gr
            .df
            .lazy()
            .with_column(col(exon_number.as_str()).cast(DataType::UInt32))
            .select([all().sort_by(
                [
                    col(seqname).cast(DataType::Categorical(None)),
                    col(strand).cast(DataType::Categorical(None)),
                    col(transcript_id).cast(DataType::Categorical(None)),
                    col(exon_number.as_str()),
                ],
                [false],
            )])
            .collect()?;

        // well done!
        fc.fix(exon_gr.df(), false)?;
        exon_gr.field_columns = fc;
        Ok(exon_gr)
    }

    /// extend each genomic feature by a given length from the start, end, or both sides.
    pub fn extend(
        &mut self,
        length: i64,
        extend_option: &ExtendOption,
        ignore_strand: bool,
    ) -> anyhow::Result<()> {
        self.validate(false, true)?;
        let start_s = self.get_column_name("start", true)?;
        let start = start_s.as_str();
        let end_s = self.get_column_name("end", true)?;
        let end = end_s.as_str();
        let strand_s = self.get_column_name("strand", true)?;
        let strand = strand_s.as_str();

        // if contains null value in strand, we cannot do strand-specific extension
        if (!ignore_strand) & (extend_option != &ExtendOption::Both)
            && self.column(strand)?.is_null().any()
                | (!self
                    .column(strand)?
                    .unique()?
                    .is_in(&Series::new("valid stands", VALIDSTRANDS))?
                    .all())
        {
            bail!("The strand column contains values other than {:?}. Please remove them first or set ignore_strand to true.", VALIDSTRANDS)
        }

        // if both, then we extend both sides
        if let ExtendOption::Both = extend_option {
            self.df
                .with_column(self.df.column(start)?.clone() - length)?;
            self.df.with_column(self.df.column(end)?.clone() + length)?;
            return Ok(());
        }

        if ignore_strand {
            match extend_option {
                ExtendOption::Start => {
                    self.df
                        .with_column(self.df.column(start)?.clone() - length)?;
                    return Ok(());
                }
                ExtendOption::End => {
                    self.df.with_column(self.df.column(end)?.clone() + length)?;
                    return Ok(());
                }
                _ => {}
            }
        } else {
            let mut df = self.df().select([start, end, strand])?;
            df = df
                .lazy()
                .with_columns([
                    // we first consider the start site
                    // when the strand is + and extend from start, or the strand is - and extend from end, we extend the start site
                    when(
                        col(strand)
                            .eq(lit("+"))
                            .eq(lit(extend_option == &ExtendOption::Start))
                            .or(col(strand)
                                .eq(lit("-"))
                                .eq(lit(extend_option == &ExtendOption::End))),
                    )
                    .then(col(start).sub(lit(length)))
                    .otherwise(col(start))
                    .alias(start),
                    // then the end site
                    // when the strand is - and extend from start, or the strand is + and extend from end, we extend the end site
                    when(
                        col(strand)
                            .eq(lit("-"))
                            .eq(lit(extend_option == &ExtendOption::Start))
                            .or(col(strand)
                                .eq(lit("+"))
                                .eq(lit(extend_option == &ExtendOption::End))),
                    )
                    .then(col(end).add(lit(length)))
                    .otherwise(col(end))
                    .alias(end),
                ])
                .collect()?;

            // Then we update df
            self.df.with_column(df.column(start)?.clone())?;
            self.df.with_column(df.column(end)?.clone())?;
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
        self.validate(false, true)?;
        let start_s = self.get_column_name("start", true)?;
        let start = start_s.as_str();
        let end_s = self.get_column_name("end", true)?;
        let end = end_s.as_str();
        let strand_s = self.get_column_name("strand", true)?;
        let strand = strand_s.as_str();

        let df = self
            .df()
            .clone()
            .lazy()
            .with_column(
                when(options.ignore_strand)
                    .then(lit(true))
                    .otherwise(col(strand).eq(lit("-")).neq(lit(options.start)))
                    .alias("start_flags_temp"),
            )
            .with_column(
                // when both is true
                when(options.both)
                    .then(
                        // when start_flag is true
                        when(col("start_flags_temp").eq(lit(true)))
                            .then(col(start) - lit(width).abs())
                            // when start_flag is false
                            .otherwise(col(end) - lit(width).abs() + lit(1)),
                    )
                    // when both is false
                    .otherwise(
                        // if width >= 0:
                        when(width >= 0)
                            .then(
                                // tstart = all_starts[idx] - abs(width) if sf else all_ends[idx] + 1
                                when(col("start_flags_temp").eq(lit(true)))
                                    .then(col(start) - lit(width))
                                    .otherwise(col(end) + lit(1)),
                            )
                            .otherwise(
                                // tstart = all_starts[idx] if sf else all_ends[idx] + abs(width) + 1
                                when(col("start_flags_temp").eq(lit(true)))
                                    .then(col(start))
                                    .otherwise(col(end) + lit(width) + lit(1)),
                            ),
                    )
                    .alias(start),
            )
            .select([
                // everything except end and start_flags
                all().exclude([end, "start_flags_temp"]),
                // new_ends.append(tstart + (width * (2 if both else 1) - 1))
                col(start)
                    .add(
                        (lit(width)
                            .abs()
                            .mul(when(lit(options.both)).then(lit(2)).otherwise(lit(1))))
                        .sub(lit(1)),
                    )
                    .alias(end),
            ])
            .select(
                self.df()
                    .get_column_names()
                    .iter()
                    .map(|x| col(x))
                    .collect::<Vec<Expr>>(),
            )
            .collect()?;
        Grangers::new(df, self.seqinfo.clone(), self.misc.clone(), self.interval_type, self.field_columns.clone(), false)
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
    pub fn gaps<T: AsRef<str>>(
        &self,
        by: &[T],
        ignore_strand: bool,
        slack: Option<usize>,
        keep_columns: Option<&[&str]>,
    ) -> anyhow::Result<Grangers> {
        // merge returns a sorted and merged Grangers object
        let mut gr = self.merge(by, ignore_strand, slack, keep_columns)?;

        gr.df = gr.apply(by, None, ignore_strand, apply_gaps, keep_columns)?;
        Ok(gr)
    }

    pub fn add_order(
        &mut self,
        by: Option<&[&str]>,
        name: &str,
        offset: Option<u32>,
        multithreaded: bool,
    ) -> anyhow::Result<()> {
        self.validate(false, true)?;
        if let Some(by) = by {
            let mut by_col = Vec::new();
            for b in by.iter() {
                by_col.push(col(self.get_column_name_str(b, true)?))
            }

            let strand_s = self.get_column_name("strand", false)?;
            let strand = strand_s.as_str();
            let start_s = self.get_column_name("start", true)?;
            let start = start_s.as_str();

            self.df = self
                .df()
                .clone()
                .lazy()
                .with_column(
                    when(col(strand).first().eq(lit("+")))
                        .then(
                            col(start)
                                .arg_sort(SortOptions {
                                    descending: false,
                                    nulls_last: false,
                                    multithreaded,
                                })
                                .add(lit(1)),
                        )
                        .otherwise(
                            col(start)
                                .arg_sort(SortOptions {
                                    descending: true,
                                    nulls_last: false,
                                    multithreaded,
                                })
                                .add(lit(1)),
                        )
                        .over(by)
                        .cast(DataType::Utf8)
                        .alias(name),
                )
                .collect()?;
        } else {
            self.df = self.df.with_row_count(name, offset)?;
        }
        Ok(())
    }

    /// drop rows inplace that contain missing values in the given columns.
    /// if `cols` is None, all columns will be checked.
    pub fn drop_nulls(&mut self, fields: Option<&[&str]>) -> anyhow::Result<()> {
        self.validate(false, true)?;
        // check the validity of the column names
        let cols = match fields {
            Some(names) => self.columns(names)?.iter().map(|s| s.name()).collect(),
            None => self.df.get_column_names(),
        };

        *self.df_mut() = self.df().drop_nulls(Some(&cols))?;
        Ok(())
    }

    /// merge the features by the given columns via the `by` argument to generate a new Grangers object.
    /// The `by` columns cannot have any missing value. If yours' do, run something like `gr.drop_nulls(&vec!["by_col1".to_string(), "by_col2".to_string()])` first.
    /// ### Argument:
    /// Merge Option: a struct containing the following fields:
    /// - `by`: a vector of string representing which group(s) to merge by. Each string should be a valid column name.
    /// To make sense, one should provide the ID column of genes or transcripts. The range of each group will be the union of the ranges of all features in each group.\
    /// - `ignore_strand`: whether to ignore the strand information when merging.
    /// - `slack`: the maximum distance between two features to be merged. Slack should be a non-negative integer. For example we have three intervals, [1,2], [3,4] and [4,5]
    ///     - slack = 1 means [1,2] and [3,4] will be merged as [1,4]. This is the default (and desired) behavior unless you do not want to merge adjacent intervals.
    ///     - slack = 0 means [1,2] and [2,3] will stay separated although they are adjacent.
    ///     - slack=2 means[1,2] and [4,5] will be merged though there is a base [3,3] in between that separate them. The slack here is equivalent to the `min.gapwidth` argument in the GenomicRanges::reduce() function in R.
    pub fn merge<T: AsRef<str>>(
        &self,
        by: &[T],
        ignore_strand: bool,
        slack: Option<usize>,
        keep_columns: Option<&[&str]>,
    ) -> anyhow::Result<Grangers> {
        self.validate(false, true)?;
        let df = self.apply(by, slack, ignore_strand, apply_merge, keep_columns)?;

        Grangers::new(
            df,
            self.seqinfo.clone(),
            self.misc.clone(),
            IntervalType::default(),
            self.field_columns.clone(),
            false,
        )
    }

    fn apply<F, T: AsRef<str>>(
        &self,
        by: &[T],
        slack: Option<usize>,
        ignore_strand: bool,
        apply_fn: F,
        keep_columns: Option<&[&str]>,
    ) -> anyhow::Result<DataFrame>
    where
        F: Fn(Series, i64) -> Result<Option<polars::prelude::Series>, PolarsError>
            + Copy
            + std::marker::Send
            + std::marker::Sync
            + 'static,
    {
        self.validate(false, true)?;

        // these are all valid after validateion
        let df = self.df();
        let fc = self.field_columns();
        let seqname = fc.seqname();
        let start = fc.start();
        let end = fc.end();
        let strand = fc.strand();

        // this makes sure that field_column fields and their corresponding columns appear only once
        let mut by_hash: HashSet<&str> = HashSet::with_capacity(by.len());
        for name in by.iter() {
            let name = self.get_column_name_str(name.as_ref(), true)?;
            by_hash.insert(name);
        }

        // make sure that essential fields are not in the by hash

        if by_hash.take(start).is_some() | by_hash.take(end).is_some() {
            bail!("The provided `by` vector cannot contain the start or end column")
        };

        let slack = if let Some(s) = slack {
            if s < 1 {
                warn!("It usually does not make sense to set slack as zero.")
            }
            s as i64
        } else {
            1
        };

        if ignore_strand {
            if by_hash.take(strand).is_some() {
                warn!("Remove `strand` from the provided `by` vector as the ignored_strand flag is set.")
            }
        } else {
            by_hash.insert(strand);
        }

        // add chromosome name and strand if needed
        if by_hash.insert(seqname) {
            debug!("Added `seqname` to the `by` vector as it is required.")
        };
        let by: Vec<&str> = by_hash.into_iter().collect();

        // we take the selected columns and add two more columns: start and end
        let mut selected = by.to_vec();
        if !selected.contains(&seqname) {
            selected.push(seqname);
        }
        if !selected.contains(&start) {
            selected.push(start);
        }
        if !selected.contains(&end) {
            selected.push(end);
        }
        if !ignore_strand && !selected.contains(&strand) {
            selected.push(strand);
        }

        // let's see polars' way of checking missing values saying df.isna().sum()
        if self.any_nulls(&selected, true, false)? {
            warn!("Found null value(s) in the selected columns -- {:?}. As null will be used for grouping, we recommend dropping all null values by calling gr.drops_nulls() beforehand.", selected)
        }

        // we want to sort the dataframe by first by columns, then the essential columns
        let mut sorted_by_exprs_essential = vec![col(seqname), col(start), col(end)];
        // we sort start in ascending order and end in descending order so that in each group,
        let mut sorted_by_desc_essential = vec![false, false, true];
        if !ignore_strand {
            sorted_by_exprs_essential.push(col(strand));
            sorted_by_desc_essential.push(false);
        }

        let mut sorted_by_exprs: Vec<Expr> = by
            .iter()
            .filter(|&n| !sorted_by_exprs_essential.contains(&col(n)))
            .map(|n| col(n))
            .collect();

        let mut sorted_by_desc = vec![false; sorted_by_exprs.len()];
        sorted_by_exprs.extend(sorted_by_exprs_essential);
        sorted_by_desc.extend(sorted_by_desc_essential);

        // the lazy API of polars takes the ownership of a dataframe
        // we want to keep the keep_columns
        if let Some(keep_columns) = keep_columns {
            for c in keep_columns {
                if !selected.contains(&self.get_column_name_str(c, false)?) {
                    selected.push(c);
                }
            }
        }

        let mut df = df.select(&selected)?;

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
            .groupby(by.iter().map(|s| col(s)).collect::<Vec<Expr>>())
            .agg([
                all().exclude([start, end]).first(),
                // process two columns at once
                // Notice the df is sorted
                as_struct(&[col(start), col(end)])
                    .apply(
                        move |s| apply_fn(s, slack),
                        GetOutput::from_type(DataType::List((DataType::Int64).into())),
                    )
                    .alias("start_end_list-temp-nobody-will-use-this-name-right"),
            ])
            .explode(["start_end_list-temp-nobody-will-use-this-name-right"])
            // with_columns returns all columns and adds extra
            // as we can't drop a non-existing column, we need to add a dummy column
            .with_columns([
                col("start_end_list-temp-nobody-will-use-this-name-right")
                    .list()
                    .get(lit(0))
                    .alias(start),
                col("start_end_list-temp-nobody-will-use-this-name-right")
                    .list()
                    .get(lit(1))
                    .alias(end),
                lit(".").alias(if ignore_strand {
                    strand
                } else {
                    "ignore_strand-temp-nobody-will-use-this-name-right"
                }),
            ])
            .drop_nulls(Some(vec![cols([start, end])]))
            .with_column(
                lit(".")
                    .cast(DataType::Utf8)
                    .alias("ignore_strand-temp-nobody-will-use-this-name-right"),
            )
            .select([all().exclude([
                "start_end_list-temp-nobody-will-use-this-name-right",
                "ignore_strand-temp-nobody-will-use-this-name-right",
            ])])
            // rearrange the columns
            .select([
                col(seqname),
                col(start),
                col(end),
                col(strand),
                all().exclude([seqname, start, end, strand]),
            ])
            // groupby is multithreaded, so the order do not preserve
            .sort_by_exprs(sorted_by_exprs, sorted_by_desc, false)
            .collect()?;

        Ok(df)
    }
}

// lappers
impl Grangers {

    /// merge the features by the given columns via the `by` argument to generate a new Grangers object.
    /// *** Argument:
    /// - `by`: a vector of string representing which group(s) to merge by. Each string should be a valid column name.
    /// - `ignore_strand`: whether to ignore the strand information when merging.
    /// - `slack`: the maximum distance between two features to be merged.
    // TODO: add slack to this function
    // pub fn lapper_merge(&mut self, by: Vec<String>, _slack: i64) -> anyhow::Result<Grangers> {
    //     // self.build_lapper(&by[..])?;

    //     let mut lappers = if let Some(lappers) = self.lappers.clone() {
    //         lappers
    //     } else {
    //         bail!("Could not find the lappers that was just built; Please report this bug!")
    //     };

    //     lappers.merge_overlaps();

    //     if !lappers.overlaps_merged {
    //         info!("Did not find overlapping features. Nothing to merge.")
    //     }

    //     let lapper_vecs = LapperVecs::new(&lappers, &by);
    //     let df = lapper_vecs.into_df()?;
    //     Ok(Grangers {
    //         df,
    //         seqinfo: self.seqinfo.clone(),
    //         misc: self.misc.clone(),
    //         lappers: Some(lappers),
    //         interval_type: self.interval_type,
    //         field_columns: self.field_columns.clone(),
    //     })
    // }

    /// Build a Lapper data structure for each `seqname` or (`seqname`,`strand`) for interval search and overlap detection. Lappers only work with non-negative start and end positions.
    /// This function takes three parameters:
    /// 1. `ignore_strand`: If true, build a lapper for each unique `seqname`. Otherwise, build a lapper for each unique (`seqname`,`strand`) pair.
    /// 2. `ignore_invalid`: If true, ignore invalid features. Otherwise, return error if found invalid features. A feature is valid if it has a positive start and end site, a valid `seqname`, and a valid strand if `ignore_strand` is not set.
    /// 3. `overwrite`: If true, overwrite if the lapper field exist. Otherwise, return without doing anything.
    /// **Note** that as this functionality depends on an external crate assuming right-exclusive intervals, i.e., [start,end), you will see that the end of each lapper interval is 1 base greater than its corresponding features in the Grangers object, which relies on inclusive intervals, i.e., [start,end]. Therefore, if you want to use the lapper outise of Grangers function, please be aware of the difference.
    /// Building the lapper for a Grangers object will provide you the following functions (from [rust-lapper's doc](https://docs.rs/rust-lapper/latest/rust_lapper/struct.Lapper.html#method.count)):
    /// 1. find: Find all intervals in the lapper that overlap a given interval
    /// 2. seek: Seek all intervals in the lapper that overlap a given interval. It uses a linear search from the last query instead of a binary search. A reference to a cursor must be passed in. This reference will be modified and should be reused in the next query. This allows seek to not need to make the lapper object mutable, and thus use the same lapper accross threads.
    /// 3. count: Count all intervals in the lapper that overlap a given interval. This performs two binary search in order to find all the excluded elements, and then deduces the intersection from there. See BITS for more details.
    /// 4. cov: Get the number of positions covered by the intervals in Lapper. This provides immutable access if it has already been set, or on the fly calculation.

    pub fn build_lapper<T: AsRef<str>>(&mut self, ignore_invalid: bool, ignore_strand: bool, overwrite: bool, group_by: &[T]) -> anyhow::Result<()> {
        // rust-lapper
        let start_time = std::time::Instant::now();
        // validate the Grangers object
        self.validate(false, true)?;

        // first we want to check if the lappers have already been built
        if self.lappers.is_some() && !overwrite {
            info!("The lappers has already been built. Nothing to do.");
            return Ok(());
        }

        let mut by = Vec::new();
        for b in group_by.iter() {
            by.push(self.get_column_name_str(b.as_ref(), true)?);
        }
        // get the column names
        let start_s = self.get_column_name("start", true)?;
        let start = start_s.as_str();
        let end_s = self.get_column_name("end", true)?;
        let end = end_s.as_str();
        let seqname_s = self.get_column_name("seqname", true)?;
        let seqname = seqname_s.as_str();
        let strand_s = self.get_column_name("strand", true)?;
        let strand = strand_s.as_str();

        let selected = [start, end, seqname, strand];
        let df = self.df();

        // we build a vector indicating if the start and end of features are valid  
        let valid_rows_df = df.select([start, end, strand])?.lazy()
            .select([
                col(start).gt(lit(0)),
                col(end).gt(lit(0)),
                col(strand).eq(lit("+")).or(col(strand).eq(lit("-"))),
            ])
            .select([col(start).and(col(end)).alias("pos_valid"), col(start).and(col(end)).alias("pos_strand_valid")])
            .collect()?;
        let valid_pos = valid_rows_df.column("pos_valid")?;
        let valid_pos_strand = valid_rows_df.column("pos_strand_valid")?;

        // Then we bail if ignore_invalid is false but we found invalid features
        if !ignore_invalid {
            if ignore_strand {
                // we need to make sure start and end are valid
                if valid_pos.iter().any(|v| v != AnyValue::Boolean(true)) {
                    bail!("Found features with non-positive start/end position. Please remove them first or set ignore_invalid to true.")
                }
            } else {
                // we need to make sure start, end and strand are all valid
                if valid_pos_strand.iter().any(|v| v != AnyValue::Boolean(true)) {
                    bail!("Found features with non-positive start/end/strand position. Please remove them first or set ignore_invalid to true.")
                }
            }
        }

        // [start, stop) Inclusive start, exclusive of stop
        // we define start and end as u64, and we use Vec<String> to store group_by column values
        type Iv = Interval<u64, (usize,Vec<String>)>;

        let mut by_iters = df
            .columns(by)?
            .iter()
            .map(|s| s.iter())
            .collect::<Vec<_>>();

        let mut ess_iters = self.df()
            .columns(&selected)?
            .iter()
            .map(|s| s.iter())
            .collect::<Vec<_>>();

        let valid_rows = if ignore_strand {
            valid_pos
        } else {
            valid_pos_strand
        };

        // we first build the vectors, and then build the lappers using that
        let mut lapper_tree_vec_hm = HashMap::new();

        for (rid, is_valid) in valid_rows.iter().enumerate() {
            // we skip invalid rows as we have bailed already if needed
            if is_valid != AnyValue::Boolean(true) {
                continue;
            }

            // then, we take the start and end
            let s = ess_iters[0]
                .next()
                .expect("should have as many iterations as rows")
                .cast(&DataType::Int64)?
                .try_extract::<i64>()? as u64;

            // we add 1 to the end because rust-lappers uses right-exclusive intervals
            let e = ess_iters[1]
                .next()
                .expect("should have as many iterations as rows")
                .cast(&DataType::Int64)?
                .try_extract::<i64>()? as u64 + 1 ;

            // we take the seqname and strand
            let seqn = ess_iters[2]
                .next()
                .expect("should have as many iterations as rows")
                .to_string();

            let strd = if ignore_strand {
                String::from(".")
            } else {
                ess_iters[3]
                    .next()
                    .expect("should have as many iterations as rows")
                    .to_string()
            };
            
            // we take the by columns
            let mut by_vec = Vec::new();
            for it in by_iters.iter_mut() {
                by_vec.push(it.next().expect("should have as many iterations as rows").to_string());
            }
            let lapper_tree_vec = lapper_tree_vec_hm.entry([seqn.clone(), strd.clone()]).or_insert(Vec::new());
            lapper_tree_vec.push(Iv {
                start: s,
                stop: e,
                val: (rid, by_vec),
            });
        }

        // we build the lappers
        let mut lappers = HashMap::new();

        for (key, lapper_tree_vec) in lapper_tree_vec_hm.into_iter() {
            let lapper = Lapper::new(lapper_tree_vec);
            lappers.insert(key, lapper);
        }


        self.lappers = Some(lappers);
        self.lappers_ignore_strand = Some(ignore_strand);
        let duration: std::time::Duration = start_time.elapsed();
        debug!("build rust-lappers in {:?}", duration);
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

    pub fn write_transcript_sequences<T: AsRef<Path>>(
        &mut self,
        ref_path: T,
        out_path: T,
        exon_name: Option<&str>,
        multithreaded: bool,
        append: bool,
    ) -> anyhow::Result<()> {
        let out_path = out_path.as_ref();

        // if the file exists and append is true, we append to the file
        let out_file = if out_path.try_exists()? && append {
            std::fs::OpenOptions::new()
                .append(true)
                .open(out_path)
                .with_context(|| {
                    format!("Could not open the output file {:?}", out_path.as_os_str())
                })?
        } else {
            // create the folder if it doesn't exist
            fs::create_dir_all(out_path.parent().with_context(|| {
                format!(
                    "Could not get the parent directory of the given output file path {:?}",
                    out_path.as_os_str()
                )
            })?)?;

            // we prepare a fasta writer
            std::fs::File::create(out_path).with_context(|| {
                format!(
                    "Could not create the output file {:?}",
                    out_path.as_os_str()
                )
            })?
        };

        self.validate(false, true)?;
        // get exon_gr
        // exons() ensures that all exon records are valid,
        // and they have a valid exon number
        let mut exon_gr = self.exons(exon_name, multithreaded)?;
        // we build an essential gr for avoiding copying unused columns
        exon_gr.df = exon_gr.df().select([
            exon_gr.get_column_name_str("seqname", true)?,
            exon_gr.get_column_name_str("start", true)?,
            exon_gr.get_column_name_str("end", true)?,
            exon_gr.get_column_name_str("strand", true)?,
            exon_gr.get_column_name_str("transcript_id", true)?,
            exon_gr.get_column_name_str("exon_number", true)?,
        ])?;

        let mut fc = exon_gr.field_columns().clone();
        fc.fix(exon_gr.df(), false)?;
        exon_gr.field_columns = fc;

        let fc = exon_gr.field_columns();
        let seqname = fc.seqname();
        let end = fc.end();
        let transcript_id = fc.transcript_id().unwrap();

        // Now, we read the fasta file and process each reference sequence at a time
        let mut reader = std::fs::File::open(ref_path)
            .map(BufReader::new)
            .map(noodles::fasta::Reader::new)?;
        // we also create a fasta writer
        let mut writer = noodles::fasta::Writer::new(out_file);

        // we iterate the fasta reader. For each fasta reacord (usually chromosome), we do
        // 1. subset the dataframe by the seqname (chromosome name)
        // 2. for each gene, we get the sequence of all its exons
        // 3. for each transcript, we join the transcripts' exon sequences to get the sequence of the transcript
        for result in reader.records() {
            let record = result?;

            let chr_name = record.name().strip_suffix(' ').unwrap_or(record.name());

            let chr_gr = exon_gr.filter(seqname, &[chr_name])?;

            if chr_gr.df().height() == 0 {
                continue;
            }

            // check if exons are in the range of the reference sequence
            if let Some(end_max) = chr_gr.df().column(end)?.i64()?.max() {
                if end_max > record.sequence().len() as i64 {
                    bail!("Found exons that exceed the length of the reference sequence. Cannot proceed")
                }
            } else {
                bail!("Could not get the maximum end value of the exons. Cannot proceed")
            }

            // we get the sequence of a chromosome at a time
            let chr_seq_vec = chr_gr.get_sequences_fasta_record(&record, &OOBOption::Skip)?;

            // we make sure that there is no invalid exon sequences
            if chr_seq_vec
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
            let mut tx_id_iter = chr_gr
                .df()
                .column(transcript_id)?
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

            for (tx_id, seq) in tx_id_iter.zip(chr_seq_vec.into_iter()) {
                if let (Some(tx_id), Some(seq)) = (tx_id, seq) {
                    // first we want to check if the transcript id is the same as the previous one
                    if tx_id == curr_tx {
                        // if it is the same, we extend the exon_vec with the current sequence
                        exon_u8_vec.extend(seq.as_ref().iter());
                    } else {
                        // // if it is not the same, we create a Sequence and push it to seq_vec
                        let definition = Definition::new(curr_tx.clone(), None);
                        let sequence = Sequence::from_iter(exon_u8_vec.clone());

                        writer
                            .write_record(&Record::new(definition, sequence))
                            .with_context(|| {
                                format!(
                                "Could not write the sequence of transcript {} to the output file",
                                curr_tx
                            )
                            })?;
                        exon_u8_vec.clear();
                        exon_u8_vec.extend(seq.as_ref().iter());
                        // update the current transcript id

                        curr_tx = tx_id.to_string();
                    }
                } else {
                    bail!("Found null transcript id or empty exon sequence. This should not happen, please report this bug.")
                }
            }

            // Don't forget our remaining transcript
            // // if it is not the same, we create a Sequence and push it to seq_vec
            let definition = Definition::new(curr_tx.clone(), None);
            let sequence = Sequence::from_iter(exon_u8_vec.clone());

            writer
                .write_record(&Record::new(definition, sequence))
                .with_context(|| {
                    format!(
                        "Could not write the sequence of transcript {} to the output file",
                        curr_tx
                    )
                })?;
            exon_u8_vec.clear();
        }

        Ok(())
    }

    /// Extract the sequence of the features in the Grangers object from the provided reference file.
    /// Currently only fasta file is supported. This function four field columns: seqname, start, end, and strand.
    /// Arguments:
    /// - `genome_path`: the path to the reference genome file.
    /// - `file_format`: the format of the reference genome file. Currently only fasta is supported.
    /// - `oob_option`: the option for out-of-boundary positions. It can be either `Truncate` or `Skip`. If `Truncate`, the out-of-boundary positions will be truncated to the start or end of the sequence. If `Skip`, a None will be returned for features with OOB positions
    /// The function outputs the extracted sequence as a vector of `Option<Sequence>`. If the feature has OOB positions and the oob_option is set as `Skip`, the corresponding element in the vector will be None. The order of the vector follows the row order of the dataframe of the Grangers object.
    pub fn _write_sequences<T: AsRef<Path>>(
        &mut self,
        ref_path: T,
        out_path: T,
        ignore_strand: bool,
        name_column: Option<&str>,
        oob_option: &OOBOption,
    ) -> anyhow::Result<()> {
        let out_path = out_path.as_ref();

        // create the folder if it doesn't exist
        fs::create_dir_all(out_path.parent().with_context(|| {
            format!(
                "Could not get the parent directory of the given output file path {:?}",
                out_path.as_os_str()
            )
        })?)?;

        // we prepare a fasta writer
        let out_file = std::fs::File::create(out_path).with_context(|| {
            format!(
                "Could not create the output file {:?}",
                out_path.as_os_str()
            )
        })?;

        self.validate(false, true)?;

        // if name is invalid, ignore
        let name_column = if let Some(name_column) = name_column {
            if self.get_column_name_str(name_column, true).is_ok() {
                self.get_column_name(name_column, false)?
            } else {
                warn!("The provided name column {:?} for naming the extracted sequences is not in the dataframe. Row order will will be used instead.", name_column);
                "row_order".to_owned()
            }
        } else {
            info!("No name column is provided. The extracted sequences will be named by the row order.");
            "row_order".to_owned()
        };

        let mut fc = self.field_columns().clone();
        // we need to map the sequence back to the original row order of the dataframe
        // So, we need to have a minimum copy of the dataset, which contains only the essential fields,
        // and add one more column representing the row order of the original dataframe
        let selection = [fc.seqname(), fc.start(), fc.end(), fc.strand()];

        let mut df = self.df.select(selection)?;

        df.with_column(Series::new(
            "row_order",
            (0..df.height() as u32).collect::<Vec<u32>>(),
        ))?;

        // if ignore strand, set the strand to +
        if ignore_strand {
            df.with_column(Series::new(fc.strand(), vec!["+"; df.height()]))?;
        }

        fc.fix(&df, false)?;

        let essential_gr = Grangers::new(df, None, None, IntervalType::default(), fc, false)?;

        let seqname = essential_gr.get_column_name_str("seqname", true)?;

        let reader = std::fs::File::open(ref_path).map(BufReader::new)?;
        let mut reader = noodles::fasta::Reader::new(reader);

        let mut writer = noodles::fasta::Writer::new(out_file);
        let mut empty_counter = 0;

        // we iterate the fasta reader. For each fasta reacord (usually chromosome), we do
        // 1. subset the dataframe by the chromosome name
        // 2. get the sequence of the features in the dataframe on that fasta record
        // 3. insert the sequence into the sequence vector according to the row order
        for result in reader.records() {
            let record = result?;

            let chr_name = record.name().strip_suffix(' ').unwrap_or(record.name());
            let chr_gr = essential_gr.filter(seqname, &[chr_name])?;

            if chr_gr.df().height() == 0 {
                continue;
            }
            let name_vec = chr_gr
                .df()
                .column(name_column.as_str())?
                .utf8()?
                .into_iter()
                .map(|s| s.unwrap())
                .collect::<Vec<_>>();
            // we get the sequence of a chromosome at a time
            let chr_seq_vec = chr_gr.get_sequences_fasta_record(&record, oob_option)?;

            // we push seuqence to the correct position
            for (name, sequence) in name_vec.into_iter().zip(chr_seq_vec.into_iter()) {
                let definition = Definition::new(name, None);
                if let Some(sequence) = sequence {
                    writer
                        .write_record(&Record::new(definition, sequence))
                        .with_context(|| {
                            format!(
                                "Could not write sequence {} to the output file; Cannot proceed.",
                                name
                            )
                        })?;
                } else {
                    empty_counter += 1;
                }
            }
        }
        if empty_counter > 0 {
            warn!("Extracted empty sequence from {} records. They are usually caused by out of boundary features.", empty_counter)
        }

        Ok(())
    }

    pub fn write_sequences<T: AsRef<Path>>(
        &mut self,
        ref_path: T,
        out_path: T,
        ignore_strand: bool,
        name_column: Option<&str>,
        oob_option: OOBOption,
        append: bool,
    ) -> anyhow::Result<()> {
        let out_path = out_path.as_ref();

        // if the file exists and append is true, we append to the file
        let out_file = if out_path.try_exists()? && append {
            std::fs::OpenOptions::new()
                .append(true)
                .open(out_path)
                .with_context(|| {
                    format!("Could not open the output file {:?}", out_path.as_os_str())
                })?
        } else {
            // create the folder if it doesn't exist
            fs::create_dir_all(out_path.parent().with_context(|| {
                format!(
                    "Could not get the parent directory of the given output file path {:?}",
                    out_path.as_os_str()
                )
            })?)?;

            // we prepare a fasta writer
            std::fs::File::create(out_path).with_context(|| {
                format!(
                    "Could not create the output file {:?}",
                    out_path.as_os_str()
                )
            })?
        };

        self.validate(false, true)?;

        // if name is invalid, ignore
        let name_column = if let Some(name_column) = name_column {
            if self.get_column_name_str(name_column, true).is_ok() {
                self.get_column_name(name_column, false)?
            } else {
                warn!("The provided name column {:?} for naming the extracted sequences is not in the dataframe. Row order will will be used instead.", name_column);
                "row_order".to_owned()
            }
        } else {
            info!("No name column is provided. The extracted sequences will be named by the row order.");
            "row_order".to_owned()
        };

        let mut fc = self.field_columns().clone();
        // we need to map the sequence back to the original row order of the dataframe
        // So, we need to have a minimum copy of the dataset, which contains only the essential fields,
        // and add one more column representing the row order of the original dataframe
        let selection = [
            fc.seqname(),
            fc.start(),
            fc.end(),
            fc.strand(),
            name_column.as_str(),
        ];

        let mut df = if name_column.as_str() == "row_order" {
            self.df
                .with_row_count("row_order", None)?
                .select(selection)?
        } else {
            self.df.select(selection)?
        };

        // if ignore strand, set the strand to +
        if ignore_strand {
            df.with_column(Series::new(fc.strand(), vec!["+"; df.height()]))?;
        }

        fc.fix(&df, false)?;

        let essential_gr = Grangers::new(df, None, None, IntervalType::default(), fc, false)?;

        let seqname_s = essential_gr.get_column_name("seqname", true)?;
        let seqname = seqname_s.as_str();

        let reader = std::fs::File::open(ref_path).map(BufReader::new)?;
        let mut reader = noodles::fasta::Reader::new(reader);

        let mut writer = noodles::fasta::Writer::new(out_file);
        let mut empty_counter = 0;

        // we iterate the fasta reader. For each fasta reacord (usually chromosome), we do
        // 1. subset the dataframe by the chromosome name
        // 2. get the sequence of the features in the dataframe on that fasta record
        // 3. insert the sequence into the sequence vector according to the row order
        for result in reader.records() {
            let record = result?;

            let chr_name = record.name().strip_suffix(' ').unwrap_or(record.name());
            let chr_gr = essential_gr.filter(seqname, &[chr_name])?;

            if chr_gr.df().height() == 0 {
                continue;
            }

            let name_vec = chr_gr
                .df()
                .column(name_column.as_str())?
                .utf8()?
                .into_iter()
                .map(|s| s.unwrap())
                .collect::<Vec<_>>();

            let chrsi = ChrRowSeqIter::new(&chr_gr, &record, oob_option)?;

            for (feat_name, chrsi_rec) in name_vec.into_iter().zip(chrsi) {
                if let Ok(sequence) = chrsi_rec {
                    let definition = Definition::new(feat_name, None);

                    // we write if the sequence is not empty
                    writer
                        .write_record(&Record::new(definition, sequence))
                        .with_context(|| {
                            format!(
                                "Could not write sequence {} to the output file; Cannot proceed.",
                                feat_name
                            )
                        })?;
                } else {
                    empty_counter += 1;
                }
            }
        }
        if empty_counter > 0 {
            warn!("Unable to extract sequence for {} records. They are usually caused by out of boundary features or an invalid alphabet.", empty_counter)
        }
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
        &mut self,
        fasta_path: T,
        exon_name: Option<&str>,
        multithreaded: bool,
    ) -> anyhow::Result<Vec<noodles::fasta::Record>> {
        self.validate(false, true)?;
        // get exon_gr
        // exons() ensures that all exon records are valid,
        // and they have a valid exon number
        let mut exon_gr = self.exons(exon_name, multithreaded)?;
        // we build an essential gr for avoiding copying unused columns
        exon_gr.df = exon_gr.df().select([
            exon_gr.get_column_name_str("seqname", true)?,
            exon_gr.get_column_name_str("start", true)?,
            exon_gr.get_column_name_str("end", true)?,
            exon_gr.get_column_name_str("strand", true)?,
            exon_gr.get_column_name_str("transcript_id", true)?,
            exon_gr.get_column_name_str("exon_number", true)?,
        ])?;

        let mut fc = exon_gr.field_columns().clone();
        fc.fix(exon_gr.df(), false)?;
        exon_gr.field_columns = fc;

        let fc = exon_gr.field_columns();
        let seqname = fc.seqname();
        let end = fc.end();
        let transcript_id = fc.transcript_id().unwrap();

        // Now, we read the fasta file and process each reference sequence at a time
        let reader = std::fs::File::open(fasta_path).map(BufReader::new)?;
        let mut reader = noodles::fasta::Reader::new(reader);
        // let mut seq_vec: Vec<Option<Sequence>> = vec![None; exon_gr.df().height()];
        let mut transcript_seq_vec: Vec<Record> =
            Vec::with_capacity(self.df().column(transcript_id)?.unique()?.len());

        // we iterate the fasta reader. For each fasta reacord (usually chromosome), we do
        // 1. subset the dataframe by the seqname (chromosome name)
        // 2. for each gene, we get the sequence of all its exons
        // 3. for each transcript, we join the transcripts' exon sequences to get the sequence of the transcript
        for result in reader.records() {
            let record = result?;

            let chr_name = record.name().strip_suffix(' ').unwrap_or(record.name());
            let chr_gr = exon_gr.filter(seqname, &[chr_name])?;

            if chr_gr.df().height() == 0 {
                continue;
            }

            // check if exons are in the range of the reference sequence
            if let Some(end_max) = chr_gr.df().column(end)?.i64()?.max() {
                if end_max > record.sequence().len() as i64 {
                    bail!("Found exons that exceed the length of the reference sequence. Cannot proceed")
                }
            } else {
                bail!("Could not get the maximum end value of the exons. Cannot proceed")
            }

            // we get the sequence of a chromosome at a time
            let chr_seq_vec = chr_gr.get_sequences_fasta_record(&record, &OOBOption::Skip)?;

            // we make sure that there is no invalid exon sequences
            if chr_seq_vec
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
            let mut tx_id_iter = chr_gr
                .df()
                .column(transcript_id)?
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

            for (tx_id, seq) in tx_id_iter.zip(chr_seq_vec.into_iter()) {
                if let (Some(tx_id), Some(seq)) = (tx_id, seq) {
                    // first we want to check if the transcript id is the same as the previous one
                    if tx_id == curr_tx {
                        // if it is the same, we extend the exon_vec with the current sequence
                        exon_u8_vec.extend(seq.as_ref().iter());
                    } else {
                        // if it is not the same, we create a Sequence and push it to seq_vec
                        let definition = Definition::new(curr_tx.clone(), None);
                        let sequence = Sequence::from_iter(exon_u8_vec.clone());
                        transcript_seq_vec.push(Record::new(definition, sequence));
                        exon_u8_vec.clear();
                        exon_u8_vec.extend(seq.as_ref().iter());
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
    pub fn _get_sequences<T: AsRef<Path>>(
        &mut self,
        fasta_path: T,
        ignore_strand: bool,
        name: Option<&str>,
        oob_option: &OOBOption,
    ) -> anyhow::Result<Vec<Option<Record>>>
// anyhow::Result<Vec<fasta::record::Sequence>>
    {
        self.validate(false, true)?;

        // if name is invalid, ignore
        let name = if name.is_some() && self.get_column_name(name.unwrap(), true).is_ok() {
            warn!(
                "The provided name column {:?} is not in the dataframe. Ignored.",
                name
            );
            Some(self.get_column_name(name.unwrap(), false)?)
        } else {
            None
        };

        let mut fc = self.field_columns().clone();
        // we need to map the sequence back to the original row order of the dataframe
        // So, we need to have a minimum copy of the dataset, which contains only the essential fields,
        // and add one more column representing the row order of the original dataframe
        let selection = [fc.seqname(), fc.start(), fc.end(), fc.strand()];

        let df = self
            .df
            .select(selection)?
            .lazy()
            .with_row_count("row_order", None)
            .with_column(if ignore_strand {
                lit("+").alias("strand")
            } else {
                col("strand")
            })
            .collect()?;

        fc.fix(&df, false)?;

        let essential_gr = Grangers::new(df, None, None, IntervalType::default(), fc, false)?;

        let seqname = essential_gr.get_column_name_str("seqname", true)?;

        let reader = std::fs::File::open(fasta_path).map(BufReader::new)?;
        let mut reader = noodles::fasta::Reader::new(reader);

        let mut seq_vec: Vec<Option<Record>> = vec![None; essential_gr.df().height()];
        // we iterate the fasta reader. For each fasta reacord (usually chromosome), we do
        // 1. subset the dataframe by the chromosome name
        // 2. get the sequence of the features in the dataframe on that fasta record
        // 3. insert the sequence into the sequence vector according to the row order
        for result in reader.records() {
            let record = result?;

            let chr_name = record.name().strip_suffix(' ').unwrap_or(record.name());
            let chr_gr = essential_gr.filter(seqname, &[chr_name])?;

            if chr_gr.df().height() == 0 {
                continue;
            }
            let name_vec = if let Some(name) = &name {
                chr_gr
                    .df()
                    .column(name)?
                    .utf8()?
                    .into_iter()
                    .map(|s| s.unwrap())
                    .collect::<Vec<_>>()
            } else {
                Vec::new()
            };
            // we get the sequence of a chromosome at a time
            let chr_seq_vec = chr_gr.get_sequences_fasta_record(&record, oob_option)?;

            // we push seuqence to the correct position
            for (idx, seq) in chr_gr
                .df()
                .column("row_order")?
                .u32()?
                .into_iter()
                .zip(chr_seq_vec.into_iter())
            {
                let idx: usize = idx.unwrap() as usize;
                let seq_name = if name.is_some() {
                    name_vec[idx].to_string()
                } else {
                    idx.to_string()
                };

                let definition = Definition::new(seq_name, None);

                seq_vec[idx] = seq.map(|seq| Record::new(definition, seq));
            }
        }

        Ok(seq_vec)
    }

    pub fn get_sequences<T: AsRef<Path>>(
        &mut self,
        ref_path: T,
        ignore_strand: bool,
        name_column: Option<&str>,
        oob_option: OOBOption,
    ) -> anyhow::Result<Vec<Option<Record>>> {
        self.validate(false, true)?;

        // if name is invalid, ignore
        let name_column = if let Some(name_column) = name_column {
            if self.get_column_name_str(name_column, true).is_ok() {
                self.get_column_name(name_column, false)?
            } else {
                warn!("The provided name column {:?} for naming the extracted sequences is not in the dataframe. Row order will will be used instead.", name_column);
                "row_order".to_owned()
            }
        } else {
            info!("No name column is provided. The extracted sequences will be named by the row order.");
            "row_order".to_owned()
        };

        let mut fc = self.field_columns().clone();
        // we need to map the sequence back to the original row order of the dataframe
        // So, we need to have a minimum copy of the dataset, which contains only the essential fields,
        // and add one more column representing the row order of the original dataframe
        let selection = [
            fc.seqname(),
            fc.start(),
            fc.end(),
            fc.strand(),
            name_column.as_str(),
        ];

        let df = self
            .df
            .select(selection)?
            .lazy()
            .with_row_count("row_order", None)
            .with_column(if ignore_strand {
                lit("+").alias("strand")
            } else {
                col("strand")
            })
            .collect()?;

        // we only use a subset of the columns, so fix fc
        fc.fix(&df, false)?;

        let essential_gr = Grangers::new(df, None, None, IntervalType::default(), fc, false)?;
        let seqname = essential_gr.get_column_name_str("seqname", true)?;
        let mut reader = std::fs::File::open(ref_path)
            .map(BufReader::new)
            .map(noodles::fasta::Reader::new)?;
        // let mut reader = noodles::fasta::Reader::new(reader);

        let mut seq_vec: Vec<Option<Record>> = vec![None; essential_gr.df().height()];

        // we iterate the fasta reader. For each fasta reacord (usually chromosome), we do
        // 1. subset the dataframe by the chromosome name
        // 2. get the sequence of the features in the dataframe on that fasta record
        // 3. insert the sequence into the sequence vector according to the row order
        for result in reader.records() {
            let record = result?;

            let chr_name = record.name().strip_suffix(' ').unwrap_or(record.name());
            let chr_gr = essential_gr.filter(seqname, &[chr_name])?;

            if chr_gr.df().height() == 0 {
                continue;
            }

            let name_vec_iter = chr_gr
                .df()
                .column(name_column.as_str())?
                .cast(&DataType::Utf8)?
                .utf8()?
                .into_iter()
                .map(|s| {
                    s.expect(
                        "The name column contains null values. Please report this bug on GitHub.",
                    )
                    .to_string()
                })
                .collect::<Vec<String>>();

            let row_order_iter = chr_gr
                .df()
                .column("row_order")?
                .u32()?
                .into_iter()
                .map(|s| s.expect("Could not get row order. Please report this bug on GitHub."));

            let chrsi = ChrRowSeqIter::new(&chr_gr, &record, oob_option)?;

            for ((idx, feat_name), sequence) in row_order_iter.zip(name_vec_iter).zip(chrsi) {
                let sequence = match sequence {
                    Ok(sequence) => sequence,
                    Err(e) => {
                        warn!("Failed to get sequence for feature {} at row {}. The error message was {:?}", feat_name, idx, e);

                        seq_vec[idx as usize] = None;
                        continue;
                    }
                };

                let definition = Definition::new(feat_name, None);
                let record = Record::new(definition, sequence);

                seq_vec[idx as usize] = Some(record);
            }
        }

        let empty_counter: usize = seq_vec
            .iter()
            .fold(0usize, |acc, s| acc + s.is_none() as usize);
        if empty_counter > 0 {
            warn!("Unable to extract sequence for {} records. They are usually caused by out of boundary features or an invalid alphabet.", empty_counter)
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
    pub(crate) fn get_sequences_fasta_record(
        &self,
        record: &fasta::record::Record,
        oob_option: &OOBOption,
    ) -> anyhow::Result<Vec<Option<Sequence>>> {
        self.validate(true, true)?;
        let df = self.df();
        let seqname = self.get_column_name("seqname", true)?;
        let start = self.get_column_name("start", true)?;
        let end = self.get_column_name("end", true)?;
        let strand = self.get_column_name("strand", true)?;

        // initialize seq vector
        if df.column(&seqname)?.unique()?.len() > 1 {
            bail!("The dataframe contains more than one reference name. Please filter the dataframe by the reference name first.")
        }

        let mut seq_vec = Vec::with_capacity(df.height());
        let ses = df.columns([start, end, strand])?;
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

struct ChrRowSeqIter<'a> {
    iters: Vec<polars::series::SeriesIter<'a>>,
    record: &'a Record,
    oob_option: OOBOption,
    seqlen: usize,
}

impl<'a> ChrRowSeqIter<'a> {
    pub fn new(
        grangers: &'a Grangers,
        record: &'a Record,
        oob_option: OOBOption,
    ) -> anyhow::Result<Self> {
        let fc = grangers.field_columns();
        let iters: Vec<polars::series::SeriesIter> = vec![
            grangers.df.column(fc.start())?.iter(),
            grangers.df.column(fc.end())?.iter(),
            grangers.df.column(fc.strand())?.iter(),
        ];
        let seqlen = record.sequence().len();
        Ok(Self {
            iters,
            record,
            oob_option,
            seqlen,
        })
    }
}

impl<'a> Iterator for ChrRowSeqIter<'a> {
    type Item = anyhow::Result<Sequence>;
    fn next(&mut self) -> Option<Self::Item> {
        // first we check if we can extract value or not
        if let (Some(start), Some(end), Some(strand)) = (
            self.iters[0].next(),
            self.iters[1].next(),
            self.iters[2].next(),
        ) {
            // the second if check if the start, end and strand are non-null
            let sequence = if let (
                AnyValue::Int64(start),
                AnyValue::Int64(end),
                AnyValue::Utf8(strand),
            ) = (start, end, strand)
            {
                // we need to convert the start and end to trunacated one if oob_option is Truncate
                let (start, end) = if self.oob_option == OOBOption::Truncate {
                    (
                        noodles::core::Position::new(std::cmp::max(1, start as usize)),
                        noodles::core::Position::new(std::cmp::min(self.seqlen, end as usize)),
                    )
                } else {
                    (
                        noodles::core::Position::new(start as usize),
                        noodles::core::Position::new(end as usize),
                    )
                };
                // the third if check if the start and end are non-negative
                if let (Some(start), Some(end)) = (start, end) {
                    let seq = self.record.sequence().get(start..=end);
                    // the fourth if check if the sequence is valid
                    if let Some(seq) = seq {
                        let mut sequence = Ok(Sequence::from_iter(seq.iter().copied()));
                        if strand == "-" {
                            sequence = sequence.unwrap().complement().rev().collect::<Result<_, _>>().with_context(||"Could not get the reverse complement of a sequence. Please check if the alphabet is valid.");
                        };
                        sequence
                    } else {
                        // we can't get a slice from the start to the end
                        Err(anyhow::anyhow!("Could not get the sequence from the start to the end. Please check if the start and end are valid."))
                    }
                } else {
                    // if start or end is negative, we return None
                    Err(anyhow::anyhow!("Found invalid start or end. Please check if the start and end are within boundary."))
                }
            } else {
                // if we can't get valid start, end or strand, then we return None
                Err(anyhow::anyhow!("Found null start, end or strand. Please check if the start, end and strand are valid."))
            };
            Some(sequence)
        } else {
            // if they are none, then we reach the end of the iterator
            None
        }
    }
}

#[allow(dead_code)]
pub fn argsort1based<T: Ord>(data: &[T], descending: bool) -> Vec<usize> {
    let mut indices = (1..=data.len()).collect::<Vec<_>>();
    indices.sort_by_key(|&i| &data[i - 1]);
    if descending {
        indices.reverse();
    }
    indices
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
            return Result::<Option<polars::prelude::Series>, PolarsError>::Err(
                PolarsError::ComputeError(
                    "Found missing value in the start or end column. Cannot proceed.".into(),
                ),
            );

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

    use crate::grangers::reader::gtf::{AttributeMode, Attributes, GStruct};
    use noodles::core::Position;

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
            IntervalType::Inclusive(1),
            FieldColumns::default(),
            false,
        )
        .unwrap();

        if SAY {
            println!("gr: {:?}", gr.df());
        }

        // default setting
        let gr1: Grangers = gr
            .merge(&["seqname", "gene_id"], false, None, None)
            .unwrap();

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
        let gr1 = gr.merge(&["seqname"], true, None, None).unwrap();

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

        let gr1: Grangers = gr
            .merge(&["seqname", "gene_id"], false, Some(0), None)
            .unwrap();

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
        let gr1: Grangers = gr
            .merge(&["seqname", "gene_id"], false, Some(2), None)
            .unwrap();

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
            IntervalType::Inclusive(1),
            FieldColumns::default(),
            false,
        )
        .unwrap();

        if SAY {
            println!("gr: {:?}", gr.df());
        }

        // default setting
        let gr1: Grangers = gr.gaps(&["seqname", "gene_id"], false, None, None).unwrap();

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
            IntervalType::Inclusive(1),
            FieldColumns::default(),
            false,
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
            IntervalType::Inclusive(1),
            FieldColumns::default(),
            false,
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
        let exon_number = vec![1i64, 2, 1, 2];
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
            IntervalType::Inclusive(1),
            FieldColumns::default(),
            false,
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
        let exon_number = vec![1i64, 2, 1, 2];
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
    fn test_introns() {
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
            IntervalType::Inclusive(1),
            FieldColumns::default(),
            false,
        )
        .unwrap();
        if SAY {
            println!("gr: {:?}", gr.df());
        }

        // extend from both
        let gr1 = gr.introns(None, None, None, true).unwrap();
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
        let gr1 = gr.introns(None, None, None, true).unwrap();
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

    #[test]
    fn test_get_sequences_fasta_record() {
        let df = df!(
            "seqname" => ["chr1", "chr1"],
            "start" => [10i64, 15],
            "end" => [20i64, 30],
            "strand"=> ["+", "-"],
            "gene_id" => ["g1", "g2"],
        )
        .unwrap();

        let gr = Grangers::new(
            df,
            None,
            None,
            IntervalType::Inclusive(1),
            FieldColumns::default(),
            false,
        )
        .unwrap();

        let definition = Definition::new("chr1", None);
        let sequence = Sequence::from(
            b"GTAGTTCTCTGGGACCTGCAAGATTAGGCAGGGACATGTGAGAGGTGACAGGGACCTGCA".to_vec(),
        );
        let record = Record::new(definition, sequence.clone());

        let seq_vec = gr
            .get_sequences_fasta_record(&record, &OOBOption::Skip)
            .unwrap();
        let start = Position::new(10).unwrap();
        let end = Position::new(20).unwrap();
        let expected_seq1 = sequence.slice(start..=end).unwrap();
        let start = Position::new(15).unwrap();
        let end = Position::new(30).unwrap();
        let expected_seq2 = sequence
            .slice(start..=end)
            .unwrap()
            .complement()
            .rev()
            .collect::<Result<Sequence, _>>()
            .unwrap();
        assert_eq!(seq_vec[0], Some(expected_seq1));
        assert_eq!(seq_vec[1], Some(expected_seq2));
    }

    #[test]
    fn test_chrrowseqiter() {
        let df = df!(
            "seqname" => ["chr1", "chr1"],
            "start" => [10i64, 15],
            "end" => [20i64, 30],
            "strand"=> ["+", "-"],
            "gene_id" => ["g1", "g2"],
        )
        .unwrap();

        let gr = Grangers::new(
            df,
            None,
            None,
            IntervalType::Inclusive(1),
            FieldColumns::default(),
            false,
        )
        .unwrap();

        let definition = Definition::new("chr1", None);
        let sequence = Sequence::from(
            b"GTAGTTCTCTGGGACCTGCAAGATTAGGCAGGGACATGTGAGAGGTGACAGGGACCTGCA".to_vec(),
        );
        let record = Record::new(definition, sequence.clone());

        let mut chrsi = ChrRowSeqIter::new(&gr, &record, OOBOption::Skip).unwrap();

        let chrsi1 = chrsi.next().unwrap().unwrap();
        let chrsi2 = chrsi.next().unwrap().unwrap();
        assert!(chrsi.next().is_none());

        let start = Position::new(10).unwrap();
        let end = Position::new(20).unwrap();
        let expected_seq1 = Sequence::from_iter(sequence.get(start..=end).unwrap().iter().cloned());
        assert_eq!(chrsi1, expected_seq1);

        let start = Position::new(15).unwrap();
        let end = Position::new(30).unwrap();
        let expected_seq2 = Sequence::from_iter(sequence.get(start..=end).unwrap().iter().cloned())
            .complement()
            .rev()
            .collect::<Result<Sequence, _>>()
            .unwrap();
        assert_eq!(chrsi2, expected_seq2);
    }

    #[test]
    fn test_write_gtf() {
        let df = df!(
            "seqname" => ["chr1", "chr1"],
            "start" => [10i64, 15],
            "end" => [20i64, 30],
            "strand"=> ["+", "-"],
            "source" => [Some("HAVANA"), None],
            "gene_id" => ["g1", "g2"],
            "gene_name" => [Some("gene_1"), None],
        )
        .unwrap();

        let gr = Grangers::new(
            df,
            None,
            None,
            IntervalType::Inclusive(1),
            FieldColumns::default(),
            false,
        )
        .unwrap();

        if SAY {
            println!("gr1: {:?}", gr.df());
        }

        let gtf_df = gr.get_gtf_df("a_fake_file").unwrap();

        if SAY {
            println!("gtf_df: {:?}", gtf_df);
        }

        let source = vec![String::from("HAVANA"), String::from(".")];
        let gtf_df_attributes = gtf_df
            .column("attributes")
            .unwrap()
            .utf8()
            .unwrap()
            .into_iter()
            .map(|x| x.unwrap())
            .collect::<Vec<&str>>();

        // we cannot make sure that the order of the attributes are always the same
        assert!(gtf_df_attributes[0].contains("gene_name \"gene_1\";"));
        assert!(gtf_df_attributes[0].contains("gene_id \"g1\";"));
        assert!(gtf_df_attributes[1].contains("gene_id \"g2\";"));

        assert_eq!(
            gtf_df
                .column("source")
                .unwrap()
                .utf8()
                .unwrap()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<&str>>(),
            source
        );
    }


    #[test]
    fn test_build_lappers() {
        let df = df!(
            "seqname" => ["chr1", "chr1", "chr1", "chr2", "chr2", "chr2", "chr2"],
            "feature_type" => ["exon", "exon", "exon", "exon", "exon", "exon", "exon"],
            "start" => [1i64, 21, -5, 1, 51, 1, 51],
            "end" => [10i64, 30, 5, 100, 150, 100, 150],
            "strand"=> ["+", "+", "+", "+", "+", "-", "-"],
            "gene_id" => ["g1", "g1", "g1", "g2", "g2", "g2", "g2"],
        )
        .unwrap();

        let mut gr = Grangers::new(
            df,
            None,
            None,
            IntervalType::Inclusive(1),
            FieldColumns::default(),
            false,
        )
        .unwrap();

        if SAY {
            println!("gr1: {:?}", gr.df());
        }

        gr.build_lapper(true, false, true, &["gene_id"]).unwrap();

        println!("{:?}", gr.lappers);

        // assert!(false);

    }
}
