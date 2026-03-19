use crate::grangers_utils;
// TODO:
// 1. update write and get sequence functions to use the same implementation
use crate::grangers_utils::*;
use crate::options::*;
use crate::reader;
use crate::reader::fasta::SeqInfo;
use anyhow::{bail, Context};
use lazy_static::lazy_static;
pub(crate) use noodles::fasta::record::{Definition, Sequence};
use nutype::nutype;
use polars::{frame::DataFrame, lazy::prelude::*, prelude::*, series::Series};
use rust_lapper::{Interval, Lapper};
// use tracing::field;
use std::collections::{HashMap, HashSet};
use std::convert::AsRef;
use std::fs;
use std::io::{BufWriter, Read, Write};
use std::iter::IntoIterator;
use std::ops::FnMut;
use std::ops::{Add, Mul, Sub};
use std::path::Path;
use std::result::Result::Ok;
use std::sync::atomic::{AtomicU32, Ordering};
use tracing::debug;
use tracing::{info, warn};


// we give each grangers struct a unique
// program identifier which is the order in
// which it was created.
lazy_static! {
    static ref GRANGERS_COUNTER: AtomicU32 = AtomicU32::new(0);
}

type LapperType = Lapper<u64, (usize, Vec<String>)>;

#[nutype(derive(Debug, Clone, AsRef))]
/// The record ID of a Grangers record.
pub struct GrangersRecordID(u32);

/// Represents a collection of genomic sequences, each associated with a unique identifier.
///
/// This structure is used to manage and manipulate collections of genomic sequences, typically
/// derived from FASTA files. Each sequence is stored alongside its unique identifier to facilitate
/// easy retrieval and referencing throughout genomic analysis workflows.
///
/// ### Fields
///
/// * `records` - A vector of tuples, each containing a [`GrangersRecordID`] and a [`noodles::fasta::Record`].
///    The [`GrangersRecordID`] serves as a unique identifier for each genomic sequence, while the
///    [`noodles::fasta::Record`] contains the actual sequence data and related metadata as defined by
///    the `noodles` crate, a Rust library for handling bioinformatics formats.
///
/// * `signature` - A 64-bit unsigned integer used as a unique signature for the entire collection.
///    This can be used to verify the integrity of the data or to quickly compare this collection with others.
///
pub struct GrangersSequenceCollection {
    pub records: Vec<(GrangersRecordID, noodles::fasta::Record)>,
    pub signature: u64,
}

impl GrangersSequenceCollection {
    /// Creates a new [`GrangersSequenceCollection`] with a specified signature and initial capacity.
    ///
    /// This constructor initializes a new instance of [`GrangersSequenceCollection`] with an empty
    /// vector of sequence records. The vector is preallocated with the specified capacity to optimize
    /// memory allocations when the expected number of sequences is known in advance. The collection
    /// is also assigned a unique signature for identification.
    ///
    /// ### Arguments
    ///
    /// * `signature` - A 64-bit unsigned integer used as the unique signature or identifier for the collection.
    /// * `capacity` - The initial capacity of the vector holding the sequence records. This value determines
    ///   how many sequence records the vector can hold before needing to reallocate memory.
    ///
    /// ### Returns
    ///
    /// A new instance of [GrangersSequenceCollection] with the specified signature and initial capacity.
    ///
    /// ### Example
    ///
    /// ```
    /// let signature = 123456789;
    /// let capacity = 100;
    /// let sequence_collection = GrangersSequenceCollection::new_with_signature_and_capacity(signature, capacity);
    /// ```
    pub fn new_with_signature_and_capacity(signature: u64, capacity: usize) -> Self {
        GrangersSequenceCollection {
            records: Vec::<(GrangersRecordID, noodles::fasta::Record)>::with_capacity(capacity),
            signature,
        }
    }

    /// Adds a new genomic sequence record to the collection.
    ///
    /// This method appends a new sequence record, consisting of a unique identifier ([GrangersRecordID])
    /// and a FASTA sequence record ([noodles::fasta::Record]), to the end of the collection. It allows
    /// populating the collection with genomic sequence data for analysis or processing.
    ///
    /// ### Arguments
    ///
    /// * `rec_id` - The unique identifier for the new sequence record as a [`GrangersRecordID`].
    /// * `rec` - The genomic sequence information contained in a [`noodles::fasta::Record`].
    ///
    /// ### Returns
    ///
    /// This method does not return any value. It modifies the collection in place by adding the new record.
    ///
    /// ### Example
    ///
    /// ```
    /// let mut collection = GrangersSequenceCollection::new_with_signature_and_capacity(123456789, 10);
    /// let rec_id = GrangersRecordID::new(1);
    /// let fasta_record = noodles::fasta::Record::new("seq1", "ACTG".as_bytes());
    /// collection.add_record(rec_id, fasta_record);
    /// ```
    pub fn add_record(&mut self, rec_id: GrangersRecordID, rec: noodles::fasta::Record) {
        self.records.push((rec_id, rec));
    }

    /// Returns an iterator over the genomic sequence records in the collection.
    ///
    /// This method provides access to all the sequence records contained within the collection without
    /// modifying them. It can be used to iterate over the sequence data for read-only operations such
    /// as analysis or reporting.
    ///
    /// ### Returns
    ///
    /// A non-mutable iterator ([std::slice::Iter]) over the sequence records in the collection. Each
    /// item in the iterator is a reference to a tuple containing a [`GrangersRecordID`] and a
    /// [`noodles::fasta::Record`].
    ///
    /// ### Example
    ///
    /// ```
    /// let collection = GrangersSequenceCollection::new_with_signature_and_capacity(123456789, 10);
    /// for (rec_id, rec) in collection.records_iter() {
    ///     println!("Record ID: {}, Sequence: {}", rec_id, std::str::from_utf8(rec.sequence()).unwrap());
    /// }
    /// ```
    pub fn records_iter(&self) -> std::slice::Iter<'_, (GrangersRecordID, noodles::fasta::Record)> {
        self.records.iter()
    }

    /// Returns a mutable iterator over the genomic sequence records in the collection.
    ///
    /// This method allows iterating over all the sequence records in the collection with the ability
    /// to modify them. It can be used for operations that need to alter the sequence data, such as
    /// updating metadata or correcting sequences.
    ///
    /// ### Returns
    ///
    /// A mutable iterator ([std::slice::IterMut]) over the sequence records in the collection. Each
    /// item in the iterator is a mutable reference to a tuple containing a [`GrangersRecordID`] and a
    /// [`noodles::fasta::Record`].
    ///
    /// ### Example
    ///
    /// ```
    /// let mut collection = GrangersSequenceCollection::new_with_signature_and_capacity(123456789, 10);
    /// for (rec_id, rec) in collection.records_iter_mut() {
    ///     if rec_id == &GrangersRecordID::new(1) {
    ///         rec.set_sequence("GATC".as_bytes());
    ///     }
    /// }
    /// ```
    pub fn records_iter_mut(
        &mut self,
    ) -> std::slice::IterMut<'_, (GrangersRecordID, noodles::fasta::Record)> {
        self.records.iter_mut()
    }

    /* -- this would induce a move of self, which makes the collection essentially
     * useless afterward. If the user want's to have an IntoIter over the records, they
     * should destructure the seq collection.
     */
    /*
    pub fn records_into_iter(self) -> std::vec::IntoIter<(GrangersRecordID, noodles::fasta::Record)> {
        self.records.into_iter()
    }
    */

    /// Retrieves the unique signature of the sequence collection.
    ///
    /// This method returns the signature of the collection, which is a unique identifier assigned during
    /// the creation of the collection. The signature can be used to differentiate this collection from others.
    ///
    /// ### Returns
    ///
    /// A 64-bit unsigned integer (`u64') representing the signature of this
    /// [GrangersSequenceCollection]
    pub fn get_signature(&self) -> u64 {
        self.signature
    }
}

/// Represents a Grangers structure containing genomic annotations and related information.
///
/// This struct is designed to hold and manage genomic data and annotations within a Polars dataframe,
/// along with additional metadata and reference genome information. It includes mechanisms for
/// tracking changes and ensuring data integrity through a unique signature.
///
/// ### Fields
///
/// * `df`: The underlying [DataFrame] from the Polars library, recording all annotations
///   and genomic data. This serves as the primary container for the genomic information.
///
/// * `misc`: An optional [`HashMap<String, Vec<String>>`] for storing additional information
///   or metadata related to the genomic data. Each key represents a metadata category with
///   an associated list of string values.
///
/// * `seqinfo`: An optional [`SeqInfo`] struct containing reference genome information, such
///   as chromosome names and sizes. This information is critical for genomic data analyses
///   and comparisons.
///
/// * `interval_type`: An [`IntervalType`] enum specifying the type of genomic intervals represented
///   in the data frame, such as ranges of genomic coordinates.
///
/// * `field_columns`: A [`FieldColumns`] struct specifying the names of the columns in the data frame
///   that are used to identify genomic features, such as gene names or chromosome positions.
///
/// * `signature`: A [`u64`] serving as the global (process-unique) signature for this [`Grangers`]
///   struct. The upper 32 bits assign each [`Grangers`] struct a distinct (sequential) number at
///   construction, while the lower 32 bits are a version field that is incremented upon each
///   mutating operation of the frame, ensuring traceability and integrity of the data.
///
/// **Notice** that Granges uses 1-based closed intervals for the ranges.
/// If your ranges are not like this, when instantiating new Grangers,
/// you should use the `interval_type` parameter to help the builder
/// to convert the ranges to 1-based closed intervals.
#[derive(Clone)]
pub struct Grangers {
    /// The underlying Polars dataframe recording all annotations
    pub df: DataFrame,
    /// The additional information (metadata)
    pub misc: Option<HashMap<String, Vec<String>>>,
    /// The reference genome information
    pub seqinfo: Option<SeqInfo>,
    /// The interval type
    pub interval_type: IntervalType,
    /// The name of the columns that are used to identify the genomic features
    pub field_columns: FieldColumns,
    /// The global (process-unique) signature assigned to this
    /// grangers struct. It is a u64 where the upper 32-bits
    /// assign each grangers struct a distinct (sequential) number
    /// at construction, and the lower 32-bits are a version field that
    /// is incremented upon each mutating operation of the frame.
    pub signature: u64,
}

impl Grangers {
    #[inline(always)]
    fn inc_signature(&mut self) {
        self.signature += 1;
    }
}

// IO
impl Grangers {
    /// Checks for null values in specified fields of the dataframe.
    ///
    /// This method validates if the specified fields in the dataframe contain any null values.
    /// It uses the `field_columns` struct to ensure the fields are valid before performing the null check.
    /// If null values are found, the function can optionally warn the user and/or halt the operation,
    /// based on the provided boolean flags.
    ///
    /// ### Generics
    ///
    /// * `T`: A type that implements [`AsRef<str>`], allowing for flexible string references as field names.
    ///
    /// ### Arguments
    ///
    /// * `fields`: A slice of items of type `T`, representing the names of the fields to check for null values.
    /// * `is_warn`: A boolean indicating whether to issue warnings if null values are found in the specified fields.
    /// * `is_bail`: A boolean indicating whether to bail out (return an error) if null values are found in the specified fields.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<bool>`]:
    /// * [`Ok`]`(true)` if any null values are found in the specified fields.
    /// * [`Ok`]`(false)` if no null values are found in the specified fields.
    /// * [`Err`] if there is an error during the validation process or if `is_bail` is true and null values are found.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let grangers = Grangers::new(...); // Assume Grangers instance is created
    /// let result = grangers.any_nulls(&["seqname", "start", "end"], true, true)?;
    /// assert!(!result); // No nulls in the specified fields
    /// ```
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

    /// Creates a new instance of [`Grangers`] with the provided configuration.
    ///
    /// This constructor initializes a [`Grangers`] struct, performing several validation and adjustment
    /// steps to ensure the integrity and consistency of the provided data. It validates field columns,
    /// adjusts data frame intervals based on the specified interval type, and ensures that there are no
    /// null values in essential fields.
    ///
    /// ### Arguments
    ///
    /// * `df`: The primary [`DataFrame`] containing the genomic annotations and data.
    /// * `seqinfo`: Optional [`SeqInfo`] providing reference genome information.
    /// * `misc`: Optional [`HashMap<String, Vec<String>>`] for storing additional metadata.
    /// * `interval_type`: The [`IntervalType`] dictating how genomic intervals should be interpreted.
    /// * `field_columns`: [`FieldColumns`] specifying the names of essential columns in the data frame.
    /// * `verbose`: A boolean flag that, if set to true, enables verbose logging.
    ///
    /// ### Returns
    ///
    /// Returns an [`Result<Grangers>`]:
    /// * [Ok]`(Grangers)`: A new [`Grangers`] instance if all validations pass and no critical errors occur.
    /// * [Err]`(...)`: An error encapsulated within an [`anyhow::Error`] if validations fail or if adjustments
    ///    to the data frame encounter issues.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let df = DataFrame::new(vec![...]); // Assume DataFrame is created
    /// let grangers = Grangers::new(df, None, None, IntervalType::Inclusive, FieldColumns::default(), true)?;
    /// ```
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

        let gid = GRANGERS_COUNTER.fetch_add(1, Ordering::SeqCst) as u64;

        // instantiate a new Grangers struct
        let gr = Grangers {
            df,
            misc,
            seqinfo,
            interval_type,
            field_columns,
            signature: (gid << 32),
        };

        // validate
        gr.any_nulls(&gr.field_columns().essential_fields(), verbose, true)?;

        Ok(gr)
    }

    /// Constructs a [`Grangers`] instance from a [`reader::GStruct`] object, typically representing parsed genomic data.
    ///
    /// This function converts genomic data contained within a [`reader::GStruct`] (a common data structure for holding
    /// genomic information) into a [`Grangers`] instance, suitable for further analysis and processing. It involves
    /// transforming genomic attributes into a Polars `DataFrame`, setting appropriate data types, and ensuring the
    /// data aligns with the expected [`Grangers`] structure.
    ///
    /// ### Arguments
    ///
    /// * `gstruct`: A [`reader::GStruct`] containing the raw genomic data to be transformed.
    /// * `interval_type`: The [`IntervalType`] specifying how interval data (start and end positions) should be interpreted.
    ///
    /// ### Returns
    ///
    /// Returns an [`Result<Grangers>`]:
    /// * [Ok]`(Grangers)`: A new [`Grangers`] instance created from the [reader::GStruct] data.
    /// * [Err]`(...)`: An error encapsulated within an `Error` if the data frame creation or [`Grangers`] initialization fails.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let gstruct = reader::GStruct::from_gtf(...); // Assume GStruct is created from a GTF file
    /// let grangers = Grangers::from_gstruct(gstruct, IntervalType::Inclusive)?;
    /// ```
    pub fn from_gstruct(
        gstruct: reader::GStruct,
        interval_type: IntervalType,
    ) -> anyhow::Result<Grangers> {
        // create dataframe!
        // we want to make some columns categorical because of this https://docs.rs/polars/latest/polars/docs/performance/index.html
        // fields
        let mut df_vec = vec![
            Column::new("seqname".into(), gstruct.seqid),
            Column::new("source".into(), gstruct.source),
            Column::new("feature_type".into(), gstruct.feature_type),
            Column::new("start".into(), gstruct.start),
            Column::new("end".into(), gstruct.end),
            Column::new("score".into(), gstruct.score),
            Column::new("strand".into(), gstruct.strand),
            Column::new("phase".into(), gstruct.phase),
        ];

        let el = if let Some(ref extra) = gstruct.attributes.extra {
            extra.len()
        } else {
            0_usize
        };
        df_vec.reserve(df_vec.len() + gstruct.attributes.essential.len() + el);

        //for essential attributes
        for (k, v) in gstruct.attributes.essential {
            if !v.is_empty() {
                df_vec.push(Column::new(k.into(), v));
            };
        }

        // for extra attributes
        if let Some(attributes) = gstruct.attributes.extra {
            for (k, v) in attributes {
                let s = if v.is_empty() {
                    Column::full_null(k.into(), gstruct.attributes.tally, &DataType::String)
                } else {
                    Column::new(k.into(), v)
                };
                df_vec.push(s);
            }
        }

        let df = DataFrame::new(df_vec)?;
        Grangers::new(
            df,
            None,
            gstruct.misc,
            interval_type,
            FieldColumns::default(),
            true,
        )
    }

    /// Constructs a [Grangers] instance from a GTF (GFF2) file.
    ///
    /// This function reads genomic data from a GTF (Gene Transfer Format) file, converts it into a [reader::GStruct] object,
    /// and then transforms this data into a [Grangers] instance. The function allows for the exclusion of non-essential
    /// attributes from the final data structure based on the `only_essential` flag.
    ///
    /// ### Arguments
    ///
    /// * `file_path`: An [`AsRef<std::path::Path>`] specifying the location of the GTF file to be read.
    /// * `only_essential`: A boolean flag indicating whether only essential attributes should be included in the final [Grangers] object.
    ///    If `true`, only essential genomic attributes are included, reducing memory usage and potentially improving performance.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<Grangers>`]:
    /// * [Ok]`(Grangers)`: A [Grangers] instance created from the specified GTF file.
    /// * [Err]`(...)`: An error encapsulated within an [anyhow::Error] if there are issues reading the file or constructing the [Grangers] instance.
    ///
    /// ### Example
    ///
    /// ```rust
    /// // Assuming `path` is a Path to a GTF file.
    /// let grangers = Grangers::from_gtf(&path, true)?;
    /// ```
    ///
    /// ### Errors
    ///
    /// This function may return an error if:
    /// * There is an issue reading or parsing the GTF file.
    /// * The conversion process from [reader::GStruct] to [Grangers] fails, such as due to data inconsistencies or internal validation errors.
    pub fn from_gtf<P: AsRef<std::path::Path>>(
        file_path: P,
        only_essential: bool,
    ) -> anyhow::Result<Grangers> {
        let am = reader::AttributeMode::from(!only_essential);
        let gstruct = reader::GStruct::from_gtf(file_path.as_ref(), am)?;
        Grangers::from_gstruct(gstruct, IntervalType::Inclusive(1))
    }

    /// Constructs a [Grangers] instance from a GFF3 file.
    ///
    /// This function reads genomic data from a GFF3 (General Feature Format) file, converts it into a [reader::GStruct] object,
    /// and then transforms this data into a [Grangers] instance. The function allows for the exclusion of non-essential
    /// attributes from the final data structure based on the `only_essential` flag.
    ///
    /// ### Arguments
    ///
    /// * `file_path`: An [`AsRef<std::path::Path>`] specifying the location of the GTF file to be read.
    /// * `only_essential`: A boolean flag indicating whether only essential attributes should be included in the final [Grangers] object.
    ///    If `true`, only essential genomic attributes are included, reducing memory usage and potentially improving peformance.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<Grangers>`]:
    /// * [Ok]`(Grangers)`: A `Grangers` instance created from the specified GFF file.
    /// * [Err]`(...)`: An error encapsulated within an [`anyhow::Error`] if there are issues reading the file or constructing the [Grangers] instance.
    /// ### Example
    ///
    /// ```rust
    /// // Assuming `path` is a Path to a GFF file.
    /// let grangers = Grangers::from_gff(&path, true)?;
    /// ```
    ///
    /// ### Errors
    ///
    /// This function may return an error if:
    /// * There is an issue reading or parsing the GFF file.
    /// * The conversion process from [reader::GStruct] to [Grangers] fails, such as due to data inconsistencies or internal validation errors.
    pub fn from_gff<P: AsRef<std::path::Path>>(
        file_path: P,
        only_essential: bool,
    ) -> anyhow::Result<Grangers> {
        let am = reader::AttributeMode::from(!only_essential);
        let gstruct = reader::GStruct::from_gff(file_path, am)?;
        Grangers::from_gstruct(gstruct, IntervalType::Inclusive(1))
    }

    // TODO: add the part about making/taking and checking the seqinfo
    pub fn add_seqinfo<T: AsRef<Path>>(&mut self, genome_file: T) -> anyhow::Result<()> {
        self.seqinfo = Some(SeqInfo::from_fasta(genome_file)?);
        Ok(())
    }

    /// Generates a [DataFrame] in GTF format from the [Grangers] instance's current data.
    ///
    /// This method processes the internal [DataFrame], organizing and formatting it to adhere to the GTF (Gene Transfer Format) specification.
    /// It involves categorizing and handling different types of columns: existing fields, missing fields, and attribute columns.
    /// Existing fields are directly transferred, missing fields are filled with default values, and attribute columns are merged into
    /// a single 'attributes' column in GTF format.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<DataFrame>`]:
    /// * [Ok]`(DataFrame)`: A new DataFrame formatted according to the GTF specifications, ready for export or further processing.
    /// * [Err]`(...)`: An error occurred during the DataFrame transformation process.
    ///
    /// ### Example
    ///
    /// ```rust
    /// // Assuming `grangers` is an instance of `Grangers`.
    /// let gtf_df = grangers.get_gtf_df()?;
    /// ```
    ///
    /// ### Errors
    ///
    /// This function may return an error if:
    /// * There is an issue retrieving or updating column names based on the GTF specification.
    /// * There is a failure in transforming the DataFrame, such as during column selection, null filling, or data type casting.
    pub fn get_gtf_df(&self) -> anyhow::Result<DataFrame> {
        // get a copy of the dataframe
        let df = self.df();
        let mut fc = self.field_columns.clone();
        let mut attr_cols: HashSet<&str> = df
            .get_column_names()
            .into_iter()
            .map(|s| s.as_str())
            .collect();
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
                .then(
                    lit(col_name) + lit(" \"") + col(col_name).cast(DataType::String) + lit("\";"),
                )
                .otherwise(lit("")))
            .alias(col_name)
        }));

        // then, we prepare the final dataframe for polar csv writer
        let out_df = self
            .df()
            .clone()
            .lazy()
            .select(expr_vec)
            .select([
                col(fc.field("seqname").expect(
                    "Could not get the seqname field. Please report this issue via GitHub.",
                )),
                col(fc.field("source").expect(
                    "Could not get the source field. Please report this issue via GitHub.",
                )),
                col(fc.field("feature_type").expect(
                    "Could not get the feature_type field. Please report this issue via GitHub.",
                )),
                col(fc
                    .field("start")
                    .expect("Could not get the start field. Please report this issue via GitHub.")),
                col(fc
                    .field("end")
                    .expect("Could not get the end field. Please report this issue via GitHub.")),
                col(fc
                    .field("score")
                    .expect("Could not get the score field. Please report this issue via GitHub.")),
                col(fc.field("strand").expect(
                    "Could not get the strand field. Please report this issue via GitHub.",
                )),
                col(fc
                    .field("phase")
                    .expect("Could not get the phase field. Please report this issue via GitHub.")),
                concat_str(
                    attr_cols.iter().map(|&c| col(c)).collect::<Vec<_>>(),
                    "",
                    false,
                )
                .alias("attributes"),
            ])
            .fill_nan(lit("."))
            .fill_null(lit("."))
            .collect()?;

        Ok(out_df)
    }

    /// Writes the [Grangers] instance's data as a GTF formatted file to the specified path.
    ///
    /// This method exports the genomic data contained within the [Grangers] instance into a GTF (Gene Transfer Format) file.
    /// It first ensures that the output directory exists, then retrieves the internal DataFrame formatted according to GTF standards,
    /// and finally writes this data to the specified file path without including the header and using tab as the separator.
    ///
    /// ### Generics
    ///
    /// * `T`: A type that implements [`AsRef<Path>`], allowing for flexible path references to the output file.
    ///
    /// ### Arguments
    ///
    /// * `file_path`: The path where the GTF file will be written. This can be any type that implements the [`AsRef<Path>`] trait, such as `&str`, `String`, `Path`, or `PathBuf`.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<()>`](anyhow::Result) indicating the outcome of the file writing operation:
    /// * [Ok]`(())`: Successfully wrote the GTF file.
    /// * [Err]`(...)`: An error occurred during the file creation or data writing process.
    ///
    /// ### Example
    ///
    /// ```rust
    /// // Assuming `grangers` is an instance of `Grangers`.
    /// let file_path = PathBuf::from("output.gtf");
    /// grangers.write_gtf(&file_path)?;
    /// ```
    ///
    /// ### Errors
    ///
    /// This function may return an error if:
    /// * The specified output directory cannot be created or accessed.
    /// * There is an issue converting the internal DataFrame to GTF format.
    /// * There is a problem opening or writing to the specified file path.
    pub fn write_gtf<T: AsRef<Path>>(&self, file_path: T) -> anyhow::Result<()> {
        let file_path = file_path.as_ref();

        // create the folder if it doesn't exist
        fs::create_dir_all(file_path.parent().with_context(|| {
            format!(
                "Could not get the parent directory of the given output file path {:?}",
                file_path.as_os_str()
            )
        })?)?;

        let mut out_df = self.get_gtf_df()?;

        let file = std::fs::File::create(file_path)?;
        let mut file = BufWriter::with_capacity(4194304, file);
        CsvWriter::new(&mut file)
            .include_header(false)
            .with_separator(b'\t')
            .with_null_value(".".to_string())
            .finish(&mut out_df)?;

        Ok(())
    }
}

// get struct fields
impl Grangers {
    /// Provides a reference to the [FieldColumns] of the Grangers instance.
    ///
    /// This method allows access to the [FieldColumns] struct, which specifies the names of essential columns
    /// in the genomic data [DataFrame] stored within the [Grangers] instance. It enables read-only operations
    /// to inspect column names and configurations.
    ///
    /// ### Returns
    ///
    /// Returns a reference to the [`FieldColumns`] struct.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let field_columns = grangers.field_columns();
    /// println!("Current field columns: {:?}", field_columns);
    /// ```
    pub fn field_columns(&self) -> &FieldColumns {
        &self.field_columns
    }

    /// Provides a mutable reference to the [FieldColumns] of the Grangers instance.
    ///
    /// This method allows for the modification of the [FieldColumns] struct, which specifies the names of essential columns
    /// in the genomic data [DataFrame] stored within the [Grangers] instance. This enables changing column names and configurations
    /// directly.
    ///
    /// ### Returns
    ///
    /// Returns a mutable reference to the [`FieldColumns`] struct.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let field_columns = grangers.field_columns_mut();
    /// field_columns.seqname = "chromosome".to_string();
    /// ```
    pub fn field_columns_mut(&mut self) -> &FieldColumns {
        &mut self.field_columns
    }

    /// Validates and provides a reference to the [FieldColumns] of the [Grangers] instance.
    ///
    /// This method validates the current [FieldColumns] against the [DataFrame] and provides a reference to it if valid.
    /// If the columns are invalid, depending on the `is_warn` and `is_bail` flags, it either warns the user or stops the execution.
    ///
    /// ### Arguments
    ///
    /// * `is_warn`: A boolean flag to indicate if a warning should be issued when validation fails.
    /// * `is_bail`: A boolean flag to indicate if an error should be returned when validation fails.
    ///
    /// ### Returns
    ///
    /// Returns a [`Result`] containing a reference to the [`FieldColumns`] struct or an error if validation fails and `is_bail` is `true`.
    ///
    /// ### Example
    ///
    /// ```rust
    /// match grangers.field_columns_checked(true, true) {
    ///     Ok(field_columns) => println!("Valid field columns: {:?}", field_columns),
    ///     Err(e) => println!("Error validating field columns: {}", e),
    /// }
    /// ```
    pub fn field_columns_checked(
        &self,
        is_warn: bool,
        is_bail: bool,
    ) -> anyhow::Result<&FieldColumns> {
        self.field_columns().is_valid(self.df(), is_warn, is_bail)?;
        Ok(self.field_columns())
    }

    /// Provides a reference to the [`DataFrame`] stored within the [Grangers] instance.
    ///
    /// This method returns a read-only reference to the genomic data [DataFrame] held by the [Grangers] instance.
    ///
    /// ### Returns
    ///
    /// Returns a reference to the [`DataFrame`].
    ///
    /// ### Example
    ///
    /// ```rust
    /// let dataframe = grangers.df();
    /// println!("Number of rows in DataFrame: {}", dataframe.height());
    /// ```
    pub fn df(&self) -> &DataFrame {
        &self.df
    }

    /// Provides a mutable reference to the [`DataFrame`] stored within the [Grangers] instance.
    ///
    /// This method allows for modifications to the genomic data [DataFrame] held by the [Grangers] instance.
    ///
    /// ### Returns
    ///
    /// Returns a mutable reference to the [`DataFrame`].
    ///
    /// ### Example
    ///
    /// ```rust
    /// let dataframe_mut = grangers.df_mut();
    /// dataframe_mut.sort_in_place("seqname", Default::default());
    /// ```
    pub fn df_mut(&mut self) -> &mut DataFrame {
        &mut self.df
    }

    /// Provides a reference to the [`IntervalType`] used by the Grangers instance.
    ///
    /// This method returns a read-only reference to the [`IntervalType`], which defines how genomic intervals are interpreted
    /// within the [Grangers] instance.
    ///
    /// ### Returns
    ///
    /// Returns a reference to the [`IntervalType`].
    ///
    /// ### Example
    ///
    /// ```rust
    /// let interval_type = grangers.interval_type();
    /// println!("Current interval type: {:?}", interval_type);
    /// ```
    pub fn interval_type(&self) -> &IntervalType {
        &self.interval_type
    }

    /// Provides a reference to the optional [SeqInfo] stored within the [Grangers] instance.
    ///
    /// This method returns an optional shared reference to the [SeqInfo] struct, which contains reference genome information.
    ///
    /// ### Returns
    ///
    /// Returns an optional reference to the [SeqInfo].
    ///
    /// ### Example
    ///
    /// ```rust
    /// if let Some(seqinfo) = grangers.seqinfo() {
    ///     println!("Reference genome information available.");
    /// } else {
    ///     println!("No reference genome information available.");
    /// }
    /// ```
    pub fn seqinfo(&self) -> Option<&SeqInfo> {
        self.seqinfo.as_ref()
    }

    /// Provides a mutable reference to the optional [SeqInfo] stored within the [Grangers] instance.
    ///
    /// This method allows for modifications to the [SeqInfo] struct, which contains reference genome information.
    ///
    /// ### Returns
    ///
    /// Returns an optional mutable reference to the [SeqInfo].
    ///
    /// ### Example
    ///
    /// ```rust
    /// if let Some(seqinfo_mut) = grangers.seqinfo_mut() {
    ///     seqinfo_mut.set_seqnames(vec!["chr1".to_string(), "chr2".to_string()]);
    /// }
    /// ```
    pub fn seqinfo_mut(&mut self) -> Option<&mut SeqInfo> {
        self.seqinfo.as_mut()
    }

    /// Sorts the [DataFrame] within the Grangers instance based on specified columns.
    ///
    /// This method sorts the internal [DataFrame] by the columns provided in the `by` argument.
    /// You can specify sorting order via the `descending` argument and whether to maintain the original order in case of ties with `maintain_order`.
    ///
    /// ### Arguments
    ///
    /// * `by`: A slice of strings representing the names of columns to sort by.
    /// * `descending`: An implementation of [`IntoVec<bool>`] that indicates whether sorting should be in descending order for each column.
    /// * `maintain_order`: A boolean that, if set to true, maintains the original order of rows in case of ties.
    ///
    /// ### Returns
    ///
    /// Returns an [anyhow::Result<()>](anyhow::Result) indicating the success or failure of the sorting operation.
    ///
    /// ### Example
    ///
    /// ```rust
    /// grangers.sort_by(&["gene_id", "start"], vec![false, true], false)?;
    /// ```
    pub fn sort_by<T>(
        &mut self,
        by: &[&str],
        descending: impl IntoIterator<Item = bool>,
        maintain_order: bool,
        multithreaded: bool,
    ) -> anyhow::Result<()> {
        self.df = self.df.sort(
            by.iter()
                .map(|s| s.to_owned().into())
                .collect::<Vec<PlSmallStr>>(),
            SortMultipleOptions::default()
                .with_order_descending_multi(descending)
                .with_maintain_order(maintain_order)
                .with_multithreaded(multithreaded),
        )?;
        self.inc_signature();
        Ok(())
    }

    /// Filters the [DataFrame] within the [Grangers] instance based on specified values in a column.
    ///
    /// This method creates a new [Grangers] instance containing rows from the internal [DataFrame] where values in the specified column match any of the provided values.
    ///
    /// ### Arguments
    ///
    /// * `by`: The name of the column to filter by.
    /// * `values`: A slice of values to filter by.
    /// * `warn_empty`: A boolean flag indicating whether to issue a warning if the filtered dataframe is empty.
    ///
    /// ### Returns
    ///
    /// Returns a new [Grangers] instance containing only the filtered rows.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let filtered_grangers = grangers.filter("gene_type", &["protein_coding", "lncRNA"])?;
    /// ```
    pub fn filter<T: AsRef<str>>(
        &self,
        by: T,
        values: &[T],
        warn_empty: bool,
    ) -> anyhow::Result<Grangers> {
        let column = self.get_column_name(by.as_ref(), false)?;
        // Create a list Series for is_in
        let values_vec: Vec<String> = values.iter().map(|s| s.as_ref().to_string()).collect();
        let list_series = Series::new("values".into(), vec![Series::new("inner".into(), values_vec)])
            .cast(&DataType::List(Box::new(DataType::String)))?;
        let s = self.df().column(&column)?.as_materialized_series();
        let df = self.df().filter(&is_in(s, &list_series, false)?)?;

        if df.is_empty() && warn_empty {
            warn!("The filtered dataframe is empty.")
        }
        Grangers::new(
            df,
            self.seqinfo().cloned(),
            self.misc.clone(),
            IntervalType::default(),
            self.field_columns().clone(),
            true,
        )
    }

    /// Retrieves the signature of the [Grangers] instance.
    ///
    /// This method returns the unique signature of the [Grangers] instance, which is useful for tracking changes or versions of the data.
    ///
    /// ### Returns
    ///
    /// Returns a 64-bit unsigned integer (`u64`) representing the signature of the instance.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let signature = grangers.get_signature();
    /// println!("Grangers instance signature: {}", signature);
    /// ```
    pub fn get_signature(&self) -> u64 {
        self.signature
    }

    /// Sets the signature of the [Grangers] instance.
    ///
    /// This method allows setting a new signature for the [Grangers] instance. This can be useful for manual versioning or tracking specific changes.
    ///
    /// ### Arguments
    ///
    /// * `other_sig`: The new signature to set for the [Grangers] instance.
    ///
    /// ### Example
    ///
    /// ```rust
    /// grangers.set_signature(new_signature);
    /// ```
    fn set_signature(&mut self, other_sig: u64) {
        self.signature = other_sig
    }

    /// Updates a specific column in the [Grangers] instance's [DataFrame].
    ///
    /// This method updates the internal [DataFrame] by replacing or adding the specified column.
    /// It can also update the internal mapping of field columns if a field column name is provided.
    ///
    /// ### Arguments
    ///
    /// * `column`: The new [Series] to be inserted or used to replace an existing column in the [DataFrame].
    /// * `field_column`: An optional string reference indicating the field column to be updated with the new column's name.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<()>`] indicating the success or failure of the update operation.
    ///
    /// ### Example
    ///
    /// ```rust
    /// grangers.update_column(new_series, Some("new_column"))?;
    /// ```
    pub fn update_column(
        &mut self,
        column: Column,
        field_column: Option<&str>,
    ) -> anyhow::Result<()> {
        // first we warn if there are null values in the column
        if column.null_count() > 0 {
            warn!("The provided Series object contains {} null values. This might cause problems when calling Grangers methods.", column.null_count());
        }

        // if a field_column is provided, we update the field_columns object
        if let Some(field_column) = field_column {
            self.field_columns
                .update(PlSmallStr::from(field_column), column.name().clone())?;
        }

        let name = column.name().to_owned();
        self.df.with_column(column).with_context(|| {
            format!(
                "Could not update Grangers with the provided Series object: {:?}",
                name
            )
        })?;

        // we don't want to do validation here because it might
        // complain about some existing nulls before the update
        // self.validate(false, false)?;
        self.inc_signature();
        Ok(())
    }

    /// Updates the [`DataFrame`] of the [`Grangers`] instance.
    ///
    /// This method allows replacing the current [`DataFrame`] with a new one. It checks for compatibility in terms of shape and column names.
    ///
    /// ### Arguments
    ///
    /// * `df`: The new [`DataFrame`] to set.
    /// * `is_warn`: A boolean flag to indicate whether to issue warnings for any discrepancies found.
    /// * `is_bail`: A boolean flag to indicate whether to halt the operation if discrepancies are found.
    ///
    /// ### Returns
    ///
    /// Returns an `anyhow::Result<()>` indicating the success or failure of the update operation.
    ///
    /// ### Example
    ///
    /// ```rust
    /// grangers.update_df(new_df, true, true)?;
    /// ```
    pub fn update_df(&mut self, df: DataFrame, is_warn: bool, is_bail: bool) -> anyhow::Result<()> {
        // check if the dataframe has the same layout as the current one
        if df.shape() != self.df.shape() {
            bail!("The provided dataframe has a different layout as the current one. Please use Grangers::new() to instantiate a new Grangers struct.")
        }

        // check if the dataframes have the same column nake
        let self_columns = self.df.get_column_names();
        let new_columns = df.get_column_names();

        if !self_columns.iter().all(|item| new_columns.contains(item)) {
            bail!("The provided dataframe have different column names as the current one. Please use Grangers::new() to instantiate a new Grangers struct.")
        }

        self.df = df;
        self.validate(is_warn, is_bail)?;
        self.inc_signature();
        Ok(())
    }
}

impl Grangers {
    /// Retrieves the column name from the Grangers instance's [`DataFrame`] or [`FieldColumns`].
    ///
    /// This method checks whether the provided name corresponds to a column in the [`DataFrame`] or an entry in the [`FieldColumns`].
    /// It validates against null values if `bail_null` is set to true. If the name is not a recognized column or field, the method will return an error.
    ///
    /// ### Generics
    ///
    /// * `T`: A type that implements [`AsRef<str>`], allowing for flexible string input.
    ///
    /// ### Arguments
    ///
    /// * `name`: The name of the column to retrieve. This can be either a direct column name or a name defined in FieldColumns.
    /// * `bail_null`: If true, the method will bail (return an error) if the column contains null values.
    ///
    /// ### Returns
    ///
    /// Returns a reference to the column name as a string slice (`&str`) if the column exists and meets the null-check condition.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let column_name = grangers.get_column_name_str("gene_id", true)?;
    /// println!("Column name: {}", column_name);
    /// ```
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
        } else if self
            .df()
            .get_column_names()
            .contains(&&PlSmallStr::from(name))
        {
            self.df().column(name)?.name()
        } else {
            bail!("{} is neither a column in the dataframe nor a field of FieldColumns. Cannot proceed", name)
        };

        if bail_null && self.df().column(name)?.null_count() > 0 {
            bail!("The column {} contains null values. Cannot proceed.", name)
        }

        Ok(name)
    }

    /// Retrieves the column name from the Grangers instance's [`DataFrame`] or [`FieldColumns`].
    ///
    /// Similar to `get_column_name_str`, this method returns the column name as a [String]. It performs checks to ensure the column exists
    /// and optionally verifies the absence of null values depending on `bail_null`.
    ///
    /// ### Generics
    ///
    /// * `T`: A type that implements [`AsRef<str>`], enabling various types to be used as string input.
    ///
    /// ### Arguments
    ///
    /// * `name`: The name of the column to find. This can refer to either an actual column name or a key in FieldColumns.
    /// * `bail_null`: A boolean flag that, when set to true, causes the method to return an error if the column contains null values.
    ///
    /// ### Returns
    ///
    /// Returns the name of the column as a [String] if the column is found and passes the null check if applicable.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let column_name = grangers.get_column_name("expression_level", false)?;
    /// println!("Column name: {}", column_name);
    /// ```
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
        } else if self
            .df()
            .get_column_names()
            .contains(&&PlSmallStr::from(name))
        {
            self.df().column(name)?.name()
        } else {
            bail!("{} is neither a column in the dataframe nor a field of FieldColumns. Cannot proceed", name)
        };

        if bail_null && self.df().column(name)?.null_count() > 0 {
            bail!("The column {} contains null values. Cannot proceed.", name)
        }

        Ok(name.to_string())
    }

    /// Retrieves a reference to a Series from the Grangers instance's [`DataFrame`] based on the column name.
    ///
    /// This method attempts to find a Series in the [`DataFrame`] directly or through the [`FieldColumns`] mapping.
    /// If the column is not found directly, it checks FieldColumns for a corresponding entry.
    ///
    /// ### Generics
    ///
    /// * `T`: A type that supports [`AsRef<str>`], enabling the use of various string types as input.
    ///
    /// ### Arguments
    ///
    /// * `name`: The name of the column to retrieve. This could be an actual [`DataFrame`] column name or a key from FieldColumns.
    ///
    /// ### Returns
    ///
    /// Returns a reference to the Series (&[`Series`]) corresponding to the specified column name if found.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let series = grangers.column("total_reads")?;
    /// println!("Total reads: {:?}", series);
    /// ```
    pub fn column<T: AsRef<str>>(&self, name: T) -> anyhow::Result<&Column> {
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

    /// Retrieves a list of Series from the Grangers instance's [`DataFrame`] for the specified columns.
    ///
    /// This method collects references to [Series] for each column name provided in the `names` array.
    /// It returns an error if any specified column name does not exist or if there is a problem retrieving the Series.
    ///
    /// ### Generics
    ///
    /// * `T`: A type that implements [`AsRef<str>`], allowing for different string input types.
    ///
    /// ### Arguments
    ///
    /// * `names`: An array of names (as references to strings) corresponding to the columns to retrieve.
    ///
    /// ### Returns
    ///
    /// Returns a vector of references to the Series objects for the requested column names.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let columns = grangers.columns(&["gene_id", "expression_level"])?;
    /// println!("Retrieved columns: {:?}", columns);
    /// ```
    pub fn columns<T: AsRef<str>>(&self, names: &[T]) -> anyhow::Result<Vec<&Column>> {
        let mut cols = Vec::new();
        for name in names {
            cols.push(self.column(name)?);
        }
        Ok(cols)
    }

    /// Retrieves the 'seqname' Series from the Grangers instance's [`DataFrame`].
    ///
    /// This method accesses the 'seqname' column defined in the FieldColumns of the Grangers instance,
    /// returning an error if the column does not exist or cannot be retrieved.
    ///
    /// ### Returns
    ///
    /// Returns a reference to the 'seqname' Series from the [`DataFrame`].
    ///
    /// ### Example
    ///
    /// ```rust
    /// let seqname_series = grangers.seqname()?;
    /// println!("Seqname Series: {:?}", seqname_series);
    /// ```
    pub fn seqname(&self) -> anyhow::Result<&Column> {
        self.column(self.field_columns.seqname())
    }

    /// Retrieves the 'start' Series from the Grangers instance's [`DataFrame`].
    ///
    /// This method accesses the 'start' column defined in the [FieldColumns] of the [Grangers] instance,
    /// returning an error if the column does not exist or cannot be retrieved.
    ///
    /// ### Returns
    ///
    /// Returns a reference to the 'start' Series from the [`DataFrame`].
    ///
    /// ### Example
    ///
    /// ```rust
    /// let start_series = grangers.start()?;
    /// println!("Start Series: {:?}", start_series);
    /// ```
    pub fn start(&self) -> anyhow::Result<&Column> {
        self.column(self.field_columns.start())
    }

    /// Retrieves the 'end' Series from the Grangers instance's [`DataFrame`].
    ///
    /// This method accesses the 'end' column defined in the [FieldColumns] of the [Grangers] instance,
    /// returning an error if the column does not exist or cannot be retrieved.
    ///
    /// ### Returns
    ///
    /// Returns a reference to the 'end' Series from the [`DataFrame`].
    ///
    /// ### Example
    ///
    /// ```rust
    /// let end_series = grangers.end()?;
    /// println!("End Series: {:?}", end_series);
    /// ```
    pub fn end(&self) -> anyhow::Result<&Column> {
        self.column(self.field_columns.end())
    }

    /// Retrieves the 'strand' Series from the Grangers instance's [`DataFrame`].
    ///
    /// This method accesses the 'strand' column defined in the [FieldColumns] of the [Grangers] instance,
    /// returning an error if the column does not exist or cannot be retrieved.
    ///
    /// ### Returns
    ///
    /// Returns a reference to the 'strand' [`Series`] from the [`DataFrame`].
    ///
    /// ### Example
    ///
    /// ```rust
    /// let strand_series = grangers.strand()?;
    /// println!("Strand Series: {:?}", strand_series);
    /// ```
    pub fn strand(&self) -> anyhow::Result<&Column> {
        self.column(self.field_columns.strand())
    }

    /// Retrieves the 'score' [`Series`] from the [Grangers] instance's [DataFrame].
    ///
    /// This method accesses the 'score' column defined in the [FieldColumns] of the [Grangers] instance.
    /// It returns an error if the 'score' column does not exist or cannot be retrieved from the [DataFrame].
    ///
    /// ### Returns
    ///
    /// Returns a reference to the 'score' [Series] from the [DataFrame].
    ///
    /// ### Example
    ///
    /// ```rust
    /// let score_series = grangers.score()?;
    /// println!("Score Series: {:?}", score_series);
    /// ```
    pub fn score(&self) -> anyhow::Result<&Column> {
        self.column(
            self.field_columns
                .score()
                .with_context(|| "Could not get the score column from the dataframe.")?,
        )
    }

    /// Retrieves the 'phase' Series from the [Grangers] instance's [DataFrame].
    ///
    /// This method accesses the 'phase' column defined in the FieldColumns of the Grangers instance.
    /// It returns an error if the 'phase' column does not exist or cannot be retrieved from the [DataFrame].
    ///
    /// ### Returns
    ///
    /// Returns a reference to the 'phase' Series from the [DataFrame].
    ///
    /// ### Example
    ///
    /// ```rust
    /// let phase_series = grangers.phase()?;
    /// println!("Phase Series: {:?}", phase_series);
    /// ```
    pub fn phase(&self) -> anyhow::Result<&Column> {
        self.column(
            self.field_columns
                .phase()
                .with_context(|| "Could not get the score column from the dataframe.")?,
        )
    }

    /// Retrieves the 'feature_type' [Series] from the [Grangers] instance's [DataFrame].
    ///
    /// This method accesses the 'feature_type' column defined in the [FieldColumns] of the [Grangers] instance.
    /// It returns an error if the 'feature_type' column does not exist or cannot be retrieved from the [DataFrame].
    ///
    /// ### Returns
    ///
    /// Returns a reference to the 'feature_type' [Series] from the [DataFrame].
    ///
    /// ### Example
    ///
    /// ```rust
    /// let feature_type_series = grangers.feature_type()?;
    /// println!("Feature Type Series: {:?}", feature_type_series);
    /// ```
    pub fn feature_type(&self) -> anyhow::Result<&Column> {
        self.column(
            self.field_columns
                .feature_type()
                .with_context(|| "Could not get the score column from the dataframe.")?,
        )
    }

    /// Retrieves a [DataFrame] representing the genomic range from the [Grangers] instance.
    ///
    /// This method selects the 'start', 'end', and 'strand' columns from the [Grangers] instance's [DataFrame]
    /// to construct a new [DataFrame] that represents the range of genomic features.
    ///
    /// ### Returns
    ///
    /// Returns a new [DataFrame] containing only the 'start', 'end', and 'strand' columns.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let range_df = grangers.range()?;
    /// println!("Range DataFrame: {:?}", range_df);
    /// ```
    pub fn range(&self) -> anyhow::Result<DataFrame> {
        let range = self.df.select([
            self.field_columns().start(),
            self.field_columns().end(),
            self.field_columns().strand(),
        ])?;
        Ok(range)
    }

    /// Checks if a column exists in the [Grangers] instance's [DataFrame].
    ///
    /// This method verifies whether a column with the specified name exists within the [DataFrame].
    ///
    /// ### Generics
    ///
    /// * `T`: A type that implements [`AsRef<str>`], allowing for different string input types.
    ///
    /// ### Arguments
    ///
    /// * `name`: The name of the column to check for existence.
    ///
    /// ### Returns
    ///
    /// Returns `true` if the column exists, otherwise returns `false`.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let exists = grangers.is_column("gene_id");
    /// println!("Does the 'gene_id' column exist? {}", exists);
    /// ```
    pub fn is_column<T: AsRef<str>>(&self, name: T) -> bool {
        self.column(name).is_ok()
    }
}

// validate Grangers
impl Grangers {
    /// Validates the [Grangers] instance's [DataFrame] and field columns.
    ///
    /// This method performs several checks: it ensures the [DataFrame] is not empty, validates the field columns, and checks for null values in essential fields.
    /// It issues warnings or errors based on the `is_warn` and `is_bail` flags.
    ///
    /// ### Arguments
    ///
    /// * `is_warn`: A boolean flag that, if set to true, will cause the method to issue warnings for validation issues.
    /// * `is_bail`: A boolean flag that, if set to true, will cause the method to return an error and halt execution if validation issues are detected.
    ///
    /// ### Returns
    ///
    /// Returns `Ok(true)` if the [DataFrame] and field columns pass all validation checks. Returns [Ok]`(false)` if there are validation issues but `is_bail` is set to false.
    /// Returns an error if `is_bail` is set to true and there are validation issues.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let is_valid = grangers.validate(true, true)?;
    /// println!("Is the Grangers instance valid? {}", is_valid);
    /// ```
    pub fn validate(&self, is_warn: bool, is_bail: bool) -> anyhow::Result<bool> {
        if self.df().height() == 0 {
            if is_bail {
                bail!("The dataframe is empty. Cannot proceed.")
            } else {
                if is_warn {
                    warn!("The dataframe is empty.")
                }
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

    /// Fixes the field columns in the [Grangers] instance based on the current [DataFrame].
    ///
    /// This method attempts to repair issues with the field columns, such as missing or incorrect column names, by adjusting them to match the current [DataFrame] structure.
    /// It will issue warnings if `is_warn` is set to true and if there are discrepancies between the field columns and the [DataFrame].
    ///
    /// ### Arguments
    ///
    /// * `is_warn`: A boolean flag that, if set to true, causes the method to issue warnings when discrepancies are found and fixed.
    ///
    /// ### Returns
    ///
    /// Returns [Ok]`(())` if the field columns were successfully fixed or if no issues were found. Returns an error if the field columns cannot be fixed.
    ///
    /// ### Example
    ///
    /// ```rust
    /// grangers.fix_field_columns(true)?;
    /// println!("Field columns have been fixed.");
    /// ```
    pub fn fix_field_columns(&mut self, is_warn: bool) -> anyhow::Result<()> {
        let mut field_columns = self.field_columns().clone();
        field_columns.fix(self.df(), is_warn)?;
        self.field_columns = field_columns;
        Ok(())
    }
}

// implement GenomicFeatures for Grangers
impl Grangers {
    /// Computes the intronic regions from the exon annotations in the [Grangers] instance.
    ///
    /// This method calculates the genomic regions corresponding to introns based on exon records. It allows customization of the aggregation level (e.g., by transcript or gene),
    /// the specific exon feature to filter by, and additional columns to keep in the resulting [Grangers] instance.
    ///
    /// ### Arguments
    ///
    /// * `by`: Optional reference to a string specifying the column by which to group exons before computing introns. Commonly set to "transcript_id" or "gene_id".
    /// * `exon_feature`: Optional reference to a string specifying the exon feature type to consider. If None, all exon features are considered.
    /// * `keep_columns`: Optional slice of string references specifying additional columns to keep in the output.
    /// * `multithreaded`: Boolean flag indicating whether to use multithreading for performance improvement.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<Grangers>`] containing the calculated intronic regions.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let introns = grangers.introns(Some("transcript_id"), None, None, false)?;
    /// ```
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

        exon_gr.gaps(&[by_str], false, None, Some(&kc), multithreaded)
    }

    /// Computes the boundary regions for genes from the exon annotations in the [Grangers] instance.
    ///
    /// This method calculates the genomic boundary regions for genes based on exon records. It allows filtering by a specific exon feature and determines the boundaries based on the "gene_id" column.
    ///
    /// ### Arguments
    ///
    /// * `exon_feature`: Optional reference to a string specifying the exon feature type to consider. If None, all exon features are considered.
    /// * `multithreaded`: Boolean flag indicating whether to use multithreading for performance improvement.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<Grangers>`] containing the calculated gene boundary regions.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let gene_boundaries = grangers.genes(None, false)?;
    /// ```
    pub fn genes(
        &self,
        exon_feature: Option<&str>,
        multithreaded: bool,
    ) -> anyhow::Result<Grangers> {
        self.validate(false, true)?;
        let gene_id = self.get_column_name_str("gene_id", true)?;
        self.boundary(gene_id, exon_feature, None, multithreaded)
    }

    /// Computes the boundary regions for transcripts from the exon annotations in the [Grangers] instance.
    ///
    /// This method calculates the genomic boundary regions for transcripts based on exon records. It allows filtering by a specific exon feature and determines the boundaries based on the "transcript_id" column.
    ///
    /// ### Arguments
    ///
    /// * `exon_feature`: Optional reference to a string specifying the exon feature type to consider. If [None], all exon features are considered.
    /// * `multithreaded`: Boolean flag indicating whether to use multithreading for performance improvement.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<Grangers>`] containing the calculated transcript boundary regions.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let transcript_boundaries = grangers.transcripts(None, false)?;
    /// ```
    pub fn transcripts(
        &self,
        exon_feature: Option<&str>,
        extra_group_by_columns: Option<&[&str]>,
        multithreaded: bool,
    ) -> anyhow::Result<Grangers> {
        let transcript_id = self.get_column_name_str("transcript_id", true)?;
        self.boundary(transcript_id, exon_feature, extra_group_by_columns, multithreaded)
    }

    /// Computes the boundary regions of genomic features (like genes or transcripts) based on their exons.
    ///
    /// This method calculates the start and end positions for each genomic feature by aggregating exon information.
    /// It groups exons by a specified field (usually 'gene_id' or 'transcript_id'), then calculates the minimum start
    /// and maximum end positions to determine the boundary of each feature.
    ///
    /// ### Arguments
    ///
    /// * `by`: The name of the field by which to group exons before calculating boundaries (e.g., 'gene_id' or 'transcript_id').
    /// * `exon_feature`: Optional reference to a string specifying the exon feature type to consider. If [None], all exon features are considered.
    /// * `multithreaded`: Boolean flag indicating whether to use multithreading for performance improvement.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<Grangers>`] containing the genomic features with their calculated boundaries.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let gene_boundaries = grangers.boundary("gene_id", None, false)?;
    /// ```
    pub fn boundary(
        &self,
        by: &str,
        exon_feature: Option<&str>,
        extra_group_by_columns: Option<&[&str]>,
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

        let mut group_by_columns = vec![seqname, by, strand];
        if let Some(extra_group_by_columns) = extra_group_by_columns {
            for c in extra_group_by_columns {
                group_by_columns.push(self.get_column_name_str(c, true)?);
            }
        }

        // check if genes are well defined: all features of a gene should have a valid seqname, start, end, and strand
        let any_invalid = exon_gr
            .df()
            .select(group_by_columns.clone())?
            .lazy()
            .group_by([by])
            .agg([
                col(seqname)
                    .unique()
                    .count()
                    .neq(lit(1))
                    .alias("seqname_any"),
                col(strand).unique().count().neq(lit(1)).alias("strand_any"),
            ])
            .select([col("seqname_any").any(true), col("strand_any").any(true)]) // true: drop nulls
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
            .group_by(group_by_columns)
            .agg([col(start).min(), col(end).max()])
            .collect()?;

        exon_gr.fix_field_columns(false)?;

        Ok(exon_gr)
    }

    /// Filters the [Grangers] instance to only include exon features, optionally filtered by a specific feature type.
    ///
    /// This method reduces the current [Grangers] dataset to only include records classified as exons, optionally filtering
    /// them by a specific exon feature type. This is primarily used as a preparatory step for other analyses such as
    /// calculating introns or feature boundaries.
    ///
    /// ### Arguments
    ///
    /// * `exon_feature`: Optional reference to a string specifying the exon feature type to consider. If [None], all exon features are considered.
    /// * `multithreaded`: Boolean flag indicating whether to use multithreading for performance improvement.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<Grangers>`] containing only exon records, optionally filtered by the specified feature type.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let exons_only = grangers.exons(Some("coding"), false)?;
    /// ```
    pub fn exons(
        &self,
        exon_feature: Option<&str>,
        multithreaded: bool,
    ) -> anyhow::Result<Grangers> {
        // validate itself
        // essential fields have no null values
        // all fields correspond to a column in the dataframe
        self.validate(false, true)?;

        let exon_feature = exon_feature.unwrap_or("exon");

        // feature_type can have null values, they will be ignored
        let feature_type = self.get_column_name_str("feature_type", false)?;

        if self.column(feature_type)?.null_count() > 0 {
            warn!("Found rows with a null `{}` value. These rows will be ignored when selecting exon records.", feature_type)
        }

        // polars way to subset
        let mut exon_gr = self.filter(feature_type, &[exon_feature], true)?;

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
        let valid_strands_vec: Vec<String> = VALIDSTRANDS.iter().map(|s| s.to_string()).collect();
        let list_series = Series::new("values".into(), vec![Series::new("inner".into(), valid_strands_vec)])
            .cast(&DataType::List(Box::new(DataType::String)))?;
        if !is_in(&exon_gr.column(strand)?.as_materialized_series().unique()?, &list_series, false)?
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
            .group_by([seqname, transcript_id])
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
                    col(seqname).cast(DataType::Categorical(None, CategoricalOrdering::Lexical)),
                    col(strand).cast(DataType::Categorical(None, CategoricalOrdering::Lexical)),
                    col(transcript_id)
                        .cast(DataType::Categorical(None, CategoricalOrdering::Lexical)),
                    col(exon_number.as_str()),
                ],
                SortMultipleOptions::default().with_multithreaded(multithreaded),
            )])
            .collect()?;

        // well done!
        fc.fix(exon_gr.df(), false)?;
        exon_gr.field_columns = fc;
        Ok(exon_gr)
    }

    /// Extends genomic intervals in the dataframe based on specified options and strand information.
    ///
    /// This method modifies the start and/or end positions of intervals in the dataframe by a specified length.
    /// The extension can be applied to the start, end, or both sides of each interval, and can be strand-specific or not.
    ///
    /// ### Arguments
    ///
    /// * `length`: The length by which to extend the intervals. This value can be positive or negative.
    /// * `extend_option`: Specifies which end(s) of the intervals should be extended (`Start`, `End`, or `Both`).
    /// * `ignore_strand`: If `true`, extends intervals without considering strand information; if `false`, extensions are strand-specific.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<()>`](anyhow::Result). The function modifies the Grangers instance in place and does not return a value.
    ///
    /// ### Example
    ///
    /// ```rust
    /// grangers.extend(1000, &ExtendOption::Both, false)?;
    /// ```
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
                | {
                    let valid_strands_vec: Vec<String> = VALIDSTRANDS.iter().map(|s| s.to_string()).collect();
                    let list_series = Series::new("values".into(), vec![Series::new("inner".into(), valid_strands_vec)])
                        .cast(&DataType::List(Box::new(DataType::String)))?;
                    !is_in(&self.column(strand)?.as_materialized_series().unique()?, &list_series, false)?
                        .all()
                }
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

    /// Generates flanking regions for genomic intervals in the dataframe.
    ///
    /// This method calculates new genomic intervals representing flanking regions. The width and side of the flank
    /// (upstream or downstream) can be specified, as well as whether to consider both sides or ignore strand information.
    ///
    /// ### Arguments
    ///
    /// * `width`: The width of the flanking regions. Positive values generate regions upstream; negative values generate downstream regions.
    /// * `options`: A [FlankOptions] struct that specifies how the flanking regions should be determined, including whether to ignore strand information, and whether to generate flanks on both sides.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<Grangers>`] containing the genomic intervals of the flanking regions.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let flank_options = FlankOptions { start: true, both: false, ignore_strand: false };
    /// let flanking_regions = grangers.flank(500, flank_options)?;
    /// ```
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
                    .map(|&x| col(x.to_owned()))
                    .collect::<Vec<Expr>>(),
            )
            .collect()?;
        Grangers::new(
            df,
            self.seqinfo.clone(),
            self.misc.clone(),
            self.interval_type,
            self.field_columns.clone(),
            false,
        )
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
        todo!("not yet implemented");
    }

    /// this function turns the seqinfo of the Grangers object into a boundary Grangers object.
    pub fn seqinfo_to_bounary(&self) -> anyhow::Result<Grangers> {
        todo!("not yet implemented");
    }

    /// Calculates gaps between genomic intervals grouped by specified columns.
    ///
    /// This method identifies gaps between intervals in the dataframe when they are grouped by specified keys.
    /// The resulting dataframe will only include these gap intervals. This can be used, for example, to find
    /// intergenic regions when the input is exon intervals.
    ///
    /// ### Generics
    ///
    /// * `T`: A type that implements [`AsRef<str>`], allowing for flexible string references.
    ///
    /// ### Arguments
    ///
    /// * `by`: Columns by which to group intervals before identifying gaps.
    /// * `ignore_strand`: If `true`, ignores strand information when identifying gaps.
    /// * `slack`: Optional slack size to reduce the size of gaps; useful for filtering out small gaps.
    /// * `keep_columns`: Optional additional columns to keep in the output.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<Grangers>`] containing a new [Grangers] instance with the identified gaps.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let gaps = grangers.gaps(&["gene_id"], false, None, Some(&["seqname", "source"]))?;
    /// ```
    pub fn gaps<T: AsRef<str>>(
        &self,
        by: &[T],
        ignore_strand: bool,
        slack: Option<usize>,
        keep_columns: Option<&[&str]>,
        multithreaded: bool,
    ) -> anyhow::Result<Grangers> {
        // merge returns a sorted and merged Grangers object
        let mut gr = self.merge(by, ignore_strand, slack, keep_columns, multithreaded)?;

        gr.df = gr.apply(
            by,
            None,
            ignore_strand,
            apply_gaps,
            keep_columns,
            multithreaded,
        )?;
        Ok(gr)
    }

    /// Adds an order column to the dataframe based on the sorting of another column.
    ///
    /// This method sorts the dataframe based on a specified 'start' column and adds a new column indicating
    /// the order. This can be used, for example, to assign exon numbers within transcripts.
    ///
    /// ### Arguments
    ///
    /// * `by`: Optional columns by which to group data before ordering. If [None], the entire dataframe is ordered.
    /// * `name`: Name of the new order column to be added.
    /// * `offset`: Optional starting value for the order (default is 1). When adding exon_number for example, the offset should be 1.
    /// * `multithreaded`: If `true`, sorting is done in parallel, improving performance on large datasets.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<()>`](anyhow::Result). The method modifies the Grangers instance in place.
    ///
    /// ### Example
    ///
    /// ```rust
    /// grangers.add_order(Some(&["gene_id", "transcript_id"]), "exon_number", None, true)?;
    /// ```
    pub fn add_order(
        &mut self,
        by: Option<&[&str]>,
        name: &str,
        offset: Option<u32>,
        multithreaded: bool,
    ) -> anyhow::Result<()> {
        self.validate(false, true)?;

        // we make the default offset 1
        let offset = offset.unwrap_or(1);

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
                                .arg_sort(SortOptions::default().with_multithreaded(multithreaded)),
                        )
                        .otherwise(
                            col(start).arg_sort(
                                SortOptions::default()
                                    .with_order_descending(true)
                                    .with_multithreaded(multithreaded),
                            ),
                        )
                        .add(Expr::Literal(LiteralValue::Scalar(offset.into())))
                        .over(by)
                        .cast(DataType::String)
                        .alias(name),
                )
                .collect()?;
        } else {
            self.df.with_row_index(name.into(), Some(offset))?;
        }
        Ok(())
    }

    /// Drops rows with null values in specified fields of the [Grangers] instance's dataframe.
    ///
    /// This method allows selective removal of rows where null values occur in specified fields. If no fields
    /// are specified, it removes rows where any null values are present.
    ///
    /// ### Arguments
    ///
    /// * `fields`: Optional array of field names to check for null values. If [None], all fields are checked.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<()>`](anyhow::Result). The method modifies the [Grangers] instance in place by removing rows
    /// with null values in the specified fields.
    ///
    /// ### Example
    ///
    /// ```rust
    /// grangers.drop_nulls(Some(&["seqname", "start", "end"]))?;
    /// ```
    pub fn drop_nulls(&mut self, fields: Option<&[&str]>) -> anyhow::Result<()> {
        self.validate(false, true)?;
        // check the validity of the column names
        let cols: Vec<String> = match fields {
            Some(names) => self
                .columns(names)?
                .iter()
                .map(|s| s.name().to_owned().into_string())
                .collect(),
            None => self
                .df
                .get_column_names()
                .into_iter()
                .map(|s| s.to_owned().into_string())
                .collect(),
        };

        *self.df_mut() = self.df().drop_nulls(Some(&cols))?;
        Ok(())
    }

    /// Merges overlapping or adjacent genomic intervals in the [Grangers] instance.
    ///
    /// This method merges intervals that are either overlapping or adjacent within specified groupings.
    /// The merge process can consider strand information and include a slack region between intervals.
    /// Optionally, specific columns can be retained in the merged result.
    ///
    /// ### Generics
    ///
    /// * `T`: A type that implements [`AsRef<str>`], allowing for flexible string references.
    ///
    /// ### Arguments
    ///
    /// * `by`: Columns by which to group intervals before merging.
    /// * `ignore_strand`: If `true`, ignores strand information during the merge process.
    /// * `slack`: Optional slack size allows for merging intervals that are within a certain distance apart.
    /// * `keep_columns`: Optional additional columns to keep in the output.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<Grangers>`] containing a new [Grangers] instance with merged intervals.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let merged_grangers = grangers.merge(&["gene_id"], false, Some(10), Some(&["seqname", "source"]))?;
    /// ```
    pub fn merge<T: AsRef<str>>(
        &self,
        by: &[T],
        ignore_strand: bool,
        slack: Option<usize>,
        keep_columns: Option<&[&str]>,
        multithreaded: bool,
    ) -> anyhow::Result<Grangers> {
        self.validate(false, true)?;
        let df = self.apply(
            by,
            slack,
            ignore_strand,
            apply_merge,
            keep_columns,
            multithreaded,
        )?;

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
        multithreaded: bool,
    ) -> anyhow::Result<DataFrame>
    where
        F: Fn(Column, i64) -> Result<Option<polars::prelude::Column>, PolarsError>
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
            .filter(|&&n| !sorted_by_exprs_essential.contains(&col(n)))
            .map(|&n| col(n))
            .collect();

        let mut sorted_by_desc = vec![false; sorted_by_exprs.len()];
        sorted_by_exprs.extend(sorted_by_exprs_essential);
        sorted_by_desc.extend(sorted_by_desc_essential);

        // the lazy API of polars takes the ownership of a dataframe
        // we want to keep the keep_columns
        if let Some(keep_columns) = keep_columns {
            for &c in keep_columns {
                if !selected.contains(&self.get_column_name_str(c, false)?) {
                    selected.push(c);
                }
            }
        }

        let mut df = df.select(selected)?;

        // we will do the following
        // 1. sort the dataframe by the `by` columns + start and end columns
        // 2. group by the `by` columns
        // 3. build the new start and end columns by applying the `apply_fn` function
        // 4. explode the new start and end columns (makes the columns tall instead of wide)
        // 5. drop the intermediate columns
        df = df
            .lazy()
            // TODO: This can be replaced by select([all().sort(essentials).over(groups)]). Not sure if it is faster
            .sort_by_exprs(
                &sorted_by_exprs,
                SortMultipleOptions::default()
                    .with_order_descending_multi(sorted_by_desc.clone())
                    .with_multithreaded(multithreaded),
                // &sorted_by_desc,
                // false, /*nulls last*/
                // false, /*force stable sort*/
            )
            .group_by(by.iter().map(|&s| col(s)).collect::<Vec<Expr>>())
            .agg([
                all().exclude([start, end]).first(),
                // process two columns at once
                // Notice the df is sorted
                as_struct([col(start), col(end)].to_vec())
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
                    .get(lit(0), false)
                    .alias(start),
                col("start_end_list-temp-nobody-will-use-this-name-right")
                    .list()
                    .get(lit(1), false)
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
                    .cast(DataType::String)
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
            .sort_by_exprs(
                sorted_by_exprs,
                SortMultipleOptions::default()
                    .with_order_descending_multi(sorted_by_desc.clone())
                    .with_maintain_order(multithreaded),
                // sorted_by_desc,
                // false, /*nulls last*/
                // false, /*force stable sort*/
            )
            .collect()?;

        Ok(df)
    }
}

// lappers
impl Grangers {
    /// Constructs lapper data structures (efficient interval trees) for fast interval queries on genomic data.
    ///
    /// This method builds a set of lappers, which are efficient data structures for interval overlap queries,
    /// based on the genomic intervals represented in the [Grangers] instance. It can optionally ignore invalid
    /// intervals and can group intervals by specified attributes before building the lappers.
    ///
    /// ### Generics
    ///
    /// * `T`: A type that implements [`AsRef<str>`], enabling flexible string references for grouping criteria.
    ///
    /// ### Arguments
    ///
    /// * `ignore_invalid`: If `true`, invalid intervals (e.g., negative start positions) are ignored instead of causing an error.
    /// * `ignore_strand`: If `true`, strand information is not considered when building lappers, which can be useful for unstranded data.
    /// * `group_by`: Columns used to group intervals before constructing individual lappers. Each group results in a separate lapper.
    ///
    /// ### Returns
    ///
    /// Returns an `anyhow::Result<HashMap<[String; 2], LapperType>>` where each key is a pair of strings (typically representing sequence name and strand)
    /// and each value is a Lapper data structure containing the grouped intervals. The LapperType is typically `Lapper<u64, (usize, Vec<String>)>`.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let lappers = grangers.build_lappers(false, false, &["gene_id"])?;
    /// ```
    pub fn build_lappers<T: AsRef<str>>(
        &mut self,
        ignore_invalid: bool,
        ignore_strand: bool,
        group_by: &[T],
    ) -> anyhow::Result<HashMap<[String; 2], LapperType>> {
        // rust-lapper
        let start_time = std::time::Instant::now();
        // validate the Grangers object
        self.validate(false, true)?;

        let mut by = Vec::new();
        for b in group_by.iter() {
            let name = b.as_ref();
            if !self.field_columns().gtf_fields().contains(&name) {
                warn!("The provided `by` vector contains a non-attribute column. There should be a strong reason of doing so")
            }
            by.push(self.get_column_name_str(name, true)?);
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
        let valid_rows_df = df
            .select([start, end, strand])?
            .lazy()
            .select([
                col(start).gt(lit(0)),
                col(end).gt(lit(0)),
                col(strand).eq(lit("+")).or(col(strand).eq(lit("-"))),
            ])
            .select([
                col(start).and(col(end)).alias("pos_valid"),
                col(start)
                    .and(col(end))
                    .and(col(strand))
                    .alias("pos_strand_valid"),
            ])
            .collect()?;
        let valid_pos = valid_rows_df.column("pos_valid")?.as_materialized_series();
        let valid_pos_strand = valid_rows_df
            .column("pos_strand_valid")?
            .as_materialized_series();

        // Then we bail if ignore_invalid is false but we found invalid features
        if !ignore_invalid {
            if ignore_strand {
                // we need to make sure start and end are valid
                if valid_pos.iter().any(|v| v != AnyValue::Boolean(true)) {
                    bail!("Found features with non-positive start/end position. Please remove them first or set ignore_invalid to true.")
                }
            } else {
                // we need to make sure start, end and strand are all valid
                if valid_pos_strand
                    .iter()
                    .any(|v| v != AnyValue::Boolean(true))
                {
                    bail!("Found features with non-positive start/end/strand position. Please remove them first or set ignore_invalid to true.")
                }
            }
        }

        // [start, stop) Inclusive start, exclusive of stop
        // we define start and end as u64, and we use Vec<String> to store group_by column values
        type Iv = Interval<u64, (usize, Vec<String>)>;

        let mut by_iters = df
            .columns(by)?
            .iter()
            .map(|s| s.as_materialized_series().iter())
            .collect::<Vec<_>>();

        let mut ess_iters = self
            .df()
            .columns(selected)?
            .iter()
            .map(|s| s.as_materialized_series().iter())
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
                println!("Found invalid row at row {}: {}", rid, is_valid);
                // we pop the invalid row
                ess_iters[0].next();
                ess_iters[1].next();
                ess_iters[2].next();
                ess_iters[3].next();
                continue;
            }

            // then, we take the start and end
            let s = ess_iters[0]
                .next()
                .expect("should have as many iterations as rows")
                .cast(&DataType::Int64)
                .try_extract::<i64>()? as u64;

            // we add 1 to the end because rust-lappers uses right-exclusive intervals
            let e = ess_iters[1]
                .next()
                .expect("should have as many iterations as rows")
                .cast(&DataType::Int64)
                .try_extract::<i64>()? as u64
                + 1;

            // we take the seqname and strand
            let seqn = if let AnyValue::String(t) = ess_iters[2]
                .next()
                .expect("should have as many iterations as rows")
            {
                t.to_string()
            } else {
                bail!("Could not get the seqname of the feature")
            };

            let strd = if ignore_strand {
                String::from(".")
            } else if let AnyValue::String(t) = ess_iters[3]
                .next()
                .expect("should have as many iterations as rows")
            {
                t.to_string()
            } else {
                bail!("Could not get the strand of the feature")
            };

            // we take the by columns
            let mut by_vec = Vec::new();
            for it in by_iters.iter_mut() {
                let v = if let AnyValue::String(t) =
                    it.next().expect("should have as many iterations as rows")
                {
                    t.to_string()
                } else {
                    bail!("Could not get the strand of the feature")
                };

                by_vec.push(v);
            }
            let lapper_tree_vec = lapper_tree_vec_hm
                .entry([seqn.clone(), strd.clone()])
                .or_insert(Vec::new());
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

        let duration: std::time::Duration = start_time.elapsed();
        debug!("build rust-lappers in {:?}", duration);
        Ok(lappers)
    }

    /// Find the features that overlap with the given interval in lapper.
    /// Notice that here the interval is inclusive, i.e., [start, end], which is the same type of interval used in Grangers but different with the right-exclusive interval type used in rust_lapper.
    /// Also note that as Grangers builds a lapper data structure for each (seqname, strand) pair (or each seqname if ignore_strand), you need to provide the seqname and (optioanl) strand of the interval.
    /// This function takes three parameters:
    /// 1. `start`: the (inclusive) start position of the interval
    /// 2. `end`: the (inclusive) end position of the interval
    /// 3. `seqname`: the seqname of the interval
    /// 4. `strand`: the strand of the interval. This argumenet should match the `ignore_strand` variable when building the lapper. If `ignore_strand` is true, this argument must be a Some variant. Otherwise, it must be a None variant.
    pub fn _lapper_find<T: AsRef<str>>(
        &self,
        _interval: InclusiveInterval,
        _seqname: T,
        _strand: Option<Strand>,
    ) -> anyhow::Result<()> {
        todo!("not yet implemented");
        // // we first check if the lappers have been built
        // let lappers = if let Some(lappers) = &self.lappers {
        //     lappers
        // } else {
        //     bail!("Could not find the lappers field. Please call build_lappers() first")
        // };

        // // Then we check if we can get the lapper for the given seqname and strand
        // // we get the ignore_strand used for building the lappers
        // let ignore_strand = if let Some(ignore_strand) = self.lappers_ignore_strand {
        //     ignore_strand
        // } else {
        //     bail!("Could not determine if strand is ignore while built lappers. Please rebuild the lappers by calling the build_lappers() method. If the lappers were built with the `build_lapers` function, this should not happen. Please report this bug on GitHub!")
        // };

        // // Then we check if ignore_strand matches the strand argument
        // if (ignore_strand & strand.is_some()) |
        //     ((!ignore_strand) & strand.is_none()) {
        //         bail!("The strand argument does not match the ignore_strand flag. If lappers were built with ignore_strand=true, the strand argument should be None. If lappers were built with ignore_strand=false, the strand argument should be Some.")
        //     }

        // // get a valid strand string
        // let strand = if let Some(strand) = strand {
        //     strand.to_string()
        // } else {
        //     String::from(".")
        // };

        // // Now, we can try to get the lapper from the hashmap
        // let lapper = if let Some(lapper) = lappers.get(&[seqname.as_ref().to_string(), strand]) {
        //     lapper
        // } else {
        //     bail!("Could not find the lapper for the given seqname and strand. Please make sure that the provided seqname is a valid seqname in the Grangers object and the strand argument is either \"Some(+)\" or \"-\"")
        // };

        // // Finally, we can query the lapper
        // // we need to add 1 to the end because rust-lappers uses right-exclusive intervals
        // let start = interval.start;
        // let end = interval.end + 1;

        // let overlaps = lapper.find(start, end);

        // for overlap in overlaps {
        //     println!("{:?}", overlap);
        // }

        // Ok(())
    }
}

// implement get sequence functions for Grangers
impl Grangers {
    /// Writes the sequences of transcripts to a FASTA format output file based on exon information.
    ///
    /// This method writes sequences of transcripts, which are constructed by concatenating exon sequences,
    /// to the specified output file. The exons are filtered based on the optional `exon_name` parameter.
    ///
    /// ### Generics
    ///
    /// * `T`: A type that implements [`AsRef<Path>`], allowing for flexible path references to the reference genome file.
    /// * `W`: A type that implements [Write], designating the output stream for writing the transcript sequences.
    ///
    /// ### Arguments
    ///
    /// * `ref_path`: The path to the reference genome FASTA file.
    /// * `out_file`: The output stream to which the transcript sequences will be written.
    /// * `exon_name`: Optional parameter specifying the name of the exon feature. If not provided, default exon feature names are used.
    /// * `multithreaded`: Boolean indicating whether the operation should be performed using multiple threads.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<()>`](anyhow::Result). If successful, transcript sequences are written to the output file; otherwise, an error is returned.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let mut grangers = Grangers::new(...);
    /// let ref_path = "reference.fasta";
    /// let out_file = File::create("transcript_sequences.fasta")?;
    /// grangers.write_transcript_sequences(&ref_path, out_file, None, true)?;
    /// ```
    pub fn write_transcript_sequences<T: AsRef<Path>, W: Write>(
        &mut self,
        ref_path: T,
        out_file: W,
        exon_name: Option<&str>,
        multithreaded: bool,
    ) -> anyhow::Result<()> {
        // let null_fn =
        self.write_transcript_sequences_with_filter(
            ref_path,
            out_file,
            exon_name,
            multithreaded,
            &mut None::<fn(&noodles::fasta::Record) -> bool>,
        )
    }

    /// Writes the sequences of transcripts to a FASTA format output file based on exon information,
    /// applying an optional filter to each transcript sequence before writing.
    ///
    /// This advanced method allows for filtering of transcript sequences based on custom criteria before writing.
    /// The filter is a closure or function pointer provided by the user that accepts a `&`[`noodles::fasta::Record`]
    /// and returns a boolean indicating whether to write the sequence to the output.
    ///
    /// ### Generics
    ///
    /// * `T`: A type that implements [`AsRef<Path>`], allowing for flexible path references to the reference genome file.
    /// * `W`: A type that implements [Write], designating the output stream for writing the transcript sequences.
    /// * `F`: A type that implements [FnMut]`(&`[noodles::fasta::Record]`) -> bool`, representing the filtering function.
    ///
    /// ### Arguments
    ///
    /// * `ref_path`: The path to the reference genome FASTA file.
    /// * `out_file`: The output stream to which the filtered transcript sequences will be written.
    /// * `exon_name`: Optional parameter specifying the name of the exon feature. If not provided, default exon feature names are used.
    /// * `multithreaded`: Boolean indicating whether the operation should be performed using multiple threads.
    /// * `record_filter`: Optional mutable reference to a filter function or closure to apply to each transcript sequence.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<()>`](anyhow::Result). If successful, filtered transcript sequences are written to the output file; otherwise, an error is returned.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let mut grangers = Grangers::new(...);
    /// let ref_path = "reference.fasta";
    /// let out_file = File::create("filtered_transcript_sequences.fasta")?;
    /// let mut filter = |record: &noodles::fasta::Record| record.sequence().len() > 100;
    /// grangers.write_transcript_sequences_with_filter(&ref_path, out_file, None, true, &mut Some(filter))?;
    /// ```
    pub fn write_transcript_sequences_with_filter<T: AsRef<Path>, W: Write, F>(
        &mut self,
        ref_path: T,
        out_file: W,
        exon_name: Option<&str>,
        multithreaded: bool,
        record_filter: &mut Option<F>,
    ) -> anyhow::Result<()>
    where
        F: FnMut(&noodles::fasta::Record) -> bool,
    {
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
        let all_seqnames = exon_gr
            .seqname()?
            .unique()?
            .str()?
            .into_iter()
            .map(|s| s.unwrap().to_string())
            .collect::<HashSet<_>>();

        // we get all seqnames
        // let all_seqnames = HashSet::from_iter(exon_gr.seqname()?.unique()?.into_iter());

        let mut reader = grangers_utils::get_noodles_reader_from_path(ref_path)?;
        // we also create a fasta writer
        let out_writer = BufWriter::with_capacity(4194304, out_file);
        let mut writer = noodles::fasta::io::writer::Builder::default()
            .set_line_base_count(usize::MAX)
            .build_from_writer(out_writer);

        // we iterate the fasta reader. For each fasta reacord (usually chromosome), we do
        // 1. subset the dataframe by the seqname (chromosome name)
        // 2. for each gene, we get the sequence of all its exons
        // 3. for each transcript, we join the transcripts' exon sequences to get the sequence of the transcript
        for result in reader.records() {
            let record = result?;
            let record_name = std::str::from_utf8(record.name())?;

            let chr_name = record_name.strip_suffix(' ').unwrap_or(record_name);

            // check if the chr_name is in all_seqnames
            // if not, we skip it
            if !all_seqnames.contains(chr_name) {
                continue;
            }

            let chr_gr = exon_gr.filter(seqname, &[chr_name], false)?;

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
                .str()?
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
                        let rec = &noodles::fasta::Record::new(definition, sequence);

                        let write_record: bool;
                        // call the callback if we have one
                        if let Some(ref mut cb) = record_filter {
                            write_record = cb(rec);
                        } else {
                            write_record = true;
                        }

                        if write_record {
                            writer
                                .write_record(rec)
                                .with_context(|| {
                                    format!(
                                        "Could not write the sequence of transcript {} to the output file",
                                        curr_tx
                                    )
                                })?;
                        }
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
            let rec = &noodles::fasta::Record::new(definition, sequence);

            let write_record: bool;
            // call the callback if we have one
            if let Some(ref mut cb) = record_filter {
                write_record = cb(rec);
            } else {
                write_record = true;
            }

            if write_record {
                writer.write_record(rec).with_context(|| {
                    format!(
                        "Could not write the sequence of transcript {} to the output file",
                        curr_tx
                    )
                })?;
            }
            exon_u8_vec.clear();
        }

        Ok(())
    }

    /// Writes sequences extracted from a reference genome based on the genomic features in the Grangers DataFrame to a FASTA file.
    ///
    /// This method reads genomic features from the current instance, extracts corresponding sequences from the reference genome,
    /// and writes them to a specified output file. Sequences can be named according to a specified column in the DataFrame or by default
    /// based on their row order. The method can also ignore strand information for sequence extraction.
    ///
    /// ### Generics
    ///
    /// * `T`: A type that implements [`AsRef<Path>`], providing a flexible reference for file paths.
    ///
    /// ### Arguments
    ///
    /// * `ref_path`: The file path to the reference genome sequences, typically in FASTA format.
    /// * `out_path`: The file path for the output FASTA file where extracted sequences will be written.
    /// * `ignore_strand`: A boolean indicating whether to ignore strand information when extracting sequences.
    /// * `name_column`: An optional string specifying the column name to use for naming extracted sequences.
    ///   If the column is invalid or not provided, sequences will be named based on their row order.
    /// * `oob_option`: A reference to [OOBOption] indicating how to handle features that extend beyond the bounds of the reference sequence.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<()>`](anyhow::Result) indicating the outcome:
    /// * [Ok]`(())`: Sequences were successfully extracted and written to the output file.
    /// * [Err]`(...)`: An error occurred during the process, such as validation failure, issues reading from the reference, or writing to the output file.
    ///
    /// ### Errors
    ///
    /// This function may return an error if:
    /// * There is a problem accessing the reference sequence file or the output file cannot be created.
    /// * The Grangers instance fails validation checks.
    /// * There are issues extracting sequences due to out-of-bound features or missing data.
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
        let out_file = BufWriter::with_capacity(4194304, out_file);

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
            "row_order".into(),
            (0..df.height() as u32).collect::<Vec<u32>>(),
        ))?;

        // if ignore strand, set the strand to +
        if ignore_strand {
            df.with_column(Series::new(fc.strand().into(), vec!["+"; df.height()]))?;
        }

        fc.fix(&df, false)?;

        let mut essential_gr = Grangers::new(df, None, None, IntervalType::default(), fc, false)?;
        essential_gr.set_signature(self.get_signature());

        let seqname = essential_gr.get_column_name_str("seqname", true)?;

        let mut reader = grangers_utils::get_noodles_reader_from_path(ref_path)?;

        let mut writer = noodles::fasta::io::writer::Builder::default()
            .set_line_base_count(usize::MAX)
            .build_from_writer(out_file);
        let mut empty_counter = 0;

        // we iterate the fasta reader. For each fasta reacord (usually chromosome), we do
        // 1. subset the dataframe by the chromosome name
        // 2. get the sequence of the features in the dataframe on that fasta record
        // 3. insert the sequence into the sequence vector according to the row order
        for result in reader.records() {
            let record = result?;
            let record_name = std::str::from_utf8(record.name())?;

            let chr_name = record_name.strip_suffix(' ').unwrap_or(record_name);
            let chr_gr = essential_gr.filter(seqname, &[chr_name], false)?;

            if chr_gr.df().height() == 0 {
                continue;
            }
            let name_vec = chr_gr
                .df()
                .column(name_column.as_str())?
                .str()?
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
                        .write_record(&noodles::fasta::Record::new(definition, sequence))
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

    /// Writes the sequences from a reference genome to a FASTA format output file based on the Grangers instance.
    ///
    /// This method extracts sequences from a reference genome based on the coordinates in the Grangers instance
    /// and writes them to the specified output file. The sequences can be named based on a specified column,
    /// or by default, they are named according to their row order in the DataFrame.
    ///
    /// ### Generics
    ///
    /// * `T`: A type that implements [`AsRef<Path>`], allowing for flexible path references to the reference genome file.
    /// * `W`: A type that implements [Write], designating the output stream for writing the sequences.
    ///
    /// ### Arguments
    ///
    /// * `ref_path`: The path to the reference genome FASTA file.
    /// * `out_file`: The output stream to which the sequences will be written.
    /// * `ignore_strand`: If true, the strand information is ignored, and all sequences are extracted in the forward direction.
    /// * `name_column`: Optional parameter specifying the name of the column to use for naming the extracted sequences. Defaults to row order if not provided.
    /// * `oob_option`: Out-of-bound option indicating how to handle features that exceed the reference sequence boundaries.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<()>`](anyhow::Result). If successful, sequences are written to the output file; otherwise, an error is returned.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let mut grangers = Grangers::new(...);
    /// let ref_path = "reference.fasta";
    /// let out_file = File::create("sequences.fasta")?;
    /// grangers.write_sequences(&ref_path, out_file, false, None, OOBOption::Skip)?;
    /// ```
    pub fn write_sequences<T: AsRef<Path>, W: Write, F>(
        &mut self,
        ref_path: T,
        out_file: W,
        ignore_strand: bool,
        name_column: Option<&str>,
        oob_option: OOBOption,
    ) -> anyhow::Result<()> {
        self.write_sequences_with_filter(
            ref_path,
            out_file,
            ignore_strand,
            name_column,
            oob_option,
            &mut None::<fn(&noodles::fasta::Record) -> bool>,
        )
    }

    /// Writes the sequences from a reference genome to a FASTA format output file based on the [Grangers] instance,
    /// applying an optional filter to each sequence before writing.
    ///
    /// This advanced method allows for filtering of sequences based on custom criteria before writing.
    /// The filter is a closure or function pointer provided by the user that accepts a `&noodles::fasta::Record`
    /// and returns a boolean indicating whether to write the sequence to the output.
    ///
    /// ### Generics
    ///
    /// * `T`: A type that implements [`AsRef<Path>`], allowing for flexible path references to the reference genome file.
    /// * `W`: A type that implements [Write], designating the output stream for writing the sequences.
    /// * `F`: A type that implements [FnMut]`(&`[noodles::fasta::Record]`) -> bool`, representing the filtering function.
    ///
    /// ### Arguments
    ///
    /// * `ref_path`: The path to the reference genome FASTA file.
    /// * `out_file`: The output stream to which the filtered sequences will be written.
    /// * `ignore_strand`: If true, the strand information is ignored, and all sequences are extracted in the forward direction.
    /// * `name_column`: Optional parameter specifying the name of the column to use for naming the extracted sequences. Defaults to row order if not provided.
    /// * `oob_option`: Out-of-bound option indicating how to handle features that exceed the reference sequence boundaries.
    /// * `record_filter`: Optional mutable reference to a filter function or closure to apply to each sequence.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<()>`](anyhow::Result). If successful, filtered sequences are written to the output file; otherwise, an error is returned.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let mut grangers = Grangers::new(...);
    /// let ref_path = "reference.fasta";
    /// let out_file = File::create("filtered_sequences.fasta")?;
    /// let mut filter = |record: &noodles::fasta::Record| record.sequence().len() > 100;
    /// grangers.write_sequences_with_filter(&ref_path, out_file, false, None, OOBOption::Skip, &mut Some(filter))?;
    /// ```
    pub fn write_sequences_with_filter<T: AsRef<Path>, W: Write, F>(
        &mut self,
        ref_path: T,
        out_file: W,
        ignore_strand: bool,
        name_column: Option<&str>,
        oob_option: OOBOption,
        record_filter: &mut Option<F>,
    ) -> anyhow::Result<()>
    where
        F: FnMut(&noodles::fasta::Record) -> bool,
    {
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
                .with_row_index("row_order".into(), None)?
                .select(selection)?
        } else {
            self.df.select(selection)?
        };

        // if ignore strand, set the strand to +
        if ignore_strand {
            df.with_column(Series::new(fc.strand().into(), vec!["+"; df.height()]))?;
        }

        fc.fix(&df, false)?;

        let mut essential_gr = Grangers::new(df, None, None, IntervalType::default(), fc, false)?;
        essential_gr.set_signature(self.get_signature());

        let seqname_s = essential_gr.get_column_name("seqname", true)?;
        let seqname = seqname_s.as_str();

        let mut reader = grangers_utils::get_noodles_reader_from_path(ref_path)?;

        let out_writer = BufWriter::with_capacity(4194304, out_file);
        let mut writer = noodles::fasta::io::writer::Builder::default()
            .set_line_base_count(usize::MAX)
            .build_from_writer(out_writer);

        let mut empty_counter = 0;

        // we iterate the fasta reader. For each fasta reacord (usually chromosome), we do
        // 1. subset the dataframe by the chromosome name
        // 2. get the sequence of the features in the dataframe on that fasta record
        // 3. insert the sequence into the sequence vector according to the row order
        for result in reader.records() {
            let record = result?;
            let record_name = std::str::from_utf8(record.name())?;

            let chr_name = record_name.strip_suffix(' ').unwrap_or(record_name);
            let chr_gr = essential_gr.filter(seqname, &[chr_name], false)?;

            if chr_gr.df().height() == 0 {
                continue;
            }

            let name_vec = chr_gr
                .df()
                .column(name_column.as_str())?
                .str()?
                .into_iter()
                .map(|s| s.unwrap())
                .collect::<Vec<_>>();

            let chrsi = ChrRowSeqIter::new(&chr_gr, &record, oob_option)?;

            for (feat_name, chrsi_rec) in name_vec.into_iter().zip(chrsi) {
                if let Ok(sequence) = chrsi_rec {
                    let definition = Definition::new(feat_name, None);
                    let rec = &noodles::fasta::Record::new(definition, sequence);

                    let write_record: bool;
                    // call the callback if we have one
                    if let Some(ref mut cb) = record_filter {
                        write_record = cb(rec);
                    } else {
                        write_record = true;
                    }

                    // we write if the sequence is not empty and
                    // it passes the filter (or there is no filter)
                    if write_record {
                        writer.write_record(rec).with_context(|| {
                            format!(
                                "Could not write sequence {} to the output file; Cannot proceed.",
                                feat_name
                            )
                        })?;
                    }
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

// implement get sequence functions for [Grangers]
impl Grangers {
    /// Extracts transcript sequences based on exon information from a FASTA file and returns them as a vector of FASTA records.
    ///
    /// This method processes the exon information within the [Grangers] instance to extract sequences corresponding
    /// to each transcript from a given reference genome. It then compiles these sequences into FASTA records.
    ///
    /// ### Generics
    ///
    /// * `T`: A type that implements [`AsRef<Path>`], allowing for flexible path references to the FASTA file.
    ///
    /// ### Arguments
    ///
    /// * `fasta_path`: The path to the reference genome FASTA file.
    /// * `exon_name`: Optional parameter specifying the name of the exon feature. Defaults to `"exon"` if not provided.
    /// * `multithreaded`: Boolean flag indicating whether to use multithreading for faster processing.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<Vec<noodles::fasta::Record>>`]. If successful, a vector of transcript sequences
    /// encoded as FASTA records is returned; otherwise, an error is returned.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let mut grangers = Grangers::new(...);
    /// let fasta_path = "reference.fasta";
    /// let transcript_sequences = grangers.get_transcript_sequences(&fasta_path, None, false)?;
    /// for seq in transcript_sequences {
    ///     println!("{}", seq);
    /// }
    /// ```
    ///
    /// ### Errors
    ///
    /// This function can return an error if there are issues with reading the FASTA file,
    /// if exon records are invalid or exceed the reference sequence, or if other validation steps fail.
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
        let mut reader = grangers_utils::get_noodles_reader_from_path(fasta_path)?;

        // let mut seq_vec: Vec<Option<Sequence>> = vec![None; exon_gr.df().height()];
        let mut transcript_seq_vec: Vec<noodles::fasta::Record> =
            Vec::with_capacity(self.df().column(transcript_id)?.unique()?.len());

        // we iterate the fasta reader. For each fasta reacord (usually chromosome), we do
        // 1. subset the dataframe by the seqname (chromosome name)
        // 2. for each gene, we get the sequence of all its exons
        // 3. for each transcript, we join the transcripts' exon sequences to get the sequence of the transcript
        for result in reader.records() {
            let record = result?;
            let record_name = std::str::from_utf8(record.name())?;

            let chr_name = record_name.strip_suffix(' ').unwrap_or(record_name);
            let chr_gr = exon_gr.filter(seqname, &[chr_name], false)?;

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
                .str()?
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
                        transcript_seq_vec.push(noodles::fasta::Record::new(definition, sequence));
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

    /// Extract the sequence of the features in the [Grangers] object from the provided reference file.
    /// Currently only fasta file is supported. This function four field columns: seqname, start, end, and strand.
    /// Arguments:
    /// - `genome_path`: the path to the reference genome file.
    /// - `file_format`: the format of the reference genome file. Currently only fasta is supported.
    /// - `oob_option`: the option for out-of-boundary positions. It can be either `Truncate` or `Skip`. If `Truncate`, the out-of-boundary positions will be truncated to the start or end of the sequence. If `Skip`, a None will be returned for features with OOB positions
    ///
    /// The function outputs the extracted sequence as a vector of `Option<Sequence>`. If the feature has OOB positions and the oob_option is set as `Skip`, the corresponding element in the vector will be None. The order of the vector follows the row order of the dataframe of the [Grangers] object.
    pub fn _get_sequences<T: AsRef<Path>>(
        &mut self,
        fasta_path: T,
        ignore_strand: bool,
        name: Option<&str>,
        oob_option: &OOBOption,
    ) -> anyhow::Result<Vec<Option<noodles::fasta::Record>>> {
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
            .with_row_index("row_order", None)
            .with_column(if ignore_strand {
                lit("+").alias("strand")
            } else {
                col("strand")
            })
            .collect()?;

        fc.fix(&df, false)?;

        let mut essential_gr = Grangers::new(df, None, None, IntervalType::default(), fc, false)?;
        essential_gr.set_signature(self.get_signature());

        let seqname = essential_gr.get_column_name_str("seqname", true)?;

        let mut reader = grangers_utils::get_noodles_reader_from_path(fasta_path)?;

        let mut seq_vec: Vec<Option<noodles::fasta::Record>> =
            vec![None; essential_gr.df().height()];
        // we iterate the fasta reader. For each fasta reacord (usually chromosome), we do
        // 1. subset the dataframe by the chromosome name
        // 2. get the sequence of the features in the dataframe on that fasta record
        // 3. insert the sequence into the sequence vector according to the row order
        for result in reader.records() {
            let record = result?;
            let record_name = std::str::from_utf8(record.name())?;

            let chr_name = record_name.strip_suffix(' ').unwrap_or(record_name);
            let chr_gr = essential_gr.filter(seqname, &[chr_name], false)?;

            if chr_gr.df().height() == 0 {
                continue;
            }
            let name_vec = if let Some(name) = &name {
                chr_gr
                    .df()
                    .column(name)?
                    .str()?
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

                seq_vec[idx] = seq.map(|seq| noodles::fasta::Record::new(definition, seq));
            }
        }

        Ok(seq_vec)
    }

    /// Extracts sequences from a reference FASTA file and returns them as a [GrangersSequenceCollection].
    ///
    /// This method reads the reference genome from a FASTA file, filters features based on their presence in the genome,
    /// and extracts their sequences. These sequences are then compiled into a [GrangersSequenceCollection],
    /// which includes a unique signature and a vector of records for further processing.
    ///
    /// ### Generics
    ///
    /// * `T`: A type that implements [`AsRef<Path>`], allowing for flexible path references to the FASTA file.
    ///
    /// ### Arguments
    ///
    /// * `ref_path`: The path to the reference genome FASTA file.
    /// * `ignore_strand`: Boolean flag indicating whether to ignore the strand information during sequence extraction.
    /// * `name_column`: Optional parameter specifying the column name to use for naming the extracted sequences.
    ///                   If not provided, sequences will be named based on their row order.
    /// * `oob_option`: Specifies how to handle features that go out of the reference sequence bounds.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<GrangersSequenceCollection>`]. If successful, a collection of extracted sequences
    /// along with a unique signature is returned; otherwise, an error is returned.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let mut grangers = Grangers::new(...);
    /// let fasta_path = "reference.fasta";
    /// let sequence_collection = grangers.get_sequences(&fasta_path, false, None, OOBOption::Trim)?;
    /// println!("{:?}", sequence_collection);
    /// ```
    ///
    /// ### Errors
    ///
    /// This function can return an error if there are issues with reading the FASTA file,
    /// if there are validation errors, or if other processing steps fail.
    pub fn get_sequences<T: AsRef<Path>>(
        &mut self,
        ref_path: T,
        ignore_strand: bool,
        name_column: Option<&str>,
        oob_option: OOBOption,
    ) -> anyhow::Result<GrangersSequenceCollection> {
        self.get_sequences_from_read(
            std::fs::File::open(ref_path)?,
            ignore_strand,
            name_column,
            oob_option,
        )
    }

    /// Extracts sequences from a FASTA format reader and returns them as a [GrangersSequenceCollection].
    ///
    /// This method reads the reference genome from a given reader implementing the [Read] trait,
    /// filters features based on their presence in the genome, and extracts their sequences. These sequences
    /// are then compiled into a [GrangersSequenceCollection], which includes a unique signature and a vector
    /// of records for further processing.
    ///
    /// ### Generics
    ///
    /// * `R`: A type that implements [Read], allowing for flexible reading of the FASTA data.
    ///
    /// ### Arguments
    ///
    /// * `reader`: A reader instance from which the reference genome will be read.
    /// * `ignore_strand`: Boolean flag indicating whether to ignore the strand information during sequence extraction.
    /// * `name_column`: Optional parameter specifying the column name to use for naming the extracted sequences.
    ///                   If not provided, sequences will be named based on their row order.
    /// * `oob_option`: Specifies how to handle features that go out of the reference sequence bounds.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<GrangersSequenceCollection>`]. If successful, a collection of extracted sequences
    /// along with a unique signature is returned; otherwise, an error is returned.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let mut grangers = Grangers::new(...);
    /// let fasta_data = "reference data here...";
    /// let sequence_collection = grangers.get_sequences_from_read(fasta_data.as_bytes(), false, None, OOBOption::Trim)?;
    /// println!("{:?}", sequence_collection);
    /// ```
    ///
    /// ### Errors
    ///
    /// This function can return an error if there are issues with reading the FASTA data,
    /// if there are validation errors, or if other processing steps fail.
    pub fn get_sequences_from_read<R: Read + 'static>(
        &mut self,
        reader: R,
        ignore_strand: bool,
        name_column: Option<&str>,
        oob_option: OOBOption,
    ) -> anyhow::Result<GrangersSequenceCollection> {
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
            .with_row_index("row_order", None)
            .with_column(if ignore_strand {
                lit("+").alias("strand")
            } else {
                col("strand")
            })
            .collect()?;

        // we only use a subset of the columns, so fix fc
        fc.fix(&df, false)?;

        let mut essential_gr = Grangers::new(df, None, None, IntervalType::default(), fc, false)?;
        essential_gr.set_signature(self.get_signature());

        let seqname = essential_gr.get_column_name_str("seqname", true)?;
        let mut reader = grangers_utils::get_noodles_reader_from_reader(reader)?;

        let sig = essential_gr.get_signature();
        let num_rec = essential_gr.df().height();
        let mut seq_coll =
            GrangersSequenceCollection::new_with_signature_and_capacity(sig, num_rec);
        let mut empty_counter = 0_usize;

        // we iterate the fasta reader. For each fasta reacord (usually chromosome), we do
        // 1. subset the dataframe by the chromosome name
        // 2. get the sequence of the features in the dataframe on that fasta record
        // 3. insert the sequence into the sequence vector according to the row order
        for result in reader.records() {
            let record = result?;
            let record_name = std::str::from_utf8(record.name())?;

            let chr_name = record_name.strip_suffix(' ').unwrap_or(record_name);
            let chr_gr = essential_gr.filter(seqname, &[chr_name], false)?;

            if chr_gr.df().height() == 0 {
                continue;
            }

            let name_vec_iter = chr_gr
                .df()
                .column(name_column.as_str())?
                .cast(&DataType::String)?
                .str()?
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

                        // don't add anything to the sequence
                        // collection in this case
                        empty_counter += 1;
                        continue;
                    }
                };

                let definition = Definition::new(feat_name, None);
                let record = noodles::fasta::Record::new(definition, sequence);

                //seq_vec[idx as usize] = Some(record);
                // add this to the sequence collection
                seq_coll.add_record(GrangersRecordID::new(idx), record);
            }
        }

        if empty_counter > 0 {
            warn!("Unable to extract sequence for {} records. They are usually caused by out of boundary features or an invalid alphabet.", empty_counter)
        }
        Ok(seq_coll)
    }

    /// Returns an iterator over sequences extracted from a reference genome provided via a reader implementing the [Read] trait.
    ///
    /// This method creates an iterator that lazily reads and processes sequences from the reference genome.
    /// It filters and extracts feature sequences based on the current [Grangers] instance's configuration
    /// and returns them in a new iterator. This method is particularly useful for processing large genomes
    /// in a memory-efficient manner.
    ///
    /// ### Generics
    ///
    /// * `R`: A type that implements [Read], allowing for reading of the FASTA data from various sources.
    ///
    /// ### Arguments
    ///
    /// * `reader`: A reader instance from which the reference genome will be read.
    /// * `ignore_strand`: Boolean flag indicating whether to ignore strand information during sequence extraction.
    /// * `name_column`: Optional parameter specifying the column name to use for naming the extracted sequences.
    ///                   If not provided, sequences will be named based on their row order.
    /// * `oob_option`: Specifies how to handle features that go out of the reference sequence bounds.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<Pin<Box<GrangersSeqIter<R>>>>`]. If successful, an iterator over extracted sequences
    /// is returned; otherwise, an error is returned.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let mut grangers = Grangers::new(...);
    /// let reader = BufReader::new(File::open("reference.fasta")?);
    /// let sequence_iterator = grangers.iter_sequences_from_reader(reader, false, None, OOBOption::Trim)?;
    /// for sequence in sequence_iterator {
    ///     println!("{:?}", sequence);
    /// }
    /// ```
    ///
    /// ### Errors
    ///
    /// This function can return an error if there are validation errors, or if there are issues
    /// setting up the iterator or reading data from the provided reader.
    pub fn iter_sequences_from_reader<R: Read + 'static>(
        &mut self,
        reader: R,
        ignore_strand: bool,
        name_column: Option<&str>,
        oob_option: OOBOption,
    ) -> anyhow::Result<Pin<Box<GrangersSeqIter>>> {
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
            .with_row_index("row_order", None)
            .with_column(if ignore_strand {
                lit("+").alias("strand")
            } else {
                col("strand")
            })
            .collect()?;

        // we only use a subset of the columns, so fix fc
        fc.fix(&df, false)?;

        let mut essential_gr = Grangers::new(df, None, None, IntervalType::default(), fc, false)?;
        essential_gr.set_signature(self.get_signature());

        let seqname = essential_gr.get_column_name_str("seqname", true)?;

        let filt_opt = GrangersFilterOpts {
            seqname: seqname.to_owned(),
            name_column,
            oob_option,
        };

        Ok(GrangersSeqIter::new(reader, filt_opt, essential_gr))
    }

    /// Returns an iterator over sequences extracted from a reference genome provided via a file path.
    ///
    /// This method is a convenience wrapper around [iter_sequences_from_reader](fn@Grangers::iter_sequences_from_reader) that opens the reference genome file
    /// and creates an iterator to lazily read and process sequences. It is suitable for processing large genomes
    /// where loading all sequences into memory is not feasible.
    ///
    /// ### Generics
    ///
    /// * `T`: A type that implements [`AsRef<Path>`], allowing for flexible file path references.
    ///
    /// ### Arguments
    ///
    /// * `ref_path`: The path to the reference genome FASTA file.
    /// * `ignore_strand`: Boolean flag indicating whether to ignore strand information during sequence extraction.
    /// * `name_column`: Optional parameter specifying the column name to use for naming the extracted sequences.
    ///                   If not provided, sequences will be named based on their row order.
    /// * `oob_option`: Specifies how to handle features that go out of the reference sequence bounds.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<Pin<Box<GrangersSeqIter<std::fs::File>>>>`]. If successful, an iterator over extracted sequences
    /// is returned; otherwise, an error is returned.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let mut grangers = Grangers::new(...);
    /// let sequence_iterator = grangers.iter_sequences("reference.fasta", false, None, OOBOption::Trim)?;
    /// for sequence in sequence_iterator {
    ///     println!("{:?}", sequence);
    /// }
    /// ```
    ///
    /// ### Errors
    ///
    /// This function can return an error if the file cannot be opened, if there are validation errors, or if there are
    /// issues setting up the iterator.
    pub fn iter_sequences<T: AsRef<Path>>(
        &mut self,
        ref_path: T,
        ignore_strand: bool,
        name_column: Option<&str>,
        oob_option: OOBOption,
    ) -> anyhow::Result<Pin<Box<GrangersSeqIter>>> {
        self.iter_sequences_from_reader(
            std::fs::File::open(ref_path)?,
            ignore_strand,
            name_column,
            oob_option,
        )
    }

    /// Extracts sequences from a given FASTA record based on the feature information contained within the [Grangers] instance.
    ///
    /// This method processes a single [noodles::fasta::Record], typically representing a chromosome or scaffold,
    /// and extracts sequences according to the feature data (e.g., gene or exon locations) contained within the [Grangers] instance.
    /// It supports handling sequences that exceed reference boundaries based on the specified [OOBOption].
    ///
    /// ### Arguments
    ///
    /// * `record`: A reference to a `noodles::fasta::Record` from which sequences will be extracted.
    /// * `oob_option`: A reference to an `OOBOption` enum determining how out-of-bound sequences should be handled.
    ///                 Options include truncating sequences at the reference boundaries or skipping them entirely.
    ///
    /// ### Returns
    ///
    /// Returns an [`anyhow::Result<Vec<Option<Sequence>>>`]. Each element in the returned vector corresponds to a sequence extracted
    /// from the FASTA record. The sequences are aligned with the features in the [Grangers] instance. `None` is used to represent
    /// sequences that could not be extracted (e.g., due to out-of-bound issues).
    ///
    /// ### Example
    ///
    /// ```rust
    /// let fasta_record = noodles::fasta::Record::new(...);
    /// let sequences = grangers.get_sequences_fasta_record(&fasta_record, &OOBOption::Truncate)?;
    /// for seq_option in sequences {
    ///     match seq_option {
    ///         Some(seq) => println!("{:?}", seq),
    ///         None => println!("Sequence out of bounds"),
    ///     }
    /// }
    /// ```
    ///
    /// ### Errors
    ///
    /// This function can return an error if there are issues with data validation, conversion of sequence positions,
    /// or if the FASTA record does not match the reference name specified in the [Grangers] instance.
    /// It also fails if the dataframe contains more than one reference name, indicating that filtering by the reference
    /// name is required prior to calling this method.
    pub(crate) fn get_sequences_fasta_record(
        &self,
        record: &noodles::fasta::Record,
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
            .zip(ses[2].str()?.into_iter())
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

/// Configuration options for filtering and extracting sequences within the `Grangers` framework.
///
/// This structure contains parameters used to specify how genomic sequences should be filtered
/// and named when being extracted from a reference sequence. It is typically used in conjunction
/// with sequence extraction methods to provide additional context and control over the extraction process.
///
/// ### Fields
///
/// * `seqname`: A [String] representing the name of the sequence (e.g., chromosome name) to which the filtering and extraction will be applied.
/// * `name_column`: A [String] specifying the column in the [Grangers] instance's dataframe that contains names for the extracted sequences.
///               If this column is invalid or not present, a fallback mechanism such as row order might be used.
/// * `oob_option`: An [OOBOption] enum value that determines how out-of-bound (OOB) sequences should be handled during extraction.
///                This could include options such as truncating the sequences at the reference boundaries or skipping them entirely.
///
pub struct GrangersFilterOpts {
    seqname: String,
    name_column: String,
    oob_option: OOBOption,
}

/// Iterator for genomic sequences based on the Grangers data structure.
///
/// This struct encapsulates the functionality needed to iterate over sequences extracted from a genomic dataset,
/// with the ability to apply specific filtering and transformation criteria defined by `GrangersFilterOpts`.
///
/// # Fields
///
/// * `essential_gr`: A [Grangers] instance that holds essential data fields required for processing all target sequences.
///    This serves as the base data structure from which specific sequences are extracted.
///
/// * `chr_gr`: An optional [Grangers] instance containing only the features relevant to the current target sequence.
///    This is dynamically updated to match the current focus of sequence extraction.
///
/// * `seq_reader`: A FASTA format reader from the `noodles` crate, wrapping an underlying reader
///    that depends on wether or not the source is compressed. It is used for reading sequence data from a reference
///    genome or other source.
///
/// * `seq_record`: Represents the current sequence record being processed by the iterator.
///    It holds both the sequence identifier and the actual sequence data.
///
/// * `filt_opt`: Filter options encapsulated within a [GrangersFilterOpts] structure. These options dictate how sequences
///    should be filtered and processed during iteration, including which sequences to include and how to handle edge cases.
///
/// * `name_vec_iter`: An iterator over the names of sequences that need to be extracted based on the current dataset.
///    This typically corresponds to identifiers like transcript or gene IDs.
///
/// * `row_order_iter`: An iterator over the row indices of sequences in the dataset, providing a link between sequence data
///    and their corresponding metadata or annotations within the [Grangers] structure.
///
/// * `chr_seq_iter`: An optional internal iterator (`ChrRowSeqIter`) that handles the iteration over individual sequence features
///    for a given target, such as exons within a transcript. This allows for fine-grained processing of sequences.
///
/// * `def_buffer`: A local buffer used to hold sequence definitions temporarily. This can be used for building FASTA headers
///    or other metadata strings associated with each sequence.
///
/// # Usage
///
/// This iterator is designed to be used in scenarios where sequences need to be extracted and processed from a larger genomic dataset.
/// It is particularly useful for applications that require iterating over sequences with specific filtering criteria, such as extracting
/// all exons from a set of transcripts.
pub struct GrangersSeqIter {
    // the essential grangers struct holding
    // the required fields across *all* of the
    // target sequences
    essential_gr: Grangers,
    // the filtered Grangers struct that
    // contains only the features for the
    // current target
    chr_gr: Option<Grangers>,
    // a noodles Fasta reader for reading the
    // target sequences
    seq_reader: grangers_utils::FastaReader,
    // the current seq record
    seq_record: noodles::fasta::Record,
    // the filter options that will be applied
    filt_opt: GrangersFilterOpts,
    // the iterator over the names of the sequences
    // we need to extract.
    name_vec_iter: <Vec<String> as IntoIterator>::IntoIter,
    // the iterator over the row order (row indices) of
    // the sequences we need to extract.
    row_order_iter: <Vec<u32> as IntoIterator>::IntoIter,
    // the "inner" iterator that iterates over the sequence
    // features of an individual target.
    chr_seq_iter: Option<ChrRowSeqIter<'static>>,
    // local buffer to hold the sequence definition
    // string.
    def_buffer: String,
}

use core::pin::Pin;

impl GrangersSeqIter {
    /// Creates a new instance of the [GrangersSeqIter].
    ///
    /// This constructor initializes a new sequence iterator for processing genomic data.
    /// It sets up the necessary internal state, including a FASTA reader for reading sequences,
    /// and prepares the iterator with user-defined filter options and essential Granger data.
    ///
    /// # Arguments
    ///
    /// * `breader`: A reader implementing the [`Read`] trait. This reader is used
    ///   to stream genomic sequence data from FASTA-formatted files or other readable sources.
    ///
    /// * `filt_opt`: Filter options encapsulated within a [GrangersFilterOpts] structure. These options
    ///   dictate how sequences should be filtered and processed during iteration, including which sequences
    ///   to include and how to handle sequences that extend beyond reference boundaries (OOB).
    ///
    /// * `essential_gr`: A [Grangers] instance that holds the essential fields and dataset needed for sequence extraction.
    ///   This provides the context in which sequence data will be interpreted and processed.
    ///
    /// # Returns
    ///
    /// Returns a [`Pin<Box<GrangersSeqIter<R>>>`]: a pinned, heap-allocated instance of the iterator.
    /// This pinning is necessary to ensure the stability of the internal references due to the self-referential nature
    /// of streaming iterators.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// let file = File::open("path/to/reference.fasta")?;
    /// let breader = BufReader::new(file);
    /// let filt_opts = GrangersFilterOpts { ... };
    /// let essential_gr = Grangers::new(...);
    ///
    /// let seq_iter = GrangersSeqIter::new(breader, filt_opts, essential_gr);
    /// ```
    ///
    /// # Note
    ///
    /// The initial `seq_record` is set to a default "empty" record; actual sequence data will be populated
    /// when the iterator is advanced.
    pub fn new<R: Read + 'static>(
        r: R,
        filt_opt: GrangersFilterOpts,
        essential_gr: Grangers,
    ) -> Pin<Box<Self>> {
        let reader = grangers_utils::get_noodles_reader_from_reader(r)
            .expect("couldn't create reader from input reader");
        let definition = Definition::new("empty", None);
        let sequence = Sequence::from(b"A".to_vec());
        let curr_record = noodles::fasta::Record::new(definition, sequence);
        let v: Vec<String> = vec![];
        let o: Vec<u32> = vec![];
        Box::pin(GrangersSeqIter {
            essential_gr,
            chr_gr: None,
            seq_reader: reader,
            seq_record: curr_record,
            filt_opt,
            name_vec_iter: v.into_iter(),
            row_order_iter: o.into_iter(),
            chr_seq_iter: None,
            def_buffer: String::new(),
        })
    }
}

impl Iterator for GrangersSeqIter {
    type Item = (GrangersRecordID, noodles::fasta::Record);

    #[allow(clippy::question_mark)]
    /// Advances the iterator and returns the next genomic sequence.
    ///
    /// This method iterates through genomic sequences extracted based on the Grangers framework.
    /// It processes sequences chromosome by chromosome and feature by feature according to the configured filters.
    ///
    /// # Returns
    ///
    /// Returns [`Some((GrangersRecordID, noodles::fasta::Record))`] when a new sequence is successfully extracted, where:
    /// - [GrangersRecordID] is an identifier corresponding to the row order of the sequence in the DataFrame.
    /// - [noodles::fasta::Record] is the actual sequence extracted and formatted according to FASTA standards.
    ///
    /// Returns [None] when all sequences have been iterated over or when the iterator encounters an error from which it cannot recover.
    ///
    /// # Panics
    ///
    /// This method may panic if:
    /// - There are issues reading definitions or sequences from the FASTA file.
    /// - There are issues parsing the sequence definition.
    /// - The expected genomic feature columns are missing or contain invalid data.
    ///
    /// # Examples
    ///
    /// Assuming `grangers_seq_iter` is an instance of `GrangersSeqIter<R>`:
    /// ```ignore
    /// while let Some((id, record)) = grangers_seq_iter.next() {
    ///     println!("ID: {:?}, Sequence: {:?}", id, record);
    /// }
    /// ```
    ///
    /// # Safety
    ///
    /// This implementation uses unsafe code to extend the lifetime of references to [Grangers] and [noodles::fasta::Record].
    /// It is crucial that these references remain valid for the duration of the iterator's usage.
    /// Misuse may lead to undefined behavior, such as use-after-free errors.
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // check if the inner iterator exists and, if so, if we have
            // elements to yield from it.
            if let Some(ref mut chr_seq_iter) = self.chr_seq_iter {
                if let (Some(chr_seq_rec), Some(feat_name), Some(row_idx)) = (
                    chr_seq_iter.next(),
                    self.name_vec_iter.next(),
                    self.row_order_iter.next(),
                ) {
                    if let Ok(sequence) = chr_seq_rec {
                        return Some((
                            GrangersRecordID::new(row_idx),
                            noodles::fasta::Record::new(Definition::new(feat_name, None), sequence),
                        ));
                    }
                    // if we don't have a sequence (this was empty), we want to go to the next
                    // iteration of this loop immediately.
                    continue;
                }
                // NOTE: we should assert somewhere that `chr_seq_iter` and `name_vec_iter` have
                // the same size. Mostly make sense only if they have exact size hints.
                // If we exhausted the iterator, then we want to set the chr_seq_iter to none and
                // go to the top of the loop.
                self.chr_seq_iter = None;
            } else {
                // in this branch, chr_seq_iter was None, so either we exhausted the previous
                // inner iterator, or we haven't created it yet.
                // we iterate the fasta reader. For each fasta reacord (usually chromosome), we do
                // 1. subset the dataframe by the chromosome name
                // 2. get the sequence of the features in the dataframe on that fasta record
                // 3. yield the iterator over that data frame
                loop {
                    self.def_buffer.clear();
                    let def_bytes = self
                        .seq_reader
                        .read_definition(&mut self.def_buffer)
                        .expect("GrangersSeqIter: could not read definition from reference file");

                    // if we reached the end of the file, exhaust the outer iterator
                    if def_bytes == 0 {
                        return None;
                    }

                    let definition = match self.def_buffer.parse() {
                        Ok(d) => d,
                        Err(e) => panic!("could not parse sequence definition: error {}", e),
                    };

                    let mut seq_buffer = Vec::<u8>::new();
                    let seq_bytes = self
                        .seq_reader
                        .read_sequence(&mut seq_buffer)
                        .expect("GrangersSeqIter: could not read sequence from reference file");
                    if seq_bytes == 0 {
                        warn!("GrangersSeqIter: was able to read record definition, but no sequence. This seems like a problem!");
                        return None;
                    }
                    let sequence = Sequence::from(seq_buffer);

                    // at this point we have the next sequence record
                    self.seq_record = noodles::fasta::Record::new(definition, sequence);
                    let record_name = std::str::from_utf8(self.seq_record.name())
                        .expect("GrangersSeqIter: could not convert record name to utf8");

                    let chr_name = record_name.strip_suffix(' ').unwrap_or(record_name);

                    self.chr_gr = Some(
                        self.essential_gr
                            .filter(
                                self.filt_opt.seqname.clone(),
                                &[chr_name.to_string()],
                                false,
                            )
                            .expect("GrangersSeqIter: cannot filter essential_gr"),
                    );

                    if self.chr_gr.as_ref().unwrap().df().height() == 0 {
                        self.chr_gr = None;
                        continue;
                    }

                    self.name_vec_iter = self
                        .chr_gr
                        .as_ref()
                        .unwrap()
                        .df()
                        .column(self.filt_opt.name_column.as_str())
                        .expect("GrangersSeqIter: cannot get name_column")
                        .str()
                        .expect("GrangersSeqIter: cannot convert name_vec to str")
                        .into_iter()
                        .map(|s| s.unwrap().to_owned())
                        .collect::<Vec<_>>()
                        .into_iter();

                    self.row_order_iter = self
                        .chr_gr
                        .as_ref()
                        .unwrap()
                        .df()
                        .column("row_order")
                        .expect("GrangersSeqIter: cannot get row_order column")
                        .u32()
                        .expect("GrangerSeqIter: cannot convert row_order entries to u32")
                        .into_iter()
                        .map(|s| {
                            s.expect("Could not get row order. Please report this bug on GitHub.")
                        })
                        .collect::<Vec<_>>()
                        .into_iter();

                    // convert the reference to a static reference
                    // this is generally highly unsafe and can lead to
                    // use after frees, but if chr_seq_iter is never moved out
                    // of this struct it shold be ok
                    // also we must use Pin on the struct so that the references
                    // are never invalidated
                    let ref_grangers = unsafe {
                        core::mem::transmute::<&Grangers, &'static Grangers>(
                            self.chr_gr.as_ref().unwrap(),
                        )
                    };
                    let ref_record = unsafe {
                        core::mem::transmute::<
                            &noodles::fasta::Record,
                            &'static noodles::fasta::Record,
                        >(&self.seq_record)
                    };
                    self.chr_seq_iter = Some(
                        ChrRowSeqIter::new(ref_grangers, ref_record, self.filt_opt.oob_option)
                            .expect("cannot create ChrRowSeqIter"),
                    );
                    break;
                } // loop over getting the next chr_seq_iter

                // if we got to this point, and we weren't able to fill in
                // self.chr_seq_iter, then the iterator should be exhausted
                if self.chr_seq_iter.is_none() {
                    return None;
                }
            }
        }
    }
}

/// Iterator for traversing over genomic feature sequences within a single chromosome or sequence record.
///
/// This struct holds the iterators for different genomic feature columns from a [DataFrame],
/// along with a reference to the current FASTA record representing the chromosome or sequence segment.
/// It's designed to be used for extracting sequences for features like exons or genes, with consideration
/// for how out-of-bounds sequences should be handled.
///
/// # Lifetime
///
/// * `'a`: The lifetime parameter `'a` ties [ChrRowSeqIter] to the lifetime of the FASTA record from which
///   it extracts sequences, ensuring that the record remains valid for the duration of the iterator's use.
///
/// # Fields
///
/// * `iters`: A vector of [`SeriesIter<'a>`], where each [SeriesIter] is an iterator over a column from a DataFrame
///   associated with the genomic features (e.g., start and end positions, strand). These iterators are used to traverse
///   the feature data and extract corresponding sequences.
///
/// * `record`: A reference to a [noodles::fasta::Record], which contains the sequence of the current chromosome
///   or sequence segment. This is the reference sequence from which genomic feature sequences are extracted.
///
/// * `oob_option`: An instance of [OOBOption] determining how to handle genomic features that extend beyond the bounds
///   of the `record` sequence. This could involve truncating or skipping such out-of-bounds features.
///
/// * `seqlen`: The length of the sequence within the current `record`. This is used to validate feature positions
///   and handle out-of-bound situations according to `oob_option`.
///
struct ChrRowSeqIter<'a> {
    starts: Vec<Option<i64>>,
    ends: Vec<Option<i64>>,
    is_reverse: Vec<Option<bool>>,
    curr_idx: usize,
    record: &'a noodles::fasta::Record,
    oob_option: OOBOption,
    seqlen: usize,
}

impl<'a> ChrRowSeqIter<'a> {
    /// Creates a new instance of [ChrRowSeqIter].
    ///
    /// This constructor sets up an iterator for processing genomic features within a single chromosome or sequence
    /// segment, based on the data contained within a [Grangers] instance and a [noodles::fasta::Record].
    ///
    /// # Arguments
    ///
    /// * `grangers`: A reference to a [Grangers] instance containing the genomic feature data for extraction,
    ///   including necessary columns like start, end, and strand of each feature.
    ///
    /// * `record`: A reference to a [noodles::fasta::Record] representing the sequence from which the genomic features
    ///   will be extracted. This typically corresponds to a single chromosome or scaffold.
    ///
    /// * `oob_option`: An [OOBOption] determining how to handle genomic features that extend beyond the sequence boundaries
    ///   defined by `record`. Options may include skipping such features, truncating them to fit within bounds, or other behaviors.
    ///
    /// # Returns
    ///
    /// Returns an [`anyhow::Result<Self>`](anyhow::Result):
    /// * [Ok]`(Self)`: If the iterator is successfully created.
    /// * [Err]`(...)`: If there is an error initializing the iterator, such as missing columns in the [Grangers] dataframe or issues accessing the sequence.
    ///
    /// # Examples
    ///
    /// Assuming `grangers` is an initialized [Grangers] instance and `fasta_record` is a [noodles::fasta::Record] for the relevant sequence:
    ///
    /// ```ignore
    /// let chr_row_seq_iter = ChrRowSeqIter::new(&grangers, &fasta_record, OOBOption::Skip)?;
    /// ```
    ///
    /// # Note
    ///
    /// This method collects iterators from the specified [Grangers] instance for the start, end, and strand columns. These iterators
    /// are then used to traverse the feature data and extract corresponding sequences based on the provided FASTA record.
    pub fn new(
        grangers: &'a Grangers,
        record: &'a noodles::fasta::Record,
        oob_option: OOBOption,
    ) -> anyhow::Result<Self> {
        let fc = grangers.field_columns();

        // Retrieve Series and rechunk/cast as needed to ensure safety and correctness
        // We collect to Vecs to avoid lifetime issues with local rechunked Series
        let s_start = grangers.df().column(fc.start())?.as_materialized_series().rechunk();
        let s_end = grangers.df().column(fc.end())?.as_materialized_series().rechunk();
        let s_strand = grangers.df().column(fc.strand())?.as_materialized_series().rechunk();

        // Collect starts and ends
        // We use i64() helper which returns ChunkedArray<Int64Type>
        let starts: Vec<Option<i64>> = s_start.i64()?.into_iter().collect();
        let ends: Vec<Option<i64>> = s_end.i64()?.into_iter().collect();

        // Handle strand
        // Cast to String to ensure we can read it as string, then collect boolean flag
        let s_strand_str = s_strand.cast(&DataType::String)?;
        let ca_strand = s_strand_str.str()?;
        let is_reverse: Vec<Option<bool>> = ca_strand.into_iter().map(|opt_s| {
            opt_s.map(|s| s == "-")
        }).collect();

        let seqlen = record.sequence().len();
        Ok(Self {
            starts,
            ends,
            is_reverse,
            curr_idx: 0,
            record,
            oob_option,
            seqlen,
        })
    }
}

#[allow(clippy::needless_lifetimes)]
impl<'a> Iterator for ChrRowSeqIter<'a> {
    type Item = anyhow::Result<Sequence>;

    /// Advances the iterator and returns the next genomic sequence.
    ///
    /// This method sequentially processes genomic feature data to extract corresponding sequences
    /// from the associated FASTA record. It respects the out-of-bound (OOB) handling strategy
    /// specified during initialization and accounts for feature orientation by handling reverse-complement
    /// sequences as needed.
    ///
    /// # Returns
    ///
    /// Returns [`Some(Ok(Sequence))`] when a new sequence is successfully extracted and complies
    /// with the provided feature data and OOB strategy.
    ///
    /// Returns [`Some(Err(...))`] when there is an issue extracting a sequence, such as invalid
    /// start/end positions, null field values, or sequence indexing errors.
    ///
    /// Returns [None] when all feature sequences have been iterated over, indicating the end
    /// of the genomic features in the DataFrame.
    ///
    /// # Examples
    ///
    /// Assuming `chr_row_seq_iter` is an instance of [`ChrRowSeqIter<'a>`]:
    /// ```ignore
    /// while let Some(result) = chr_row_seq_iter.next() {
    ///     match result {
    ///         Ok(sequence) => println!("Extracted sequence: {:?}", sequence),
    ///         Err(e) => println!("Error extracting sequence: {}", e),
    ///     }
    /// }
    /// ```
    ///
    /// # Note
    ///
    /// This method performs several checks to ensure the integrity of the extracted sequence,
    /// including validating start/end positions and handling strand-specific sequence extraction.
    /// It leverages the [`OOBOption`] settings to decide how sequences extending beyond the reference
    /// are treated, either truncating them to fit or skipping them entirely.
    fn next(&mut self) -> Option<Self::Item> {
        if self.curr_idx >= self.starts.len() {
            return None;
        }

        let start_opt = self.starts[self.curr_idx];
        let end_opt = self.ends[self.curr_idx];
        let is_rev_opt = self.is_reverse[self.curr_idx];
        self.curr_idx += 1;

        if let (Some(start), Some(end), Some(is_rev)) = (start_opt, end_opt, is_rev_opt) {
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

            if let (Some(start), Some(end)) = (start, end) {
                let seq = self.record.sequence().get(start..=end);
                if let Some(seq) = seq {
                    let mut sequence = Ok(Sequence::from_iter(seq.iter().copied()));
                    if is_rev {
                         sequence = sequence.unwrap().complement().rev().collect::<Result<_, _>>().with_context(||"Could not get the reverse complement of a sequence. Please check if the alphabet is valid.");
                    }
                    Some(sequence)
                } else {
                    Some(Err(anyhow::anyhow!("Could not get the sequence from the start to the end. Please check if the start and end are valid.")))
                }
            } else {
                Some(Err(anyhow::anyhow!("Found invalid start or end. Please check if the start and end are within boundary.")))
            }
        } else {
            Some(Err(anyhow::anyhow!("Found null start, end or strand. Please check if the start, end and strand are valid.")))
        }
    }
}

#[allow(dead_code)]
/// Sorts the indices of a slice based on the slice's values, returning 1-based index order.
///
/// This function calculates the sorted order of indices for a given data slice, adjusting the
/// indices to start from 1 instead of 0 to conform with one-based indexing systems. This can be
/// particularly useful in contexts where array indices are expected to start from 1 (e.g., in
/// certain mathematical or data analysis applications).
///
/// # Type Parameters
///
/// * `T`: The type of the elements in `data`. Must implement the [std::cmp::Ord] trait to enable sorting.
///
/// # Arguments
///
/// * `data`: A slice of data of type `T`. The elements in this slice are used to determine the
///   sorted order of their corresponding indices.
/// * `descending`: A boolean flag indicating the desired sorting order. If `true`, the indices
///   will be sorted according to the corresponding values in descending order. If `false`, the
///   sorting will be in ascending order.
///
/// # Returns
///
/// Returns a [`Vec<usize>`] containing the sorted indices of the slice's elements, starting from 1.
/// For example, if the highest value is at the first position of the input slice, and `descending`
/// is `true`, the first element of the returned vector will be 1.
///
/// # Examples
///
/// ```rust
/// let data = vec![10, 20, 30, 20];
/// let ascending_indices = argsort1based(&data, false);
/// assert_eq!(ascending_indices, vec![1, 2, 4, 3]);
///
/// let descending_indices = argsort1based(&data, true);
/// assert_eq!(descending_indices, vec![3, 2, 4, 1]);
/// ```
///
/// # Note
///
/// The function adjusts for one-based indexing by initializing the indices vector to start
/// from 1 to `data.len()`, thereby aligning with mathematical conventions where arrays are
/// often one-indexed.
pub fn argsort1based<T: Ord>(data: &[T], descending: bool) -> Vec<usize> {
    let mut indices = (1..=data.len()).collect::<Vec<_>>();
    indices.sort_by_key(|&i| &data[i - 1]);
    if descending {
        indices.reverse();
    }
    indices
}

/// Merges genomic intervals based on their start and end positions with specified slack.
///
/// This function takes a Series of structured data containing genomic intervals (start and end positions)
/// and merges intervals that are overlapping or adjacent within a specified slack distance. The function
/// assumes the intervals are sorted by their start positions.
///
/// # Arguments
///
/// * `s`: A [Series] of structured type, expected to contain at least two fields: the start and end positions
///   of genomic intervals. These fields should be of integer type.
///
/// * `slack`: An `i64` value representing the allowed distance between intervals to consider them for merging.
///   If the distance between the end of one interval and the start of the next is less than or equal to this value,
///   the intervals are merged.
///
/// # Returns
///
/// Returns a [`Result<Option<polars::prelude::Column>, PolarsError>`]:
/// * [Ok]`(Some(Series))`: A new `Series` where each element is a merged interval if any merging occurs.
///   The merged intervals are represented as a Series of lists, each containing the start and end of the merged interval.
/// * [Ok]`(None)`: If the input Series is empty or only contains null values.
/// * [Err]`(PolarsError)`: If there are missing values in the start or end columns or other processing errors.
///
/// # Examples
///
/// Assuming `intervals` is a Polars [DataFrame] with a column "intervals" containing structured series
/// with "start" and "end" fields:
///
/// ```rust
/// let merged_intervals = apply_merge(intervals.column("intervals").unwrap().clone(), 10)?;
/// ```
///
/// # Note
///
/// This function requires that the input [Series] is sorted by the start positions of the intervals and contains
/// no null values in the start and end fields. It's designed specifically for genomic data processing where
/// intervals might need to be merged based on their proximity or overlap.
fn apply_merge(s: Column, slack: i64) -> Result<Option<polars::prelude::Column>, PolarsError> {
    // get the two columns from the struct
    let ca: StructChunked = s.struct_()?.clone();

    // get the start and end series
    let start_series = &ca.fields_as_series()[0];
    let end_series = &ca.fields_as_series()[1];

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
            return Result::<Option<polars::prelude::Column>, PolarsError>::Err(
                PolarsError::ComputeError(
                    "Found missing value in the start or end column. Cannot proceed.".into(),
                ),
            );

            // return Result::<Option<polars::prelude::Column>, polars::prelude::PolarsError>::Ok(Some(
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
            return Result::<Option<polars::prelude::Column>, polars::prelude::PolarsError>::Err(
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
                PlSmallStr::from_string(id.to_string()),
                [window_start, window_end],
            ));
            // }
            window_start = curr_start;
            window_end = curr_end;
        }
    }

    // Dont forget the last group
    out_list.push(Series::new("one more".into(), [window_start, window_end]));

    let ls = Column::new("pos".into(), out_list);
    Result::<Option<polars::prelude::Column>, PolarsError>::Ok(Some(ls))
}

/// Identifies gaps between adjacent genomic features based on their start and end positions.
///
/// This function processes a structured [Series] containing genomic intervals and calculates the gaps,
/// i.e., regions that do not overlap with any feature, between these intervals. It is particularly
/// useful for identifying regions like intergenic spaces or introns in genomic datasets.
///
/// # Arguments
///
/// * `s`: A [Series] containing structured data with at least two fields: the start and end positions
///   of genomic features. These features should be sorted by their start positions.
///
/// * `_slack`: Currently unused. Reserved for future use to possibly adjust the definition of gaps based
///   on a slack parameter.
///
/// # Returns
///
/// Returns a [`Result<Option<polars::prelude::Column>, PolarsError>`]:
/// * [Ok]`(Some(Column))`: A new `Column` where each element represents a gap identified between features.
///   The elements are formatted as intervals (start and end positions of the gaps) if any gaps exist.
/// * [Ok]`(None)`: If the input `Column` contains only one feature or is otherwise incapable of forming gaps.
/// * [Err]`(PolarsError)`: If there are missing values in the start or end columns or other issues encountered
///   during processing.
///
/// # Examples
///
/// Assuming `features` is a Polars DataFrame with a column "intervals" containing structured series
/// with "start" and "end" fields representing genomic features:
///
/// ```rust
/// let gaps = apply_gaps(features.column("intervals").unwrap().clone(), 0)?;
/// ```
///
/// # Note
///
/// The function assumes that the input [Series] is sorted by the start positions of the intervals. It is
/// important that there are no null values in the start and end fields for accurate computation. The gaps
/// are defined as regions starting from one feature's end position plus one to the next feature's start
/// position minus one.
// TODO: The implementation is now assuming the intervals are inclusive. This should be changed to be more flexible.
fn apply_gaps(s: Column, _slack: i64) -> Result<Option<polars::prelude::Column>, PolarsError> {
    // get the two columns from the struct
    let ca: StructChunked = s.struct_()?.clone();
    // get the start and end series
    let start_series = &ca.fields_as_series()[0];
    let end_series = &ca.fields_as_series()[1];

    // if we have only one feature, we return an empty list
    if start_series.len() == 1 {
        // return an empty list
        return Result::<Option<polars::prelude::Column>, PolarsError>::Ok(None);

        // return Result::<Option<polars::prelude::Column>, PolarsError>::Ok(Some(
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
            return Result::<Option<polars::prelude::Column>, polars::prelude::PolarsError>::Err(
                polars::prelude::PolarsError::ComputeError(
                    "Found missing value in the start or end column. This should not happen."
                        .into(),
                ),
            );
        };

        out_list.push(Series::new(
            PlSmallStr::from(id.to_string()),
            [gap_start, gap_end],
        ));
    }

    let ls = Column::new("pos".into(), out_list).cast(&DataType::List((DataType::Int64).into()))?;
    Result::<Option<polars::prelude::Column>, PolarsError>::Ok(Some(ls))
}

#[cfg(test)]
mod tests {
    // use polars::prelude::*;
    use super::*;

    use crate::reader::gtf::{AttributeMode, Attributes, GStruct};
    use noodles::core::Position;

    use crate::grangers_utils::FileFormat;

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
            .merge(&["seqname", "gene_id"], false, None, None, true)
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
        let gr1 = gr.merge(&["seqname"], true, None, None, true).unwrap();

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
            .merge(&["seqname", "gene_id"], false, Some(0), None, true)
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
            .merge(&["seqname", "gene_id"], false, Some(2), None, true)
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
        let gr1: Grangers = gr
            .gaps(&["seqname", "gene_id"], false, None, None, true)
            .unwrap();

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
                .str()
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
        let record = noodles::fasta::Record::new(definition, sequence.clone());

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
        let record = noodles::fasta::Record::new(definition, sequence.clone());

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

        let gtf_df = gr.get_gtf_df().unwrap();

        if SAY {
            println!("gtf_df: {:?}", gtf_df);
        }

        let source = vec![String::from("HAVANA"), String::from(".")];
        let gtf_df_attributes = gtf_df
            .column("attributes")
            .unwrap()
            .str()
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
                .str()
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

        // build grangers
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

        // build lapper
        let lappers = gr.build_lappers(true, false, &["gene_id"]).unwrap();

        // In a high level, the hashmap should contains 3 keys:
        // chr1 positive strand
        // chr2 positive strand
        // chr2 negative strand
        let chr1p = lappers.get(&["chr1".to_string(), "+".to_string()]).unwrap();
        if SAY {
            println!("chr1 positive lapper: {:?}", chr1p);
        }

        let chr2p = lappers.get(&["chr2".to_string(), "+".to_string()]).unwrap();
        if SAY {
            println!("chr2 positive lapper: {:?}", chr2p);
        }

        let chr2n = lappers.get(&["chr2".to_string(), "+".to_string()]).unwrap();
        if SAY {
            println!("chr2 negative lapper: {:?}", chr2n);
        }

        // we check if the lappers are built correctly
        let chr1p_o = chr1p.find(11, 15);
        assert!(chr1p_o.count() == 0);

        let chr1p_o = chr1p.find(10, 15);
        assert!(chr1p_o.count() == 1);
    }

    #[test]
    fn test_update_column() {
        let df = df!(
            "seqname" => ["chr1", "chr1", "chr1", "chr2", "chr2", "chr2", "chr2"],
            "feature_type" => ["exon", "exon", "exon", "exon", "exon", "exon", "exon"],
            "start" => [1i64, 21, -5, 1, 51, 1, 51],
            "end" => [10i64, 30, 5, 100, 150, 100, 150],
            "strand"=> ["+", "+", "+", "+", "+", "-", "-"],
            "gene_id" => ["g1", "g1", "g1", "g2", "g2", "g2", "g2"],
        )
        .unwrap();

        // build grangers
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

        // update a non-field column
        let new_col = Column::new("new_col".into(), &[1i64, 2, 3, 4, 5, 6, 7]);
        gr.update_column(new_col.clone(), None).unwrap();
        assert_eq!(gr.column("new_col").unwrap(), &new_col);

        // update an existing field column of the same name
        let gene_id_col = Column::new("gene_id".into(), &["g", "g", "g", "g", "g", "g", "g"]);
        gr.update_column(gene_id_col.clone(), None).unwrap();
        assert_eq!(gr.column("gene_id").unwrap(), &gene_id_col);

        // update an existing field column with a different name
        let gene_id_col = Column::new("gene_id_new".into(), &["g", "g", "g", "g", "g", "g", "g"]);

        gr.update_column(gene_id_col.clone(), Some("gene_id"))
            .unwrap();

        assert_eq!(gr.field_columns().gene_id(), Some("gene_id_new"));

        assert_eq!(
            gr.column(gr.field_columns().gene_id().unwrap()).unwrap(),
            &gene_id_col
        );
    }

    #[test]
    fn test_update_df() {
        let df = df!(
            "seqname" => ["chr1", "chr1", "chr1", "chr2", "chr2", "chr2", "chr2"],
            "feature_type" => ["exon", "exon", "exon", "exon", "exon", "exon", "exon"],
            "start" => [1i64, 21, -5, 1, 51, 1, 51],
            "end" => [10i64, 30, 5, 100, 150, 100, 150],
            "strand"=> ["+", "+", "+", "+", "+", "-", "-"],
            "gene_id" => ["g1", "g1", "g1", "g2", "g2", "g2", "g2"],
        )
        .unwrap();

        // build grangers
        let mut gr: Grangers = Grangers::new(
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

        // first check if the dataframe can be updated
        let df1 = df!(
            "seqname" => ["chr1111", "chr1", "chr1", "chr2", "chr2", "chr2", "chr2"],
            "feature_type" => ["exon", "exon", "exon", "exon", "exon", "exon", "exon"],
            "start" => [1i64, 21, -5, 1, 51, 1, 51],
            "end" => [10i64, 30, 5, 100, 150, 100, 150],
            "strand"=> ["+", "+", "+", "+", "+", "-", "-"],
            "gene_id" => ["g1", "g1", "g1", "g2", "g2", "g2", "g2"],
        )
        .unwrap();

        gr.update_df(df1.clone(), false, false).unwrap();
        assert_eq!(gr.df(), &df1);

        // Then we check if we will get error if the new dataframe is unexpected.
        assert!(gr
            .clone()
            .update_df(DataFrame::default(), false, false)
            .is_err());
        // first check if the dataframe can be updated
        let df2 = df!(
            "seqname" => ["chr1111", "chr1", "chr1", "chr2", "chr2", "chr2", "chr2"],
            "feature_type" => ["exon", "exon", "exon", "exon", "exon", "exon", "exon"],
            "start" => [1i64, 21, -5, 1, 51, 1, 51],
            "end" => [10i64, 30, 5, 100, 150, 100, 150],
            "strand"=> ["+", "+", "+", "+", "+", "-", "-"],
            "gene_iddddddd" => ["g1", "g1", "g1", "g2", "g2", "g2", "g2"],
        )
        .unwrap();

        assert!(gr.update_df(df2, false, false).is_err());
    }
}
