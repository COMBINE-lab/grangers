use crate::grangers_utils::FIELDCOLUMNS;
use anyhow::bail;
use polars::prelude::DataFrame;
use std::collections::HashSet;
use tracing::warn;

/// Represents an inclusive genomic interval.
///
/// This structure is used to define a range on a genome, such as a gene, a regulatory element,
/// or any other genomic feature, where both the start and end positions are considered part of the interval.
///
/// # Fields
///
/// * `start`: The starting position of the interval (1-based index). This is the first position
///   included in the interval.
/// * `end`: The ending position of the interval (1-based index). This position is also included
///   in the interval.
///
/// # Notes
///
/// In bioinformatics, it's crucial to clarify whether intervals are 0-based or 1-based, as well
/// as whether they are inclusive or exclusive. This structure uses 1-based indexing and includes
/// both start and end positions, which is common in formats like GTF or GFF.
///
/// # Examples
///
/// Creating an inclusive interval representing the first 100 bases of a chromosome:
///
/// ```rust
/// let interval = InclusiveInterval { start: 1, end: 100 };
/// assert_eq!(interval.start, 1);
/// assert_eq!(interval.end, 100);
/// ```
///
/// This structure simplifies interval handling and is integral for genomic analyses, ensuring clear
/// and accurate representation of genomic ranges.
pub struct InclusiveInterval {
    pub start: u64,
    pub end: u64,
}

pub enum Strand {
    Positive,
    Negative,
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Strand::Positive => write!(f, "+"),
            Strand::Negative => write!(f, "-"),
        }
    }
}

#[derive(Copy, Clone)]
/// Configuration options for generating flanking regions around genomic intervals.
///
/// This structure is used to specify how flanking regions should be constructed relative to a given
/// genomic interval. Flanking regions can be important for various genomic analyses, including promoter
/// studies, regulatory element identification, and more.
///
/// # Fields
///
/// * `start`: If `true`, generate a flanking region at the start (5' end) of the genomic interval.
/// * `both`: If `true`, generate flanking regions at both ends of the genomic interval. If set to `true`,
///   this option overrides `start` to apply flanking regions to both ends.
/// * `ignore_strand`: If `true`, flanking regions are generated without considering the strand orientation
///   of the genomic interval. This can be useful when the strand information is irrelevant or unavailable.
///
/// # Notes
///
/// Flanking regions are additional sequences that lie adjacent to the main interval of interest. Depending
/// on the options set, these can be added to just one side of the interval, or both. Strand orientation can
/// affect the direction in which flanks are added (5' or 3' ends) unless ignored.
///
/// # Examples
///
/// Creating flank options to get regions only at the start of the interval, considering strand orientation:
///
/// ```rust
/// let flank_options = FlankOptions {
///     start: true,
///     both: false,
///     ignore_strand: false,
/// };
/// ```
///
/// Creating options to get flanking regions on both sides, ignoring strand:
///
/// ```rust
/// let flank_options = FlankOptions {
///     start: false, // Ignored due to `both` being `true`
///     both: true,
///     ignore_strand: true,
/// };
/// ```
///
/// This structure provides a clear and flexible way to define how flanking regions should be generated
/// around genomic intervals for different types of genomic studies.
pub struct FlankOptions {
    pub start: bool,
    pub both: bool,
    pub ignore_strand: bool,
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
    /// Constructs a new `FlankOptions` instance with custom settings.
    ///
    /// Allows the user to specify whether to create flanks at the start, at both ends,
    /// and whether to consider the genomic strand orientation.
    ///
    /// # Arguments
    ///
    /// * `start`: If `true`, a flank will be generated at the start of the interval.
    /// * `both`: If `true`, flanks will be generated on both sides of the interval.
    /// * `ignore_strand`: If `true`, the flanking region will be generated without considering strand orientation.
    ///
    /// # Examples
    ///
    /// Creating a `FlankOptions` instance to generate flanking regions at both ends without strand consideration:
    ///
    /// ```rust
    /// let custom_options = FlankOptions::new(false, true, true);
    /// assert_eq!(custom_options.both, true);
    /// assert_eq!(custom_options.ignore_strand, true);
    /// ```
    ///
    /// This constructor allows for the dynamic creation of flank generation settings, facilitating tailored genomic data manipulation.
    pub fn new(start: bool, both: bool, ignore_strand: bool) -> FlankOptions {
        FlankOptions {
            start,
            both,
            ignore_strand,
        }
    }
}

/// Configuration options for merging genomic intervals.
///
/// This structure is used to define how genomic intervals should be merged based on specific criteria such as sequence name and strand,
/// while considering or ignoring certain attributes like strand orientation. Additionally, it allows setting a slack for merging intervals
/// that are not exactly adjacent but close enough to be considered part of the same feature.
///
/// # Fields
///
/// * `by`: A vector of strings representing the columns by which intervals should be merged.
/// * `slack`: The maximum distance between intervals that can still be considered for merging.
/// * `ignore_strand`: If `true`, intervals are merged without considering their strand orientation.
///
/// # Default
///
/// Implements the `Default` trait, providing a default set of options:
///
/// * `by`: Default to merging by "seqname" and "strand", which is common for genomic intervals.
/// * `slack`: Default slack of 1, allowing adjacent intervals to be merged.
/// * `ignore_strand`: Defaults to `false`, considering strand in merging process.
///
/// # Examples
///
/// Creating default merge options:
///
/// ```rust
/// let default_options = MergeOptions::default();
/// assert_eq!(default_options.by, vec!["seqname", "strand"]);
/// assert_eq!(default_options.slack, 1);
/// assert_eq!(default_options.ignore_strand, false);
/// ```
///
/// Creating custom merge options:
///
/// ```rust
/// let custom_options = MergeOptions::new(&["seqname"], true, 10)?;
/// assert!(custom_options.ignore_strand);
/// assert_eq!(custom_options.slack, 10);
/// ```
pub struct MergeOptions {
    pub by: Vec<String>,
    pub slack: i64,
    pub ignore_strand: bool,
}

impl Default for MergeOptions {
    fn default() -> MergeOptions {
        MergeOptions {
            by: vec![String::from("seqname"), String::from("strand")],
            slack: 1,
            ignore_strand: false,
        }
    }
}

impl MergeOptions {
    /// Creates a new `MergeOptions` instance with specified settings.
    ///
    /// This method initializes `MergeOptions` with a set of columns for merging criteria,
    /// an option to ignore strand information, and a slack for considering nearly adjacent intervals.
    ///
    /// # Arguments
    ///
    /// * `by`: An array of references to strings specifying the columns used for merging.
    ///   It should not include "start" or "end" columns.
    /// * `ignore_strand`: If true, intervals will be merged without regard to their strand orientation.
    /// * `slack`: The distance within which intervals are considered adjacent and thus mergeable.
    ///
    /// # Errors
    ///
    /// Returns an error if the `by` array contains "start" or "end", as these columns cannot be used for merging.
    ///
    /// # Examples
    ///
    /// Creating custom merge options:
    ///
    /// ```rust
    /// let options = MergeOptions::new(&["gene_id", "transcript_id"], false, 5)?;
    /// assert_eq!(options.by, vec!["gene_id", "transcript_id", "strand"]); // Note: "strand" is added by default
    /// assert_eq!(options.slack, 5);
    /// assert!(!options.ignore_strand);
    /// ```
    ///
    /// This method allows the user to configure the merging process to fit specific needs, enhancing flexibility in genomic data handling.
    pub fn new<T: AsRef<str>>(
        by: &[T],
        ignore_strand: bool,
        slack: i64,
    ) -> anyhow::Result<MergeOptions> {
        // avoid duplicated columns
        let mut by_hash: HashSet<String> = by.iter().map(|n| n.as_ref().to_string()).collect();

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
        if by_hash.insert(String::from("seqname")) {
            warn!("Added `seqname` to the `by` vector as it is required.")
        };

        Ok(MergeOptions {
            by: by_hash.into_iter().collect(),
            slack,
            ignore_strand,
        })
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
/// Options for extending a genomic feature.
///
/// This enum is used to specify how a genomic interval should be extended. It can be
/// extended towards the start, the end, or in both directions. This is commonly used
/// in genomic data processing and analysis where extending genomic features like genes
/// or regulatory elements is necessary for various computational tasks.
///
/// # Variants
///
/// * `Start`: Extend the genomic feature towards the start (5' direction).
/// * `End`: Extend the genomic feature towards the end (3' direction).
/// * `Both`: Extend the genomic feature on both sides, towards both start and end.
///
/// # Examples
///
/// Specifying extension towards the start of the feature:
///
/// ```rust
/// let extension_option = ExtendOption::Start;
/// ```
///
/// Specifying extension in both directions:
///
/// ```notrun
/// let extension_option = ExtendOption::Both;
/// ```
///
/// These extension options can be used in various genomic data manipulation tasks, such as
/// extending regulatory regions, adjusting gene coordinates, or creating buffers around features for
/// downstream analysis.
pub enum ExtendOption {
    /// Extend the feature to the start
    Start,
    /// Extend the feature to the end
    End,
    /// Extend the feature to both sides
    Both,
}

pub struct GetSequenceOptions {}

/// Options used for dealing with out-of-boundary features
#[derive(Clone, Copy, PartialEq, Eq)]
/// Options for handling out-of-bounds (OOB) genomic coordinates.
///
/// This enumeration defines strategies for dealing with genomic features or intervals
/// that extend beyond the boundaries of a reference sequence, such as a chromosome or contig.
/// It's commonly used in genomic data processing to determine how to handle intervals
/// when extracting sequences or extending features beyond reference sequence limits.
///
/// # Variants
///
/// * `Truncate`: Adjust out-of-bounds intervals to fit within the available sequence range.
///   This can result in shorter sequences than originally specified but ensures that all
///   returned sequences are valid within the context of the reference.
/// * `Skip`: Ignore any intervals that extend beyond the boundaries of the reference sequence.
///   This can result in some data being omitted from the results but maintains the original
///   size and integrity of the remaining sequences.
///
/// # Examples
///
/// Choosing to truncate sequences that extend beyond the reference:
///
/// ```rust
/// let oob_option = OOBOption::Truncate;
/// ```
///
/// Choosing to skip any features that extend beyond the reference boundaries:
///
/// ```rust
/// let oob_option = OOBOption::Skip;
/// ```
///
/// These handling strategies are vital for ensuring that genomic data analyses remain robust
/// and adaptable to varying data qualities and reference sequence constraints.
pub enum OOBOption {
    Truncate,
    Skip,
}

/// Structure representing the columns of genomic feature data, commonly found in GTF/GFF files.
///
/// This structure contains fields that correspond to standard column names in GTF or GFF format files,
/// facilitating the mapping between the file's columns and the expected fields used in genomic data processing.
/// The fields in this struct are utilized by various Grangers methods to interpret the genomic feature information correctly.
///
/// It is designed according to the GTF/GFF format.
/// <https://useast.ensembl.org/info/website/upload/gff.html>
/// The default values are the column names used in GTF/GFF files.
/// This will be used in almost all Grangers methods.
/// # Fields
///
/// - `seqname`: Corresponds to the reference sequence name column in GTF/GFF (required).
/// - `source`: Source of the feature, such as database or algorithm used (optional).
/// - `feature_type`: Type of genomic feature (e.g., gene, transcript, exon) (optional).
/// - `start`: Start position of the feature in the reference sequence (required).
/// - `end`: End position of the feature in the reference sequence (required).
/// - `score`: Score or value associated with the feature (optional).
/// - `strand`: Strand of the feature ('+' or '-') (required).
/// - `phase`: Phase or frame for coding features (optional).
/// - `gene_id`: Identifier for the gene associated with the feature (optional).
/// - `gene_name`: Name of the gene associated with the feature (optional).
/// - `transcript_id`: Identifier for the transcript associated with the feature (optional).
/// - `exon_number`: Order or number of the exon within its transcript (optional).
///
/// # Examples
///
/// Creating default `FieldColumns` for standard GTF/GFF processing:
///
/// ```rust
/// let field_columns = FieldColumns::default();
/// ```
///
/// Customizing `FieldColumns` for specific data processing needs:
///
/// ```rust
/// let mut custom_fields = FieldColumns::default();
/// custom_fields.gene_id = Some("gene_ident".to_string());
/// custom_fields.exon_number = Some("exon_order".to_string());
/// ```
///
/// The default implementations adhere to typical conventions, but they can be adjusted as needed
/// to fit the specific layout and requirements of the input data files.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct FieldColumns {
    /// the name of the reference sequence column, usually it is called "seqname".\
    /// This corresponds to the first column in a GTF/GFF file.\
    /// This field is required, and the corresponding column should not contain missing values (nulls).\
    /// You can drop all rows with missing values in this column by calling `df.drop_nulls(Some(["seqname"]))`.
    pub seqname: String,
    /// the name of the source column, usually it is called "source".
    /// This is the second column in a GTF/GFF file. It records the source (HAVANA, ENSEMBL, etc) of the feature.
    pub source: Option<String>,
    /// The name of the feature type column whose values should be "gene", "transcript" and "exon", etc.
    /// This should be the third column in a GTF/GFF file.
    /// This field is required for transcriptome related methods, like `introns()`, `get_transcript_sequeneces()`, etc.
    pub feature_type: Option<String>,
    /// The column name of the start position column, usually it is called "start".
    /// This is the fourth column in a GTF/GFF file.
    /// This field is required, and the corresponding column should not contain missing values (nulls).
    pub start: String,
    /// The column name of the end position column, usually it is called "end".
    /// This is the fifth column in a GTF/GFF file.
    /// This field is required, and the corresponding column should not contain missing values (nulls).
    pub end: String,
    /// The column name of the score column, usually it is called "score".
    pub score: Option<String>,
    /// The column name of the strand column, usually it is called "strand".
    /// This is the seventh column in a GTF/GFF file.\
    /// This field is required, and the corresponding column should not contain missing values (nulls) for calling most Grangers methods.\
    /// If this field is missing, you can add one by calling `df.update_column("strand", vec!['.'; df.height()])`.\
    /// If it contains missing values, you can fill them by calling `df.fill_none("strand", '.')`.
    pub strand: String,
    /// The column name of the phase column, usually it is called "frame" or "phase".
    /// This is the eighth column in a GTF/GFF file.
    pub phase: Option<String>,
    /// The column name of the gene_id column, usually it is called "gene_id".
    pub gene_id: Option<String>,
    /// The column name of the gene_name column, usually it is called "gene_id".
    pub gene_name: Option<String>,
    /// The column name of the transcript ID column, usually it is called "transcript_id".
    pub transcript_id: Option<String>,
    // The column name of the exon ID column, usually it is called "exon_id".
    // pub exon_id: Option<String>,
    /// The column name of the exon number (order) column, usually it is called "exon_number".
    /// This column is used to sort the exons of a transcript.
    /// If this column is missing, the exons will be sorted by their start positions.
    pub exon_number: Option<String>,
}

impl FieldColumns {
    /// Returns a reference to the `seqname` field.
    ///
    /// This field represents the name of the reference sequence (such as a chromosome or contig) to which the genomic features pertain.
    /// It's a crucial identifier in genomic datasets, typically corresponding to the first column in GTF/GFF files.
    ///
    /// # Returns
    /// A string slice (&str) pointing to the `seqname` value.
    pub fn seqname(&self) -> &str {
        self.seqname.as_str()
    }

    /// Returns an optional reference to the `source` field.
    ///
    /// This field denotes the origin of the genomic feature, such as the software or organization that generated the data.
    /// It aligns with the second column in standard GTF/GFF formats, though it's not mandatory for all processing tasks.
    ///
    /// # Returns
    /// An optional string slice (&str) pointing to the `source` value, if it exists.
    pub fn source(&self) -> Option<&str> {
        self.source.as_deref()
    }

    /// Returns an optional reference to the `feature_type` field.
    ///
    /// This specifies the type of genomic feature (e.g., gene, transcript, exon) and corresponds to the third column in GTF/GFF files.
    /// It's vital for differentiating between various biological entities within the genomic data.
    ///
    /// # Returns
    /// An optional string slice (&str) pointing to the `feature_type` value, if it exists.
    pub fn feature_type(&self) -> Option<&str> {
        self.feature_type.as_deref()
    }
    /// Returns a reference to the `start` field.
    ///
    /// This indicates the starting position of the genomic feature on the reference sequence, typically aligning with the fourth column in GTF/GFF files.
    /// It's essential for identifying the genomic location of the feature.
    ///
    /// # Returns
    /// A string slice (&str) pointing to the `start` value.
    pub fn start(&self) -> &str {
        self.start.as_str()
    }
    /// Returns a reference to the `end` field.
    ///
    /// This marks the ending position of the genomic feature, typically aligning with the fifth column in GTF/GFF formats.
    /// This value, along with `start`, helps define the exact genomic span of the feature.
    ///
    /// # Returns
    /// A string slice (&str) pointing to the `end` value.
    pub fn end(&self) -> &str {
        self.end.as_str()
    }
    /// Returns an optional reference to the `score` field.
    ///
    /// This field may contain a numeric score or value associated with the feature's significance or confidence level, corresponding to the sixth column in some formats.
    ///
    /// # Returns
    /// An optional string slice (&str) pointing to the `score` value, if it exists.
    pub fn score(&self) -> Option<&str> {
        self.score.as_deref()
    }
    /// Returns a reference to the `strand` field.
    ///
    /// This indicates the genomic strand (either '+' or '-') that the feature is associated with, typically found in the seventh column in GTF/GFF files.
    /// It's crucial for understanding the directional context of the feature.
    ///
    /// # Returns
    /// A string slice (&str) pointing to the `strand` value.
    pub fn strand(&self) -> &str {
        self.strand.as_str()
    }
    /// Returns an optional reference to the `phase` field.
    ///
    /// This field, often called 'frame', is relevant for coding sequences and can affect how the sequence is translated into amino acids.
    ///
    /// # Returns
    /// An optional string slice (&str) pointing to the `phase` value, if it exists.
    pub fn phase(&self) -> Option<&str> {
        self.phase.as_deref()
    }
    /// Returns an optional reference to the `gene_id` field.
    ///
    /// This field is typically used to associate features with a specific gene identifier, facilitating gene-centric analyses.
    ///
    /// # Returns
    /// An optional string slice (&str) pointing to the `gene_id` value, if it exists.
    pub fn gene_id(&self) -> Option<&str> {
        self.gene_id.as_deref()
    }

    /// Returns an optional reference to the `gene_name` field.
    ///
    /// This is akin to `gene_id` but typically contains a human-readable gene name, useful for reports and labeling.
    ///
    /// # Returns
    /// An optional string slice (&str) pointing to the `gene_name` value, if it exists.
    pub fn gene_name(&self) -> Option<&str> {
        self.gene_name.as_deref()
    }

    /// Returns an optional reference to the `transcript_id` field.
    ///
    /// This field links genomic features to a specific transcript identifier, essential for transcriptome analyses.
    ///
    /// # Returns
    /// An optional string slice (&str) pointing to the `transcript_id` value, if it exists.
    pub fn transcript_id(&self) -> Option<&str> {
        self.transcript_id.as_deref()
    }
    // get a reference to the exon_id field
    // pub fn exon_id(&self) -> Option<&str> {
    //     self.exon_id.as_deref()
    // }

    /// Returns an optional reference to the `exon_number` field.
    ///
    /// This field is used to specify the ordering of exons within a transcript. It is particularly
    /// useful for reconstructing the transcript structure from exon features.
    ///
    /// # Examples
    ///
    /// ```
    /// let field_columns = FieldColumns::default();
    /// assert!(field_columns.exon_number().is_none());
    ///
    /// let mut field_columns = FieldColumns::default();
    /// field_columns.update("exon_number", "my_exon_number");
    /// assert_eq!(field_columns.exon_number(), Some("my_exon_number"));
    /// ```
    ///
    /// # Returns
    ///
    /// An optional string slice (`&str`) pointing to the `exon_number` value, if it exists.
    /// This aids in sorting exons to align with their sequential order in the corresponding transcript.
    pub fn exon_number(&self) -> Option<&str> {
        self.exon_number.as_deref()
    }
}

impl Default for FieldColumns {
    fn default() -> Self {
        Self {
            seqname: "seqname".to_string(),
            source: Some("source".to_string()),
            feature_type: Some("feature_type".to_string()),
            start: "start".to_string(),
            end: "end".to_string(),
            score: Some("score".to_string()),
            strand: "strand".to_string(),
            phase: Some("phase".to_string()),
            gene_id: Some("gene_id".to_string()),
            gene_name: Some("gene_name".to_string()),
            transcript_id: Some("transcript_id".to_string()),
            exon_number: Some("exon_number".to_string()),
        }
    }
}

impl FieldColumns {
    /// Returns an array of optional field names present in `FieldColumns`.
    ///
    /// This method provides the names of optional fields that may not always be
    /// required for every analysis but can provide additional context or data
    /// when available. These fields include `source`, `feature_type`, `score`,
    /// `phase`, `gene_id`, `gene_name`, `transcript_id`, and `exon_number`.
    ///
    /// # Returns
    ///
    /// An array of `Option<&str>`, where each element represents the name of an
    /// optional field if it is set, otherwise `None`.
    ///
    /// # Examples
    ///
    /// ```
    /// let field_columns = FieldColumns::default();
    /// let optional_fields = field_columns.optional_fields();
    /// println!("{:?}", optional_fields);
    /// ```
    pub fn optional_fields(&self) -> [Option<&str>; 8] {
        [
            self.source(),
            self.feature_type(),
            self.score(),
            self.phase(),
            self.gene_id(),
            self.gene_name(),
            self.transcript_id(),
            self.exon_number(),
        ]
    }
    /// Returns an array of field names used specifically in the GTF file format.
    ///
    /// This method extracts the names of fields that are standard for GTF format files.
    /// These include `seqname`, `source`, `feature_type`, `start`, `end`, `score`,
    /// `strand`, and `phase`. Some of these fields might be optional in general but
    /// are typically found in GTF files.
    ///
    /// # Returns
    ///
    /// An array of `&str`, each representing a field name used in GTF files.
    ///
    /// # Examples
    ///
    /// ```
    /// let field_columns = FieldColumns::default();
    /// let gtf_fields = field_columns.gtf_fields();
    /// println!("{:?}", gtf_fields);
    /// ```
    pub fn gtf_fields(&self) -> [&str; 8] {
        [
            self.seqname(),
            self.source().unwrap_or(""),
            self.feature_type().unwrap_or(""),
            self.start(),
            self.end(),
            self.score().unwrap_or(""),
            self.strand(),
            self.phase().unwrap_or(""),
        ]
    }

    /// Returns an array of attribute field names used specifically in the GTF format.
    ///
    /// In the context of GTF files, certain attributes provide additional information
    /// about genomic features such as `gene_id`, `gene_name`, `transcript_id`, and
    /// `exon_number`. This method returns the names of these attribute fields if they
    /// have been set.
    ///
    /// # Returns
    ///
    /// An array of `Option<&str>`, where each element corresponds to one of the GTF
    /// attribute fields. If a particular attribute is not set, its value in the array
    /// will be `None`.
    ///
    /// # Examples
    ///
    /// ```
    /// let field_columns = FieldColumns::default();
    /// let gtf_attributes = field_columns.gtf_attributes();
    /// println!("{:?}", gtf_attributes);
    /// ```
    pub fn gtf_attributes(&self) -> [Option<&str>; 4] {
        [
            self.gene_id(),
            self.gene_name(),
            self.transcript_id(),
            self.exon_number(),
        ]
    }
    /// Returns an array of essential field names required for basic genomic feature analysis.
    ///
    /// This method provides the names of fields that are fundamental to genomic analyses,
    /// specifically those related to the positioning and identification of genomic features.
    /// These fields include `seqname`, `start`, `end`, and `strand`.
    ///
    /// # Returns
    ///
    /// An array of `&str`, each representing a name of an essential field.
    ///
    /// # Examples
    ///
    /// ```
    /// let field_columns = FieldColumns::default();
    /// let essential_fields = field_columns.essential_fields();
    /// println!("{:?}", essential_fields);
    /// ```
    pub fn essential_fields(&self) -> [&str; 4] {
        [self.seqname(), self.start(), self.end(), self.strand()]
    }

    /// Validates the `FieldColumns` against a given DataFrame.
    ///
    /// This method checks if the required fields specified in the `FieldColumns` exist in
    /// the provided DataFrame. It optionally prints warnings and can halt execution if
    /// critical fields are missing and cannot be automatically fixed.
    ///
    /// # Arguments
    ///
    /// * `df` - The DataFrame against which to validate the `FieldColumns`.
    /// * `is_warn` - If `true`, print warnings for missing fields.
    /// * `is_bail` - If `true`, throw an error if essential fields are missing.
    ///
    /// # Returns
    ///
    /// `Ok(bool)` indicating whether the `FieldColumns` is valid within the context of
    /// the provided DataFrame. Returns `Err` if critical fields are missing and cannot be fixed.
    ///
    /// # Examples
    ///
    /// ```
    /// let mut field_columns = FieldColumns::default();
    /// let df = DataFrame::new(vec![])?;
    /// let is_valid = field_columns.is_valid(&df, true, false)?;
    /// println!("Is valid: {}", is_valid);
    /// ```
    pub fn is_valid(&self, df: &DataFrame, is_warn: bool, is_bail: bool) -> anyhow::Result<bool> {
        let mut is_valid = true;
        // check required fields
        if df.column(self.seqname()).is_err() {
            is_valid = false;
            if is_warn {
                warn!(
                    "The dataframe does not contain the specified seqname column {}; Cannot proceed. You can add one by calling `df.update_column(\"seqname\", vec!['.'; df.height()])`",
                    self.seqname()
                )
            }
        }
        if df.column(self.start()).is_err() {
            is_valid = false;
            if is_warn {
                warn!(
                    "The dataframe does not contain the specified start column {}; Cannot proceed. You can add one by calling `df.update_column(\"start\", vec!['.'; df.height()])`",
                    self.start()
                )
            }
        }
        if df.column(self.end()).is_err() {
            is_valid = false;
            if is_warn {
                warn!(
                    "The dataframe does not contain the specified end column {}; Cannot proceed. You can add one by calling `df.update_column(\"end\", vec!['.'; df.height()])`",
                    self.end()
                )
            }
        }
        if df.column(self.strand()).is_err() {
            is_valid = false;
            if is_warn {
                warn!(
                    "The dataframe does not contain the specified strand column {}; Cannot proceed. You can add one by calling `df.update_column(\"strand\", vec!['.'; df.height()])`",
                    self.strand()
                )
            }
        }
        // check additional fields
        if let Some(s) = self.source() {
            if df.column(s).is_err() {
                is_valid = false;
                if is_warn {
                    warn!("The provided source column {} is not found in the dataframe; It will be ignored", s)
                }
            }
        }
        if let Some(s) = self.feature_type() {
            if df.column(s).is_err() {
                is_valid = false;
                if is_warn {
                    warn!("The provided feature_type column {} is not found in the dataframe; It will be ignored", s)
                }
            }
        }
        if let Some(s) = self.score() {
            if df.column(s).is_err() {
                is_valid = false;
                if is_warn {
                    warn!("The provided score column {} is not found in the dataframe; It will be ignored", s)
                }
            }
        }
        if let Some(s) = self.phase() {
            if df.column(s).is_err() {
                is_valid = false;
                if is_warn {
                    warn!("The provided phase column {} is not found in the dataframe; It will be ignored", s)
                }
            }
        }

        if let Some(s) = self.gene_id() {
            if df.column(s).is_err() {
                is_valid = false;
                if is_warn {
                    warn!("The provided gene_id column {} is not found in the dataframe; It will be ignored", s)
                }
            }
        }

        if let Some(s) = self.gene_name() {
            if df.column(s).is_err() {
                is_valid = false;
                if is_warn {
                    warn!("The provided gene_name column {} is not found in the dataframe; It will be ignored", s)
                }
            }
        }

        if let Some(s) = self.transcript_id() {
            if df.column(s).is_err() {
                is_valid = false;
                if is_warn {
                    warn!("The provided transcript_id column {} is not found in the dataframe; It will be ignored", s)
                }
            }
        }

        if let Some(s) = self.exon_number() {
            if df.column(s).is_err() {
                is_valid = false;
                if is_warn {
                    warn!("The provided exon_number column {} is not found in the dataframe; It will be ignored", s)
                }
            }
        }

        if !is_valid & is_bail {
            bail!(
                "The FieldColumns is not valid; Please try fix it by calling FieldColumns::fix()."
            )
        }

        if !is_valid & is_warn {
            warn!(
                "The FieldColumns is not valid; Please try fix it by calling FieldColumns::fix()."
            )
        }

        Ok(is_valid)
    }

    /// Attempts to correct any missing or incorrect field mappings based on a provided DataFrame.
    ///
    /// This method attempts to fix the FieldColumns instance by ensuring that all required and optional
    /// fields have corresponding columns in the given DataFrame. If a necessary field is missing or
    /// incorrect in FieldColumns but exists in the DataFrame under a standard name (like "seqname" for
    /// the sequence name), this method updates the FieldColumns instance to use the correct column name.
    /// If a standard column name does not exist in the DataFrame, an error is returned.
    ///
    /// # Arguments
    ///
    /// * `df` - The DataFrame against which to validate and fix the FieldColumns.
    /// * `is_warn` - If true, the method will log a warning message for each issue it encounters and attempts to fix.
    ///
    /// # Returns
    ///
    /// `Ok(())` if the FieldColumns instance was successfully fixed or if no fixes were needed.
    /// `Err(anyhow::Error)` if a required field could not be fixed because the corresponding column
    /// does not exist in the DataFrame.
    ///
    /// # Examples
    ///
    /// ```
    /// let mut field_columns = FieldColumns::default();
    /// let df = DataFrame::new(vec![])?; // Assume this is a populated DataFrame
    /// field_columns.fix(&df, true)?;
    /// ```
    pub fn fix(&mut self, df: &DataFrame, is_warn: bool) -> anyhow::Result<()> {
        // try fix required fields
        if df.column(self.seqname()).is_err() {
            if is_warn {
                warn!(
                    "cannot find the specified seqname column {} in the dataframe; try to fix",
                    self.seqname()
                );
            }
            if df.column("seqname").is_ok() {
                self.seqname = "seqname".to_string();
            } else {
                bail!("The dataframe does not contain the specified seqname column {} or a column named \"seqname\"; Cannot fix.", self.seqname());
            }
        }
        if df.column(self.start()).is_err() {
            if is_warn {
                warn!(
                    "cannot find the specified start column {} in the dataframe; try to fix",
                    self.start()
                );
            }
            if df.column("start").is_ok() {
                self.start = "start".to_string();
            } else {
                bail!("The dataframe does not contain the specified start column {} or a column named \"start\"; Cannot fix.", self.start());
            }
        }
        if df.column(self.end()).is_err() {
            if is_warn {
                warn!(
                    "cannot find the specified end column {} in the dataframe; try to fix",
                    self.end()
                );
            }
            if df.column("end").is_ok() {
                self.end = "end".to_string();
            } else {
                bail!("The dataframe does not contain the specified end column {} or a column named \"end\"; Cannot fix.", self.end());
            }
        }
        if df.column(self.strand()).is_err() {
            if is_warn {
                warn!(
                    "cannot find the specified strand column {} in the dataframe; try to fix",
                    self.strand()
                );
            }
            if df.column("strand").is_ok() {
                self.strand = "strand".to_string();
            } else {
                bail!("The dataframe does not contain the specified strand column {} or a column named \"strand\"; Cannot fix. If this is desired, you can add a dummy strand column by calling `df.update_column(\"strand\", vec!['.'; df.height()])`", self.strand());
            }
        }

        // try fix optional fields
        if let Some(s) = self.source() {
            if df.column(s).is_err() {
                if is_warn {
                    warn!(
                        "cannot find the specified source column {} in the dataframe; try to fix",
                        s
                    );
                }
                self.source = if df.column("source").is_ok() {
                    Some("source".to_string())
                } else {
                    None
                }
            }
        }
        if let Some(s) = self.feature_type() {
            if df.column(s).is_err() {
                if is_warn {
                    warn!("cannot find the specified feature_type column {} in the dataframe; try to fix", s);
                }
                self.feature_type = if df.column("feature_type").is_ok() {
                    Some("feature_type".to_string())
                } else {
                    None
                }
            }
        }
        if let Some(s) = self.score() {
            if df.column(s).is_err() {
                if is_warn {
                    warn!(
                        "cannot find the specified score column {} in the dataframe; try to fix",
                        s
                    );
                }
                self.score = if df.column("score").is_ok() {
                    Some("score".to_string())
                } else {
                    None
                }
            }
        }
        if let Some(s) = self.phase() {
            if df.column(s).is_err() {
                if is_warn {
                    warn!(
                        "cannot find the specified phase column {} in the dataframe; try to fix",
                        s
                    );
                }
                self.phase = if df.column("phase").is_ok() {
                    Some("phase".to_string())
                } else {
                    None
                }
            }
        }
        if let Some(s) = self.gene_id() {
            if df.column(s).is_err() {
                if is_warn {
                    warn!(
                        "cannot find the specified gene_id column {} in the dataframe; try to fix",
                        s
                    );
                }
                self.gene_id = if df.column("gene_id").is_ok() {
                    Some("gene_id".to_string())
                } else {
                    None
                }
            }
        }

        if let Some(s) = self.gene_name() {
            if df.column(s).is_err() {
                if is_warn {
                    warn!(
                        "cannot find the specified gene_name column {} in the dataframe; try to fix",
                        s
                    );
                }
                self.gene_name = if df.column("gene_name").is_ok() {
                    Some("gene_name".to_string())
                } else {
                    None
                }
            }
        }

        if let Some(s) = self.transcript_id() {
            if df.column(s).is_err() {
                if is_warn {
                    warn!("cannot find the specified transcript_id column {} in the dataframe; try to fix", s);
                }
                self.transcript_id = if df.column("transcript_id").is_ok() {
                    Some("transcript_id".to_string())
                } else {
                    None
                }
            }
        }

        if let Some(s) = self.exon_number() {
            if df.column(s).is_err() {
                if is_warn {
                    warn!("cannot find the specified exon_number column {} in the dataframe; try to fix", s);
                }
                self.exon_number = if df.column("exon_number").is_ok() {
                    Some("exon_number".to_string())
                } else {
                    None
                }
            }
        }

        Ok(())
    }

    /// Updates the value of a specific field within the `FieldColumns`.
    ///
    /// This method allows modifying the column names mapped in `FieldColumns` to better match the
    /// structure of a different DataFrame or to correct any mistakes in the initial setup. If the
    /// specified field name does not exist in `FieldColumns`, the method will return an error.
    ///
    /// # Arguments
    ///
    /// * `field` - A string slice that holds the name of the field to update (e.g., "seqname", "start").
    /// * `value` - A string slice representing the new value for the field.
    ///
    /// # Returns
    ///
    /// * `Ok(())` if the field was successfully updated.
    /// * `Err(anyhow::Error)` if the field name provided does not exist in `FieldColumns`.
    ///
    /// # Examples
    ///
    /// ```
    /// let mut field_columns = FieldColumns::default();
    /// field_columns.update("seqname", "chromosome")?;
    /// assert_eq!(field_columns.seqname(), "chromosome");
    /// ```
    pub fn update<T: AsRef<str>>(&mut self, field: T, value: T) -> anyhow::Result<()> {
        let value = value.as_ref().to_string();
        match field.as_ref() {
            "seqname" => self.seqname = value,
            "source" => self.source = Some(value),
            "feature_type" => self.feature_type = Some(value),
            "start" => self.start = value,
            "end" => self.end = value,
            "score" => self.score = Some(value),
            "strand" => self.strand = value,
            "phase" => self.phase = Some(value),
            "gene_id" => self.gene_id = Some(value),
            "gene_name" => self.gene_name = Some(value),
            "transcript_id" => self.transcript_id = Some(value),
            "exon_number" => self.exon_number = Some(value),
            _ => bail!("invalid field name: {}", field.as_ref()),
        }

        Ok(())
    }

    /// Retrieves the value of a specified field from the `FieldColumns`.
    ///
    /// This method returns the current value associated with a given field, if it exists.
    ///
    /// # Arguments
    ///
    /// * `field` - A string slice representing the name of the field to retrieve.
    ///
    /// # Returns
    ///
    /// An `Option<&str>` which is `Some` containing the value of the field if the field exists,
    /// otherwise `None`.
    ///
    /// # Examples
    ///
    /// ```
    /// let field_columns = FieldColumns::default();
    /// assert_eq!(field_columns.field("seqname"), Some("seqname"));
    /// assert_eq!(field_columns.field("nonexistent"), None);
    /// ```
    pub fn field<T: AsRef<str>>(&self, field: T) -> Option<&str> {
        match field.as_ref() {
            "seqname" => Some(self.seqname.as_str()),
            "source" => self.source(),
            "feature_type" => self.feature_type(),
            "start" => Some(self.start.as_str()),
            "end" => Some(self.end.as_str()),
            "score" => self.score(),
            "strand" => Some(self.strand.as_str()),
            "phase" => self.phase(),
            "gene_id" => self.gene_id(),
            "gene_name" => self.gene_name(),
            "transcript_id" => self.transcript_id(),
            "exon_number" => self.exon_number(),
            _ => None,
        }
    }

    /// Retrieves the value of a specified field from the `FieldColumns`, with an option to trigger an error if the field is not found.
    ///
    /// This method functions similarly to `field`, but it can optionally return an error if the field does not exist, based on the `is_bail` argument.
    ///
    /// # Arguments
    ///
    /// * `field` - A string slice representing the name of the field to retrieve.
    /// * `is_bail` - A boolean indicating whether to bail (return an error) if the field does not exist.
    ///
    /// # Returns
    ///
    /// * `Ok(Some(&str))` if the field exists.
    /// * `Ok(None)` if the field does not exist and `is_bail` is set to `false`.
    /// * `Err(anyhow::Error)` if the field does not exist and `is_bail` is set to `true`.
    ///
    /// # Examples
    ///
    /// ```
    /// let field_columns = FieldColumns::default();
    /// assert_eq!(field_columns.field_checked("seqname", false).unwrap(), Some("seqname"));
    /// assert!(field_columns.field_checked("nonexistent", true).is_err());
    /// ```
    pub fn field_checked<T: AsRef<str>>(
        &self,
        field: T,
        is_bail: bool,
    ) -> anyhow::Result<Option<&str>> {
        match field.as_ref() {
            "seqname" => Ok(Some(self.seqname.as_str())),
            "source" => Ok(self.source()),
            "feature_type" => Ok(self.feature_type()),
            "start" => Ok(Some(self.start.as_str())),
            "end" => Ok(Some(self.end.as_str())),
            "score" => Ok(self.score()),
            "strand" => Ok(Some(self.strand.as_str())),
            "phase" => Ok(self.phase()),
            "gene_id" => Ok(self.gene_id()),
            "gene_name" => Ok(self.gene_name()),
            "transcript_id" => Ok(self.transcript_id()),
            "exon_number" => Ok(self.exon_number()),
            _ => {
                if is_bail {
                    bail!(
                        "The provided field {} is not a valid field name; Cannot proceed",
                        field.as_ref()
                    );
                }

                warn!(
                    "The provided field {} is not a valid field name; It will be ignored",
                    field.as_ref()
                );
                Ok(None)
            }
        }
    }

    pub fn all_fields() -> [&'static str; 12] {
        FIELDCOLUMNS
    }
}
