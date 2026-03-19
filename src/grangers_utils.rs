// use crate::grangers_info::{Grangers, GrangersSequenceCollection};
use flate2::bufread::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use tracing::trace;
use noodles::fasta::io::Reader as NoodlesFastaReader;

/// Type alias for a noodles FASTA reader that can read from
/// a `dyn BufRead`. It is used to allow reading from either
/// a compressed or uncompressed FASTA file.
pub type FastaReader = NoodlesFastaReader<Box<dyn BufRead>>;

pub(crate) const VALIDSTRANDS: [&str; 2] = ["+", "-"];

#[derive(Copy, Clone, PartialEq, Eq)]
/// Represents supported genomic file formats.
///
/// This enumeration includes various file formats commonly used in bioinformatics and genomic data processing,
/// such as GTF, GFF, BED, SAM, BAM, FASTA, and FASTQ.
///
/// # Variants
///
/// * `GTF` - Gene Transfer Format, commonly used for annotation and gene information.
/// * `GFF` - General Feature Format, used for describing genes and other features of DNA, RNA, and protein sequences.
/// * `BED` - Browser Extensible Data, a format for describing genomic regions.
/// * `SAM` - Sequence Alignment/Map format, used for storing sequence alignments.
/// * `BAM` - Binary version of the SAM format, more efficient for storage and analysis.
/// * `FASTA` - Text-based format for representing nucleotide or peptide sequences.
/// * `FASTQ` - Similar to FASTA but also includes quality information for each sequence.
///
/// # Methods
///
/// * `get_essential`: Returns a slice of essential attribute names specific to GTF or GFF formats.
/// * `is_gtf`: Checks if the file format is GTF.
///
pub enum FileFormat {
    GTF,
    GFF,
    BED,
    SAM,
    BAM,
    FASTA,
    FASTQ,
}

impl FileFormat {
    /// Retrieves the essential attributes associated with the file format.
    ///
    /// This method provides a convenient way to access a predefined list of essential attributes
    /// for genomic file formats that have structured attributes, such as GTF and GFF. These attributes
    /// are crucial for the correct parsing and handling of these file formats in genomic analyses.
    ///
    /// # Returns
    ///
    /// A slice of string slices (&[`[&str]`]):
    /// * For [FileFormat::GTF] and [FileFormat::GFF], it returns a reference to an array containing the names
    ///   of essential attributes defined in `GXFESSENTIALATTRIBUTES`.
    /// * For all other file formats, it returns an empty slice, as they do not have a predefined set of
    ///   essential attributes.
    ///
    /// # Examples
    ///
    /// Getting essential attributes for GTF and GFF formats:
    ///
    /// ```rust
    /// let gtf_format = FileFormat::GTF;
    /// let gtf_essentials = gtf_format.get_essential();
    /// assert_eq!(gtf_essentials, ["gene_id", "transcript_id"]); // Assuming these are the essential attributes for GTF/GFF.
    ///
    /// let bed_format = FileFormat::BED;
    /// let bed_essentials = bed_format.get_essential();
    /// assert!(bed_essentials.is_empty());
    /// ```
    ///
    /// This method is particularly useful when working with genomic data parsers or validators that require
    /// knowledge of essential fields for correct processing.
    pub fn get_essential(&self) -> &[&str] {
        match self {
            FileFormat::GTF => GXFESSENTIALATTRIBUTES.as_ref(),
            FileFormat::GFF => GXFESSENTIALATTRIBUTES.as_ref(),
            _ => &[],
        }
    }

    /// Checks if the file format is GTF (Gene Transfer Format).
    ///
    /// This method allows for a quick verification to determine whether the current instance
    /// of the [FileFormat] enum is set to the GTF format, which is commonly used for gene annotations.
    ///
    /// # Returns
    ///
    /// * `true`: If the [FileFormat] instance is [FileFormat::GTF].
    /// * `false`: Otherwise.
    ///
    /// # Examples
    ///
    /// Checking the file format:
    ///
    /// ```rust
    /// let format = FileFormat::GTF;
    /// assert!(format.is_gtf()); // Returns true because the format is GTF.
    ///
    /// let another_format = FileFormat::FASTA;
    /// assert!(!another_format.is_gtf()); // Returns false because the format is not GTF.
    /// ```
    ///
    /// This method is particularly useful for conditional processing based on the file format, such as applying
    /// GTF-specific parsing logic or validation checks.
    pub fn is_gtf(&self) -> bool {
        matches!(self, FileFormat::GTF)
    }
}

impl std::str::FromStr for FileFormat {
    type Err = anyhow::Error;

    /// Converts from a [&str] to an appropriate [FileFormat] type.
    /// The result is returned in an [`anyhow::Result<FileFormat>`]
    /// and is an error variant if there is no corresponding type for
    /// the input argument `s`.
    fn from_str(s: &str) -> anyhow::Result<FileFormat> {
        let ft = match s.to_lowercase().as_str() {
            "gtf" => FileFormat::GTF,
            "gff2" => FileFormat::GTF,
            "gff" => FileFormat::GFF,
            "gff3" => FileFormat::GFF,
            "bed" => FileFormat::BED,
            "sam" => FileFormat::SAM,
            "bam" => FileFormat::BAM,
            "fasta" => FileFormat::FASTA,
            "fastq" => FileFormat::FASTQ,
            _ => anyhow::bail!("Cannot parse the file type."),
        };
        Ok(ft)
    }
}

impl std::fmt::Display for FileFormat {
    /// Print the formatted description of the current [FileFormat]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FileFormat::GTF => write!(f, "GTF"),
            FileFormat::GFF => write!(f, "GFF"),
            FileFormat::BED => write!(f, "BED"),
            FileFormat::SAM => write!(f, "SAM"),
            FileFormat::BAM => write!(f, "BAM"),
            FileFormat::FASTA => write!(f, "FASTA"),
            FileFormat::FASTQ => write!(f, "FASTQ"),
        }
    }
}

pub(crate) const GXFESSENTIALATTRIBUTES: [&str; 4] =
    ["gene_id", "gene_name", "transcript_id", "exon_number"];

pub(crate) const GXFFIELDS: [&str; 8] = [
    "seqname",
    "source",
    "feature_type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
];

/// traverse the given file to get the number of lines in the file
pub fn _file_line_count<T: AsRef<Path>>(file_path: T) -> anyhow::Result<usize> {
    let reader = BufReader::new(File::open(file_path)?);
    let mut num_lines = 0usize;
    for l in reader.lines() {
        let line = l?;
        if !(line.trim().starts_with('#') | line.trim().is_empty()) {
            num_lines += 1;
        }
    }
    Ok(num_lines)
}

// Returns `true` if the input vectors are of equal length and false otherwise.
pub fn equal_length<T, R>(vec1: &[T], vec2: &[R]) -> bool {
    vec1.len() == vec2.len()
}

pub fn _setdiff<T: Eq + Clone>(vec1: &[T], vec2: &[T]) -> Vec<T> {
    let mut diff: Vec<T> = Vec::new();
    for v in vec1.iter() {
        if !vec2.contains(v) {
            diff.push(v.to_owned());
        }
    }
    diff
}

/// Tests if the stream underlying the [BufReader] `reader` is gzipped or not by examining the
/// first 2 bytes for the magic header.  This function *requires*, but does not check, that
/// none of the stream has yet been consumed (i.e. that no read calls have yet been issued
/// to `reader`). It will fill the buffer to examine the first two bytes, but will not consume
/// them.
///
/// If the first 2 bytes could be succesfully read, this returns
/// [Ok]`(true)` if the file is a gzipped file
/// [Ok]`(false)` if it is not a gzipped file
///
/// If the first 2 bytes could not be succesfully read, then this
/// returns the relevant [std::io::Error].
///
/// Notes: implementation taken from
/// <https://github.com/zaeleus/noodles/blob/ba1b34ce22e72c2df277b20ce4c5c7b75d75a199/noodles-util/src/variant/reader/builder.rs#L131>
pub fn is_gzipped<T: BufRead>(reader: &mut T) -> std::io::Result<bool> {
    const GZIP_MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];

    let src = reader.fill_buf()?;
    if src.get(..2) == Some(&GZIP_MAGIC_NUMBER) {
        Ok(true)
    } else {
        Ok(false)
    }
}

/// Creates a [FastaReader] from the provided path. This function will automatically
/// determine if the provided path points to a gzip compressed or an uncompressed FASTA
/// file, and will return the appropriate reader accordingly.
///
/// It returns [Ok]`(`[FastaReader]`)` on success and an [anyhow::Error] on failure.
pub fn get_noodles_reader_from_path<T: AsRef<Path>>(p: T) -> anyhow::Result<FastaReader> {
    let file = std::fs::File::open(p.as_ref())?;
    let mut inner_rdr = std::io::BufReader::new(file);
    // Now, we read the fasta file and process each reference sequence at a time
    if is_gzipped(&mut inner_rdr)? {
        trace!("auto-detected gzipped FASTA file - reading via decompression");
        Ok(NoodlesFastaReader::new(Box::new(BufReader::new(
            MultiGzDecoder::new(inner_rdr),
        ))))
    } else {
        Ok(NoodlesFastaReader::new(Box::new(inner_rdr)))
    }
}

/// Creates a [FastaReader] from the provided reader. This function will automatically
/// determine if the provided reader is reading from a gzip compressed or an uncompressed FASTA
/// file, and will return the appropriate reader accordingly.
///
/// It returns [Ok]`(`[FastaReader]`)` on success and an [anyhow::Error] on failure.
///
/// **Note** : It is intended that this function *take ownership* of the underlying reader, which
/// is the reason behind the `'static` lifetime bound.
pub fn get_noodles_reader_from_reader(r: impl Read + 'static) -> anyhow::Result<FastaReader> {
    let mut inner_rdr = std::io::BufReader::new(r);
    // Now, we read the fasta file and process each reference sequence at a time
    if is_gzipped(&mut inner_rdr)? {
        trace!("auto-detected gzipped FASTA file - reading via decompression");
        Ok(NoodlesFastaReader::new(Box::new(BufReader::new(
            MultiGzDecoder::new(inner_rdr),
        ))))
    } else {
        Ok(NoodlesFastaReader::new(Box::new(inner_rdr)))
    }
}

#[derive(Clone, Copy)]
/// Represents types of intervals used in various genomic file formats.
///
/// This enum allows distinguishing between different conventions for handling interval start and end
/// points, which can vary between inclusive and exclusive, and can impact genomic data interpretation.
///
/// # Variants
///
/// * `Inclusive(i64)`: Represents intervals that include both start and end points (e.g., GTF, GFF formats).
/// * `Exclusive(i64)`: Represents intervals that exclude both start and end points (e.g., VCF format).
/// * `LeftInclusive(i64)`: Represents intervals that include the start point but exclude the end point.
/// * `RightInclusive(i64)`: Represents intervals that exclude the start point but include the end point (e.g., BED format).
///
/// The `i64` value represents an offset that might be applied to the interval (commonly 0 or 1 depending on the format).
///
/// # Default
///
/// Implements the `Default` trait, providing a default value:
///
/// * `Inclusive(1)`: Used as a safe default for most genomic formats that include both start and end points.
///
/// # Methods
///
/// * `from`: Creates an `IntervalType` based on the given file format, handling common formats like GTF, GFF, BED, SAM, and BAM.
/// * `start_offset`: Calculates the offset to be applied to the start coordinate based on the interval type.
/// * `end_offset`: Calculates the offset to be applied to the end coordinate based on the interval type.
///
/// # Examples
///
/// Determining interval type based on file format:
///
/// ```rust
/// let interval_type = IntervalType::from("BED");
/// assert_eq!(interval_type, IntervalType::RightInclusive(0));
/// ```
///
/// Calculating coordinate adjustments for interval types:
///
/// ```rust
/// let interval_type = IntervalType::Inclusive(1);
/// assert_eq!(interval_type.start_offset(), 0); // No adjustment needed for inclusive start
/// assert_eq!(interval_type.end_offset(), 0);   // No adjustment needed for inclusive end
/// ```
///
/// This enum is particularly useful when converting between genomic coordinate systems used in different file formats
/// or when implementing algorithms that require precise handling of sequence intervals.
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
    /// Provides a default interval type.
    ///
    /// This method returns a default value for the [IntervalType] enum. The default is set to
    /// [`IntervalType::Inclusive(1)`], which is a common setting for many genomic formats like GTF, GFF, and SAM,
    /// where intervals are typically inclusive, meaning they include both the start and end positions.
    ///
    /// # Returns
    ///
    /// Returns [`IntervalType::Inclusive(1)`], representing an inclusive interval with an offset of 1.
    /// This offset reflects the 1-based indexing common to certain genomic data formats.
    ///
    /// # Examples
    ///
    /// Getting a default interval type:
    ///
    /// ```rust
    /// let default_interval = IntervalType::default();
    /// assert_eq!(default_interval, IntervalType::Inclusive(1));
    /// ```
    ///
    /// This default method simplifies the initialization of interval types for common scenarios in genomic data processing,
    /// ensuring consistency and reducing the need for repetitive specification of the same interval type.
    fn default() -> Self {
        IntervalType::Inclusive(1)
    }
}

impl IntervalType {
    /// Creates an interval type based on a genomic file format.
    ///
    /// This method converts a string representing a genomic file format into an `IntervalType` instance,
    /// applying conventional interval settings for that format. This allows users to handle genomic data
    /// correctly according to the conventions used by different file formats.
    ///
    /// # Note
    ///
    /// Currently, only "GTF", "GFF", "SAM", "BAM", and "BED" formats are supported. Other formats will cause the method
    /// to panic. It is recommended to handle this potential panic in your code or ensure the format is supported
    /// before calling this method.
    ///
    /// # Type Parameters
    ///
    /// * `T`: A type that implements the [ToString] trait, which allows passing string literals,
    ///   [String] objects, or any other type that can be converted to a string representation of the file format.
    ///
    /// # Arguments
    ///
    /// * `file_type`: The file format for which to create the corresponding [IntervalType].
    ///
    /// # Returns
    ///
    /// Returns an instance of [IntervalType] corresponding to the given file format:
    /// * `IntervalType::Inclusive(1)`: For "GTF", "GFF", or "SAM" formats, which typically use inclusive intervals.
    /// * `IntervalType::RightInclusive(0)`: For "BAM" or "BED" formats, which use right-inclusive (0-based) intervals.
    ///
    /// # Panics
    ///
    /// Panics if an unsupported file format is provided. It's important to ensure that the file format
    /// is one of the supported types to prevent runtime errors.
    ///
    /// # Examples
    ///
    /// ```rust
    /// let interval_type_gtf = IntervalType::from("GTF");
    /// assert_eq!(interval_type_gtf, IntervalType::Inclusive(1));
    ///
    /// let interval_type_bed = IntervalType::from("BED");
    /// assert_eq!(interval_type_bed, IntervalType::RightInclusive(0));
    /// ```
    pub fn from<T: ToString>(file_type: T) -> Self {
        match file_type.to_string().to_lowercase().as_str() {
            "gtf" | "gff" | "sam" => IntervalType::Inclusive(1),
            "bam" | "bed" => IntervalType::RightInclusive(0),
            _ => panic!("The file type is not supported"),
        }
    }

    /// Calculates the start position offset based on the interval type.
    ///
    /// This method determines how much to adjust the start coordinate of a genomic interval
    /// to convert between different interval notations used in various file formats. The adjustment
    /// depends on the specific type of interval (inclusive, exclusive, etc.) and its inherent offset.
    ///
    /// # Returns
    ///
    /// Returns an `i64` value representing the offset to be applied to the start coordinate:
    /// * For `IntervalType::Inclusive` and `IntervalType::LeftInclusive`, returns `1 - c`, adjusting for 1-based inclusive coordinates.
    /// * For `IntervalType::RightInclusive` and `IntervalType::Exclusive`, returns `2 - c`, adjusting for exclusive start points or 0-based coordinates.
    ///
    /// Here, `c` is the custom offset associated with the interval type, typically reflecting
    /// differences in how different genomic file formats treat interval start points.
    ///
    /// # Examples
    ///
    /// Calculating start offsets for different interval types:
    ///
    /// ```rust
    /// assert_eq!(IntervalType::Inclusive(1).start_offset(), 0); // No adjustment needed for 1-based inclusive start
    /// assert_eq!(IntervalType::Exclusive(1).start_offset(), 2); // Adjusting for exclusive start, typically 0-based
    /// assert_eq!(IntervalType::RightInclusive(0).start_offset(), 2); // Same adjustment as exclusive but for right inclusive
    /// assert_eq!(IntervalType::LeftInclusive(1).start_offset(), 0); // Similar to inclusive, maintains 1-based start
    /// ```
    ///
    /// The `start_offset` method aids in the proper conversion and interpretation of genomic intervals,
    /// ensuring compatibility and correctness across different data formats and analysis workflows.
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

    /// Calculates the end position offset based on the interval type.
    ///
    /// This method determines the necessary adjustment to the end coordinate of a genomic interval
    /// to accommodate differences in interval notations (inclusive vs exclusive) across various genomic data formats.
    /// The offset is calculated based on whether the end position is considered part of the interval (inclusive)
    /// or not (exclusive).
    ///
    /// # Returns
    ///
    /// Returns an `i64` value representing the offset to be applied to the end coordinate:
    /// * For `IntervalType::Inclusive` and `IntervalType::RightInclusive`, returns `1 - c`, maintaining the end position for inclusive coordinates.
    /// * For `IntervalType::LeftInclusive` and `IntervalType::Exclusive`, returns `0 - c`, adjusting for exclusive end points or converting from 1-based to 0-based end coordinates.
    ///
    /// Here, `c` is the custom offset associated with the interval type, typically used to adapt
    /// to the differences in indexing conventions between genomic file formats.
    ///
    /// # Examples
    ///
    /// Calculating end offsets for different interval types:
    ///
    /// ```rust
    /// assert_eq!(IntervalType::Inclusive(1).end_offset(), 0); // No adjustment needed for 1-based inclusive end
    /// assert_eq!(IntervalType::Exclusive(1).end_offset(), 0); // Adjusting for exclusive end, typically reflecting removal of last position
    /// assert_eq!(IntervalType::RightInclusive(1).end_offset(), 0); // No adjustment, as end is included (typical for BED format)
    /// assert_eq!(IntervalType::LeftInclusive(1).end_offset(), 0); // Adjusting as end is exclusive, moving from 1-based to 0-based system
    /// ```
    ///
    /// The `end_offset` method is crucial for correctly handling genomic intervals during conversion
    /// and processing tasks, ensuring accurate representation and analysis of genomic data across different formats.
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

pub static FIELDCOLUMNS: [&str; 12] = [
    "seqname",
    "source",
    "feature_type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "gene_id",
    "gene_name",
    "transcript_id",
    "exon_number",
];

// --- Grangers struct related utility functionality

/*
pub struct TestIter<'a, 'b> {
    grangers: &'a Grangers,
    seq_coll: &'b GrangersSequenceCollection
}

impl<'a, 'b> Iterator for TestIter<'a, 'b> {
    type Item = (polars::frame::DataFrame::)
}

impl Grangers {
    pub fn iter_with_sequences<'a, 'b>(&'a self, seq_collection: &'b GrangersSequenceCollection) -> TestIter<'a, 'b>{
        TestIter {
            grangers: &self,
            seq_coll: seq_collection
        }
    }
}
*/
