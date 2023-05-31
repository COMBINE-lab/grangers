use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Copy, Clone)]
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
    pub fn get_essential(&self) -> &[&str] {
        match self {
            FileFormat::GTF => GTFESSENTIALATTRIBUTES.as_ref(),
            FileFormat::GFF => GFFESSENTIALATTRIBUTES.as_ref(),
            _ => &[],
        }
    }
    pub fn is_gtf(&self) -> bool {
        matches!(self, FileFormat::GTF)
    }
}
impl std::str::FromStr for FileFormat {
    type Err = anyhow::Error;

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

pub const GTFESSENTIALATTRIBUTES: [&str; 3] = ["gene_id", "gene_name", "transcript_id"];
pub const GFFESSENTIALATTRIBUTES: [&str; 4] = ["ID", "gene_id", "gene_name", "transcript_id"];
// pub const FIELDS: [&str; 8] = [
//     "seqid",
//     "source",
//     "feature_type",
//     "start",
//     "end",
//     "score",
//     "strand",
//     "phase",
// ];
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

pub fn equal_length<T, R>(vec1: &Vec<T>, vec2: &Vec<R>) -> bool {
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

/// Tests if the stream underlying the BufReader `reader` is gzipped or not by examining the
/// first 2 bytes for the magic header.  This function *requires*, but does not check, that
/// none of the stream has yet been consumed (i.e. that no read calls have yet been issued
/// to `reader`). It will fill the buffer to examine the first two bytes, but will not consume
/// them.
///
/// If the first 2 bytes could be succesfully read, this returns
/// Ok(true) if the file is a gzipped file
/// Ok(false) if it is not a gzipped file
///
/// If the first 2 bytes could not be succesfully read, then this
/// returns the relevant std::io::Error.
///
/// Notes: implementation taken from
/// https://github.com/zaeleus/noodles/blob/ba1b34ce22e72c2df277b20ce4c5c7b75d75a199/noodles-util/src/variant/reader/builder.rs#L131
pub fn is_gzipped<T: BufRead>(reader: &mut T) -> std::io::Result<bool> {
    const GZIP_MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];

    let src = reader.fill_buf()?;
    if src.get(..2) == Some(&GZIP_MAGIC_NUMBER) {
        Ok(true)
    } else {
        Ok(false)
    }
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
