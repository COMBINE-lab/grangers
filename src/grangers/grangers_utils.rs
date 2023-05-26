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
