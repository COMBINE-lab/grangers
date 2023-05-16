use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub const GTFESSENTIALATTRIBUTES: [&str; 3] = ["gene_id", "gene_name", "transcript_id"];
pub const GFFESSENTIALATTRIBUTES: [&str; 4] = ["ID","gene_id", "gene_name", "transcript_id"];

/// traverse the given file to get the number of lines in the file
pub fn _file_line_count<T: AsRef<Path>>(file_path: T) -> anyhow::Result<usize> {
    let reader = BufReader::new(File::open(file_path)?);
    let mut num_lines = 0usize;
    for l in reader.lines() {
        let line = l?;
        if !(line.trim().starts_with("#") | line.trim().is_empty()) {
            num_lines += 1;
        }
    }
    Ok(num_lines)
}
