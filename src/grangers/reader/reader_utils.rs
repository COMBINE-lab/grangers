use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub const GTFESSENTIALATTRIBUTES: [&str; 3] = ["gene_id", "gene_name", "transcript_id"];
pub const GFFESSENTIALATTRIBUTES: [&str; 4] = ["ID","gene_id", "gene_name", "transcript_id"];
pub const FIELDS: [&str;8] = ["seqid", "source", "feature_type", "start", "end", "score", "strand", "phase"];
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


pub fn equal_length<T,R>(vec1: &Vec<T>, vec2: &Vec<R>) -> bool {
    vec1.len() == vec2.len()
}

pub fn setdiff<T: Eq+Clone>(vec1: &[T],vec2: &[T]) -> Vec<T> {
    let mut diff: Vec<T> = Vec::new();
    for v in vec1.iter() {
        if !vec2.contains(v) {
            diff.push(v.to_owned());
        }
    }
    diff
}