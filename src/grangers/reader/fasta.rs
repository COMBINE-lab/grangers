use crate::grangers::grangers_utils::equal_length;
use anyhow::bail;
use noodles::fasta;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
#[derive(Clone)]
pub struct SeqInfo {
    seqname: Vec<String>,
    seqlengths: Option<Vec<usize>>,
    is_circular: Option<Vec<bool>>,
    genome: Option<String>,
    extra: Option<HashMap<String, Vec<String>>>,
}

impl SeqInfo {
    /// get the seqname of the genome/reference set
    pub fn seqname(&self) -> &Vec<String> {
        &self.seqname
    }

    /// get the seqlengths of the genome/reference set
    pub fn seqlengths(&self) -> &Option<Vec<usize>> {
        &self.seqlengths
    }

    /// get the is_circular of the genome/reference set
    pub fn is_circular(&self) -> &Option<Vec<bool>> {
        &self.is_circular
    }

    /// get the genome of the genome/reference set
    pub fn genome(&self) -> &Option<String> {
        &self.genome
    }

    /// get the extra of the genome/reference set
    pub fn extra(&self) -> &Option<HashMap<String, Vec<String>>> {
        &self.extra
    }

    /// set the seqname of the genome/reference set
    pub fn set_seqnames(&mut self, seqname: Vec<String>) {
        self.seqname = seqname;
    }

    pub fn new(
        seqname: Vec<String>,
        seqlengths: Option<Vec<usize>>,
        is_circular: Option<Vec<bool>>,
        genome: Option<String>,
        extra: Option<HashMap<String, Vec<String>>>,
    ) -> anyhow::Result<SeqInfo> {
        if let Some(v) = &seqlengths {
            if equal_length(&seqname, v) {
                bail!("seqname and seqlengths have different length; Could not create SeqInfo")
            }
        }
        if let Some(v) = &is_circular {
            if equal_length(&seqname, v) {
                bail!("seqname and is_circular have different length; Could not create SeqInfo")
            }
        }
        if let Some(hm) = &extra {
            for (k, v) in hm.iter() {
                if equal_length(&seqname, v) {
                    bail!(
                        "seqname and {} have different length; Could not create SeqInfo",
                        k
                    )
                }
            }
        }

        Ok(SeqInfo {
            seqname,
            seqlengths,
            is_circular,
            genome,
            extra: None,
        })
    }

    /// This function parse a fasta file to get the seqinfo of the genome/reference set.
    /// seqinfo can contain "seqname", "seqlengths", "isCircular", "genome", and other extra string info
    pub fn from_fasta<T: AsRef<Path>>(file_path: T) -> anyhow::Result<SeqInfo> {
        // get chromsize
        let (seqname, seqlengths) = get_chromsize(&file_path)?;
        let si = SeqInfo::new(
            seqname,
            Some(seqlengths),
            None,
            Some(file_path.as_ref().to_string_lossy().to_string()),
            None,
        );
        si
    }
}
/// traverses the given fasta file returns a tuple of seqname and seqlengths
pub fn build_fasta_reader<T: AsRef<Path>, R: BufRead>(
    file_path: T,
) -> anyhow::Result<fasta::Reader<BufReader<File>>> {
    // create reader
    let reader: fasta::Reader<BufReader<File>> = File::open(file_path)
        .map(BufReader::new)
        .map(fasta::Reader::new)?;
    Ok(reader)
}

/// traverses the given fasta file returns a tuple of seqname and seqlengths
pub fn get_chromsize<T: AsRef<Path>>(file_path: T) -> anyhow::Result<(Vec<String>, Vec<usize>)> {
    // create reader
    let mut reader = File::open(file_path)
        .map(BufReader::new)
        .map(fasta::Reader::new)?;

    // call inner function
    _get_chromsize(&mut reader)
}

fn _get_chromsize<T: BufRead>(
    rdr: &mut fasta::Reader<T>,
) -> anyhow::Result<(Vec<String>, Vec<usize>)> {
    // get chromsize
    let mut seqname: Vec<String> = Vec::new();
    let mut seqlengths: Vec<usize> = Vec::new();

    for result in rdr.records() {
        let record = result?;
        seqname.push(
            record
                .name()
                .split_once(' ')
                .unwrap_or((record.name(), ""))
                .0
                .to_string(),
        );
        seqlengths.push(record.sequence().len());
    }

    Ok((seqname, seqlengths))
}

#[cfg(test)]
mod tests {

    use super::*;
    #[test]
    fn test_get_chromsize() {
        let fasta_data = b">sq0 test\nACGT\n>sq1\nNNNN\nNNNN\nNN\n";
        let mut rdr = fasta::Reader::new(&fasta_data[..]);

        let chromsize = _get_chromsize(&mut rdr).unwrap();
        assert_eq!(chromsize.1.first().unwrap(), &4);
    }
}
