use super::reader_utils::equal_length;
use anyhow::bail;
use noodles::fasta;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
#[derive(Clone)]
pub struct SeqInfo {
    seqnames: Vec<String>,
    seqlengths: Option<Vec<usize>>,
    is_circular: Option<Vec<bool>>,
    genome: Option<String>,
    extra: Option<HashMap<String, Vec<String>>>,
}

impl SeqInfo {
    pub fn new(
        seqnames: Vec<String>,
        seqlengths: Option<Vec<usize>>,
        is_circular: Option<Vec<bool>>,
        genome: Option<String>,
        extra: Option<HashMap<String, Vec<String>>>,
    ) -> anyhow::Result<SeqInfo> {
        if let Some(v) = &seqlengths {
            if equal_length(&seqnames, v) {
                bail!("seqnames and seqlengths have different length; Could not create SeqInfo")
            }
        }
        if let Some(v) = &is_circular {
            if equal_length(&seqnames, &v) {
                bail!("seqnames and is_circular have different length; Could not create SeqInfo")
            }
        }
        if let Some(hm) = &extra {
            for (k, v) in hm.iter() {
                if equal_length(&seqnames, &v) {
                    bail!(
                        "seqnames and {} have different length; Could not create SeqInfo",
                        k
                    )
                }
            }
        }

        Ok(SeqInfo {
            seqnames,
            seqlengths,
            is_circular,
            genome,
            extra: None,
        })
    }

    /// This function parse a fasta file to get the seqinfo of the genome/reference set.
    /// seqinfo can contain "seqnames", "seqlengths", "isCircular", "genome", and other extra string info
    pub fn from_fasta<T: AsRef<Path>>(file_path: T) -> anyhow::Result<SeqInfo> {
        // get chromsize
        let (seqnames, seqlengths) = get_chromsize(&file_path)?;
        let si = SeqInfo::new(
            seqnames,
            Some(seqlengths),
            None,
            Some(file_path.as_ref().to_string_lossy().to_string()),
            None,
        );
        si
    }
}

/// traverses the given fasta file returns a tuple of seqnames and seqlengths
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
    let mut seqnames: Vec<String> = Vec::new();
    let mut seqlengths: Vec<usize> = Vec::new();

    for result in rdr.records() {
        let record = result?;
        seqnames.push(
            record
                .name()
                .split_once(' ')
                .unwrap_or((record.name(), ""))
                .0
                .to_string(),
        );
        seqlengths.push(record.sequence().len());
    }

    Ok((seqnames, seqlengths))
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
