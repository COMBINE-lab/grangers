use crate::grangers_utils::{self, equal_length};
use anyhow::bail;
use noodles::fasta;
use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

/// Represents sequence information for genomic or reference datasets.
///
/// This struct holds detailed information about sequences, including their names, lengths,
/// circularity, associated genome, and any additional custom data. It's particularly useful
/// for handling metadata of genomic datasets.
///
/// # Fields
///
/// * `seqname`: A vector of [String] containing the names of the sequences (e.g., chromosomes).
/// * `seqlengths`: An optional vector of `usize` representing the lengths of each sequence.
/// * `is_circular`: An optional vector of `bool` indicating whether each sequence is circular.
/// * `genome`: An optional [String] representing the name or identifier of the genome.
/// * `extra`: An optional [HashMap] for storing additional metadata as key-value pairs.
///
/// # Methods
///
/// ## Constructor
///
/// * `new`: Constructs a new [SeqInfo] instance from provided sequence details.
///
/// ## Accessors
///
/// * `seqname`: Returns a reference to the vector of sequence names.
/// * `seqlengths`: Returns a reference to the optional vector of sequence lengths.
/// * `is_circular`: Returns a reference to the optional vector indicating circularity.
/// * `genome`: Returns a reference to the optional genome identifier.
/// * `extra`: Returns a reference to the optional additional metadata.
///
/// ## Mutators
///
/// * `set_seqnames`: Sets the sequence names for the [SeqInfo] instance.
///
/// ## From Fasta
///
/// * `from_fasta`: Constructs a new [SeqInfo] instance by parsing a FASTA file to extract
///   sequence names and lengths.
///
/// # Examples
///
/// Creating a new [SeqInfo] instance:
///
/// ```rust
/// let seq_info = SeqInfo::new(
///     vec!["chr1".to_string(), "chr2".to_string()],
///     Some(vec![248956422, 242193529]),
///     Some(vec![false, false]),
///     Some("Human".to_string()),
///     None,
/// )?;
/// ```
///
/// # Note
///
/// Ensure that the vectors for `seqname`, `seqlengths`, and `is_circular` (if provided) are of the same length
/// to maintain consistency across sequence metadata. Mismatches in length will result in errors during instantiation.
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

    /// Creates a new instance of [SeqInfo].
    ///
    /// This function initializes a [SeqInfo] struct with provided sequence information,
    /// validating that the lengths of sequence names, sequence lengths, and circularity indicators match,
    /// ensuring data integrity and consistency.
    ///
    /// # Arguments
    ///
    /// * `seqname`: A vector of [String] representing the names of the sequences (e.g., chromosome names).
    /// * `seqlengths`: An optional vector of `usize` representing the lengths of each sequence.
    /// * `is_circular`: An optional vector of `bool` indicating whether each sequence is circular.
    /// * `genome`: An optional [String] representing the name or identifier of the genome.
    /// * `extra`: An optional [HashMap] for storing additional metadata as key-value pairs.
    ///
    /// # Returns
    ///
    /// Returns [`anyhow::Result<SeqInfo>`](anyhow::Result):
    /// * [Ok]`(`[SeqInfo]`)`: If the instance is successfully created.
    /// * [Err]`(`[anyhow::Error]`)`: If there is a mismatch in the lengths of provided vectors or other issues.
    ///
    /// # Examples
    ///
    /// ```rust
    /// let seq_info = SeqInfo::new(
    ///     vec!["chr1".to_string(), "chr2".to_string()],
    ///     Some(vec![248956422, 242193529]),
    ///     Some(vec![false, false]),
    ///     Some("Human".to_string()),
    ///     None,
    /// )?;
    /// ```
    ///
    /// # Errors
    ///
    /// An error is returned if:
    /// * The lengths of `seqname`, `seqlengths`, `is_circular`, or any key in `extra` do not match.
    /// * There are other issues in creating the [SeqInfo] instance.
    pub fn new(
        seqname: Vec<String>,
        seqlengths: Option<Vec<usize>>,
        is_circular: Option<Vec<bool>>,
        genome: Option<String>,
        extra: Option<HashMap<String, Vec<String>>>,
    ) -> anyhow::Result<SeqInfo> {
        if let Some(v) = &seqlengths {
            if equal_length(&seqname[..], &v[..]) {
                bail!("seqname and seqlengths have different length; Could not create SeqInfo")
            }
        }
        if let Some(v) = &is_circular {
            if equal_length(&seqname[..], &v[..]) {
                bail!("seqname and is_circular have different length; Could not create SeqInfo")
            }
        }
        if let Some(hm) = &extra {
            for (k, v) in hm.iter() {
                if equal_length(&seqname[..], &v[..]) {
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

    /// ## From Fasta
    ///
    /// * `from_fasta`: Constructs a new [SeqInfo] instance by parsing a FASTA file to extract
    ///   sequence names and lengths.
    ///
    /// seqinfo can contain "seqname", "seqlengths", "isCircular", "genome", and other extra string info
    ///
    /// # Examples
    ///
    /// Creating a new [SeqInfo] instance:
    ///
    /// ```rust
    /// let seq_info = SeqInfo::new(
    ///     vec!["chr1".to_string(), "chr2".to_string()],
    ///     Some(vec![248956422, 242193529]),
    ///     Some(vec![false, false]),
    ///     Some("Human".to_string()),
    ///     None,
    /// )?;
    /// ```
    ///
    /// # Note
    ///
    /// Ensure that the vectors for `seqname`, `seqlengths`, and `is_circular` (if provided) are of the same length
    /// to maintain consistency across sequence metadata. Mismatches in length will result in errors during instantiation.
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

/// Extracts sequence names and lengths from a FASTA file.
///
/// This function parses a FASTA file to extract the names (identifiers) and lengths of the sequences
/// it contains, typically representing chromosomal or contig information in genomic datasets.
///
/// # Type Parameters
///
/// * `T`: A type that can be referenced as a file path, implementing the `AsRef<Path>` trait.
///
/// # Arguments
///
/// * `file_path`: The path to the FASTA file.
///
/// # Returns
///
/// Returns `anyhow::Result<(Vec<String>, Vec<usize>)>`:
/// * `Ok((Vec<String>, Vec<usize>))`: A tuple containing two vectors, the first with sequence names and
///   the second with their corresponding lengths.
/// * `Err(anyhow::Error)`: An error if there is a problem opening the file or reading from it.
///
/// # Examples
///
/// ```rust
/// let (seqnames, seqlengths) = get_chromsize("path/to/sequence.fasta")?;
/// ```
///
/// This function is useful for generating metadata about the sequences in a FASTA file, such as before
/// creating an index or performing sequence length-dependent calculations.
pub fn get_chromsize<T: AsRef<Path>>(file_path: T) -> anyhow::Result<(Vec<String>, Vec<usize>)> {
    // create reader
    let mut reader = grangers_utils::get_noodles_reader_from_path(file_path)?;

    // call inner function
    _get_chromsize(&mut reader)
}

fn _get_chromsize<T: BufRead>(
    rdr: &mut fasta::io::Reader<T>,
) -> anyhow::Result<(Vec<String>, Vec<usize>)> {
    // get chromsize
    let mut seqname: Vec<String> = Vec::new();
    let mut seqlengths: Vec<usize> = Vec::new();

    for result in rdr.records() {
        let record = result?;

        let record_name = std::str::from_utf8(record.name())?;

        seqname.push(
            record_name
                .split_once(' ')
                .unwrap_or((record_name, ""))
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
    use std::io::Cursor;
    #[test]
    fn test_get_chromsize() {
        let fasta_data = b">sq0 test\nACGT\n>sq1\nNNNN\nNNNN\nNN\n";
        let reader = Cursor::new(fasta_data);
        let mut rdr =
            grangers_utils::get_noodles_reader_from_reader(reader).expect("can't construct reader");

        let chromsize = _get_chromsize(&mut rdr).unwrap();
        assert_eq!(chromsize.1.first().unwrap(), &4);
    }
}
