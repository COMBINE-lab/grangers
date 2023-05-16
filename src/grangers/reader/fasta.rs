use noodles::fasta;
use std::path::Path;
use std::fs::File;
use std::io::{BufReader, BufRead};
use std::collections::HashMap;

pub fn get_chromsize<T: AsRef<Path>>(file_path: T) -> anyhow::Result<HashMap<String, usize>> {
    // create reader
    let mut reader = File::open(file_path)
    .map(BufReader::new)
    .map(fasta::Reader::new)?;

    // get chromsize
    let chromsize = _get_chromsize(&mut reader)?;
    Ok(chromsize)
}

pub fn _get_chromsize<T: BufRead>(rdr: &mut fasta::Reader<T>) -> anyhow::Result<HashMap<String, usize>> {

    // get chromsize
    let mut chromsize = HashMap::new();

    for result in rdr.records() {
        let record = result?;
        // if chromosome name contains multiple parts separated by space, only use the first part
        chromsize.insert(record.name().split_once(' ').unwrap_or((record.name(),"")).0.to_string(), record.sequence().len());
    }

    Ok(chromsize)
}


#[cfg(test)]
mod tests {

    use super::*;
    #[test]
    fn test_get_chromsize() {
        let fasta_data = b">sq0 test\nACGT\n>sq1\nNNNN\nNNNN\nNN\n";
        let mut rdr = fasta::Reader::new(&fasta_data[..]);
    
        let chromsize = _get_chromsize(&mut rdr).unwrap();
        assert_eq!(chromsize.get("sq0").unwrap(), &4);
    }


}