use anyhow::{self, Ok};
use noodles::gff::record::{Phase, Strand};
use noodles::{gff, gtf};
use std::fs::File;
use std::{
    collections::HashMap,
    hash::Hash,
    path::{self, Path, PathBuf},
};

use std::time::Instant;
use tracing::{info, warn};

use reader::*;
use std::io::{BufRead, BufReader};
pub mod grangers_utils;
pub mod reader;

// This function traverses the GTF file:
// - If everything is good, count the # of transcripts for each gene and number of exons for each txp to instantiate the struct
// - If there is any invalid record, report an error and write a clean_gtf.gtf file.
//     We care only about records with type gene, transcript and exon.
//     1. Each transcript/exon record must have a valid `transcript_id` attribute.
//         - If this is not satisfied, it returns an error.
//         - Only the records with a valid transcript_id will be written to the clean_gtf.gtf.
//     2. For gene_id and gene_name attributes,
//         - If these two attributes do not exist in any record, An error will be returned.
//             At the same time, in the clean_gtf.gtf, the two fields will be imputed using the transcript_id fields.
//             However, this means each transcript will have a count because there is no transcript to gene mapping.
//         - If one of these two fields is completely missing, a warning will be generated, and the missing field will be imputed using the other.
//         - if some records miss gene_id and/or gene_name, a warning will be printed,
//             and the missing values will be imputed by the following rules:
//                 For records miss gene_id or gene_name, impute the missing one using the other one;
//                 If both are missing, impute them using transcript_id, which cannot be missing.
//    3. If there is no "transcript" or "gene" feature record, a warning will be printed.
//         Moreover, those missing records will be imputed using the "exon" feature records:
//             The Start and End site of the gene/transcript will be imputed as the bounds of their corresponding exons.
//    4. If the boundaries defined in the transcripts'/genes' feature records do not match those implied by their exons' feature records,
//         report a warning but still use transcripts'/genes' feature records to extract unspliced sequences.
//         To be specific, if some but not all transcripts/genes have their corresponding transcripts'/genes' feature records,
//         or the Start and/or End site defined in the transcript/gene feature records do not match the corresponding exons' bounds,
//         then the existing transcripts'/genes' feature records will be used to extract unspliced transcripts.
//         At the same time, in the clean_gtf.gtf, all genes/transcripts that appeared in the exon feature records
//         will have their corresponding transcripts'/genes' feature records, in which the boundaries match the corresponding exons' bounds.

/// Go through the provided file
pub fn parse<T: AsRef<Path>>(
    file_path: T,
    file_type: &FileType,
    temp_dir: T,
) -> anyhow::Result<()> {
    let mut gene_hm: HashMap<String, u32> = HashMap::new();
    let mut tx_hm: HashMap<String, u32> = HashMap::new();

    // initiate temp.gtf writer
    //

    // build GTF reader
    let f = File::open(file_path)?;
    let f = BufReader::new(f);
    let mut reader = gtf::Reader::new(f);
    let mut num_g_rec: usize = 0;
    let mut num_tx_rec: usize = 0;
    let mut num_g_rec_no_gid: usize = 0;
    let mut num_g_rec_no_gidname: usize = 0;

    //  read lines
    for (i, r) in reader.lines().enumerate() {
        let rec = r?;
        // we ignore comments and parse records
        if let gtf::line::Line::Record(r) = rec {
            match r.ty() {
                "gene" => {
                    num_g_rec += 1;
                    // get gid
                    let gid =
                        if let (Some(gid)) = r.attributes().iter().find(|x| x.key() == "gene_id") {
                            gid
                        } else if let (Some(gname)) =
                            r.attributes().iter().find(|x| x.key() == "gene_name")
                        {
                            num_g_rec_no_gid += 1;
                            gname
                        } else {
                            num_g_rec_no_gidname += 1;
                            continue;
                        };
                    gene_hm.entry(gid.value().to_string()).or_insert(1);
                }
                "transcript" => {}
                "exon" => {}
                _ => {}
            }
            println!("{:?}", r);
        }
        if i == 6 {
            break;
        }
    }

    // a gene record must have one of gene_name and gene_id

    // a transcript/exon record must have transcript_id

    // we also count the number of transcript for each gene, and # of exon for each tx

    Ok(())
}

struct GTF {}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ptr::read_unaligned;
    const GTF_RECORD: &[u8] = b"##provider: GENCODE\nchr1\tHAVANA\tgene\t29554\t31109\t.\t+\t.\tgene_id \"ENSG00000243485\"; gene_version \"5\"; gene_type \"lncRNA\"; gene_name \"MIR1302-2HG\"; level 2; hgnc_id \"HGNC:52482\"; tag \"ncRNA_host\"; havana_gene \"OTTHUMG00000000959.2\";\nchr1\tHAVANA\ttranscript\t29554\t31097\t.\t+\t.\tgene_id \"ENSG00000243485\"; gene_version \"5\"; transcript_id \"ENST00000473358\"; transcript_version \"1\"; gene_type \"lncRNA\"; gene_name \"MIR1302-2HG\"; transcript_type \"lncRNA\"; transcript_name \"MIR1302-2HG-202\"; level 2; transcript_support_level \"5\"; hgnc_id \"HGNC:52482\"; tag \"not_best_in_genome_evidence\"; tag \"dotter_confirmed\"; tag \"basic\"; havana_gene \"OTTHUMG00000000959.2\"; havana_transcript \"OTTHUMT00000002840.1\";\nchr1\tHAVANA\texon\t29554\t30039\t.\t+\t.\tgene_id \"ENSG00000243485\"; gene_version \"5\"; transcript_id \"ENST00000473358\"; transcript_version \"1\"; gene_type \"lncRNA\"; gene_name \"MIR1302-2HG\"; transcript_type \"lncRNA\"; transcript_name \"MIR1302-2HG-202\"; exon_number 1; exon_id \"ENSE00001947070\"; exon_version \"1\"; level 2; transcript_support_level \"5\"; hgnc_id \"HGNC:52482\"; tag \"not_best_in_genome_evidence\"; tag \"dotter_confirmed\"; tag \"basic\"; havana_gene \"OTTHUMG00000000959.2\"; havana_transcript \"OTTHUMT00000002840.1\";\nchr1\tHAVANA\texon\t30564\t30667\t.\t+\t.\tgene_id \"ENSG00000243485\"; gene_version \"5\"; transcript_id \"ENST00000473358\"; transcript_version \"1\"; gene_type \"lncRNA\"; gene_name \"MIR1302-2HG\"; transcript_type \"lncRNA\"; transcript_name \"MIR1302-2HG-202\"; exon_number 2; exon_id \"ENSE00001922571\"; exon_version \"1\"; level 2; transcript_support_level \"5\"; hgnc_id \"HGNC:52482\"; tag \"not_best_in_genome_evidence\"; tag \"dotter_confirmed\"; tag \"basic\"; havana_gene \"OTTHUMG00000000959.2\"; havana_transcript \"OTTHUMT00000002840.1\";\nchr1\tHAVANA\ttranscript\t30267\t31109\t.\t+\t.\tgene_id \"ENSG00000243485\"; gene_version \"5\"; transcript_id \"ENST00000469289\"; transcript_version \"1\"; gene_type \"lncRNA\"; gene_name \"MIR1302-2HG\"; transcript_type \"lncRNA\"; transcript_name \"MIR1302-2HG-201\"; level 2; transcript_support_level \"5\"; hgnc_id \"HGNC:52482\"; tag \"not_best_in_genome_evidence\"; tag \"basic\"; havana_gene \"OTTHUMG00000000959.2\"; havana_transcript \"OTTHUMT00000002841.2\";";

    #[test]
    fn test_validate() {
        let mut rdr = gtf::Reader::new(GTF_RECORD);
        let mut lines = rdr.lines();
        println!("{:?}", lines.next().transpose());
        println!("{:?}", lines.next().transpose());
        println!("{:?}", lines.next().transpose());
        println!("{:?}", lines.next().transpose());
        println!("{:?}", lines.next().transpose());
        println!("{:?}", lines.next().transpose());
        assert!(true)
    }
}
