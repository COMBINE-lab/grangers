use polars::prelude as polarsprelude;
use std::collections::HashMap;
use crate::grangers::reader::{gtf, fasta};
use std::path::Path;

/// The Grangers struct contains the following fields:
/// - df: the underlying polars dataframe
/// - comments: the comments in the GTF file
/// - chromsize: the chromosome size (set as none before calling `add_chromsize`)
/// - directives: the directives in the GFF file (set as none for GTF files)
#[derive(Clone)]
pub struct Grangers {
    file_type: gtf::FileType,
    df: polarsprelude::DataFrame,
    comments: Vec<String>,
    chromsize: Option<HashMap<String, usize>>,
    directives: Option<Vec<String>>,
}

impl Grangers {
    pub fn new(df: polarsprelude::DataFrame, comments: Vec<String>, directives: Option<Vec<String>>, file_type: gtf::FileType) -> Grangers {
        Grangers {
            file_type,
            df,
            comments,
            chromsize: None,
            directives,
        }
    }

    // Build the Grangers struct from a GTF file.\
    // Attributes except gene_id, gene_name and transcript_id 
    // will be ignored if only_essential is set to true.\
    // For a human GRh38 GTF file, the extra attributes 
    // can take more than 2GB of memory.
    pub fn from_gtf(file_path: &std::path::Path, only_essential: bool) -> anyhow::Result<Grangers> {
        let am = gtf::AttributeMode::from(!only_essential);
        let gstruct = gtf::GStruct::from_gtf(file_path, am)?;
        let gr = gstruct.to_grangers()?;
        Ok(gr)
    }

    // Build the Grangers struct from a GFF file.\
    // Attributes except ID, gene_id, gene_name and transcript_id 
    // will be ignored if only_essential is set to true.\
    // For a human GRh38 GTF file, the extra attributes 
    // can take more than 2GB of memory.
    pub fn from_gff(file_path: &std::path::Path, only_essential: bool) -> anyhow::Result<Grangers> {
        let am = gtf::AttributeMode::from(!only_essential);
        let gstruct = gtf::GStruct::from_gff(file_path, am)?;
        let gr = gstruct.to_grangers()?;
        Ok(gr)
    }

    // return the underlying polars dataframe
    pub fn df(&self) -> &polarsprelude::DataFrame {
        &self.df
    }

    // add chromsize to the Grangers struct according to a fasta file
    pub fn add_chromsize<T: AsRef<Path>>(&mut self, genome_file: T) -> anyhow::Result<()> {
        self.chromsize = Some(fasta::get_chromsize(genome_file)?);
        Ok(())
    }

    // return the chromsize
    pub fn chromsize(&self) -> Option<&HashMap<String, usize>> {
        self.chromsize.as_ref()
    }

    // return the comments
    pub fn comments(&self) -> &Vec<String> {
        &self.comments
    }

    // return the directives
    pub fn directives(&self) -> Option<&Vec<String>> {
        self.directives.as_ref()
    }
}

// implement GenomicFeatures for Grangers
impl Grangers {
    /// get the intronic sequences of each gene by the "gene_id" attribute.
    pub fn intron_by_gene(&self) -> anyhow::Result<Grangers> {

        Ok(self.clone())
    }

    pub fn intron_by_transcript(&self) -> anyhow::Result<Grangers> {

        Ok(self.clone())
    }
    pub fn merge(&mut self) -> anyhow::Result<()> {
        Ok(())
    }
}
