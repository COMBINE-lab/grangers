use crate::grangers_utils::*;
use anyhow::{self, Ok};
use noodles::gff::record::{Phase, Strand};
use noodles::{gff, gtf};
use polars::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;
use std::{
    collections::HashMap,
    hash::Hash,
    path::{self, Path, PathBuf},
};
use tracing::{info, warn};

pub enum FileType {
    GTF,
    GFF2,
    GFF3,
}

pub enum GStructMode {
    Essential,
    Full,
}

impl GStructMode {
    pub fn from(is_full: bool) -> GStructMode {
        if is_full {
            GStructMode::Full
        } else {
            GStructMode::Essential
        }
    }

    pub fn is_full(&self) -> bool {
        match self {
            GStructMode::Full => true,
            GStructMode::Essential => false,
        }
    }
}

/// This is a wrapper of noodles Attribute struct.
/// As noodles define the Attribute struct twice, for gff and gtf,
/// this is a uniform wrapper

// pub struct Attributes;
pub enum Attributes {
    GeneID,
    GeneName,
    TranscriptID,
}

impl AsRef<str> for Attributes {
    fn as_ref(&self) -> &str {
        match self {
            Self::GeneID => "gene_id",
            Self::GeneName => "gene_name",
            Self::TranscriptID => "transcript_id",
        }
    }
}

pub enum FeatureType {
    Gene,
    Transcript,
    Exon,
    Other
}

impl std::str::FromStr for FeatureType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> anyhow::Result<FeatureType> {
        let ft = match s {
            "gene" => FeatureType::Gene,
            "transcript" => FeatureType::Transcript,
            "exon" => FeatureType::Exon,
            _ => FeatureType::Other,
        };
        Ok(ft)
    }
}


/// This struct contains all information in a GTF file. it will be used to construct the
/// polars data frame. If this is no faster than generating
pub struct GStruct {
    pub mode: GStructMode,
    pub seqid: Vec<String>,
    pub source: Vec<String>,
    pub feature_type: Vec<String>,
    pub start: Vec<u64>,
    pub end: Vec<u64>,
    pub score: Vec<Option<f32>>,
    pub strand: Vec<Option<String>>,
    pub phase: Vec<Option<String>>,
    pub gene_id: Vec<Option<String>>,
    pub gene_name: Vec<Option<String>>,
    pub transcript_id: Vec<Option<String>>,
    // this is memory intensive: use 1G more memory
    pub attributes: Option<HashMap<String, Vec<Option<String>>>>,
}

impl GStruct {
    pub fn new(mode: GStructMode) -> anyhow::Result<GStruct> {
        let gr = GStruct {
            seqid: Vec::new(),
            source: Vec::new(),
            feature_type: Vec::new(),
            start: Vec::new(),
            end: Vec::new(),
            score: Vec::new(),
            strand: Vec::new(),
            phase: Vec::new(),
            gene_id: Vec::new(),
            gene_name: Vec::new(),
            transcript_id: Vec::new(),
            attributes: if mode.is_full() {Some(HashMap::new())} else {None},
            mode,
        };
        Ok(gr)
    }
    pub fn from_gtf<T: AsRef<Path>>(gtf_file: T, gm: GStructMode) -> anyhow::Result<GStruct> {
        

        // instantiate the struct with a proper size
        let mut rdr: gtf::Reader<BufReader<File>> = gtf::Reader::new(BufReader::new(File::open(gtf_file)?));
        let mut gr = GStruct::new(gm)?;
        // initiate a hashmap for attributes. This hashmap will receive all attributes in all records
        let mut attribute_hms: HashMap<String, Vec<Option<String>>> =  HashMap::new();
        // initiate a hashmap for the attributes in one record
        let mut rec_attr_hm: HashMap<String, String> = HashMap::with_capacity(100);

        // we want to record how many records have invalid fields or attributes.
        let mut invalid_strand = 0usize;
        let mut missing_gene_id = 0usize;
        let mut missing_gene_name = 0usize;  
        let mut missing_transcript_id = 0usize;

        // let mut s1 = Series::new("1", vec![Attributes::GeneID.as_ref(), Attributes::GeneID.as_ref(),Attributes::GeneName.as_ref()])
        // .cast(&DataType::Categorical(None))
        // .unwrap();
        let start = Instant::now();
        let mut recid = 0usize;
        // parse the file
        for l in rdr.lines() {
            let line = l?;
    
            // ignore comments
            if let gtf::Line::Record(r) = line {
                // receive fields
                // noodles reader will take care most of thing for us. 
                // however, we do need to impose something we required.
                // That is, we require each none-gene-type record to have a transcript_id 
                // build attribute hashmap

                // then, parse essential fields
                let seqid = r.reference_sequence_name().to_string();
                let source = r.source().to_string();
                let feature_type: String = r.ty().to_string();
                let start = r.start().get().to_owned() as u64;
                let end: u64 = r.end().get().to_owned() as u64;
                let score = r.score();
                let strand = if let Some(st) = r.strand() {
                    Some(st.as_ref().to_owned())
                } else {
                    None
                };

                // gtf frame is the same thing as gff phase
                let phase = if let Some(ph) = r.frame() {
                    Some(ph.to_string())
                } else {
                    None
                };
                rec_attr_hm.extend(r.attributes().iter().map(|a| (a.key().to_string(),a.value().to_string())));

                // parse attributes
                // for (k, v) in rec_attr_hm.iter() {
                //     // if the attribute doesn't exist, insert a vector for it
                //     if !gr.attributes.contains_key(k) {
                //         gr.attributes.insert(k.to_string(), vec![None;nl]);
                //     };

                //     // then, insert value for this record according to line number
                //     if let Some(av) = gr.attributes.get_mut(k) {
                //         // it is safe to use unwrap as we created the vec
                //         let mptr = av.get_mut(recid).unwrap();
                //         *mptr = Some(v.to_owned());
                //     }
                // }


                // // if tid does not exist, skip
                // let mut transcript_id = if let Some(tid) = rec_attr_hm.get(Attributes::TranscriptID.as_ref()) {
                //     Some(tid.to_owned())
                // } else if is_gene {
                //     missing_transcript_id += 1;
                //     continue
                // } else {
                //     None
                // };

                // // if strand is unknown, skip
                // let st = if let Some(st) = r.strand() {
                //     Some(st.as_ref().to_owned())
                // } else {
                //     invalid_strand+=1;
                //     None
                // };

                // // gtf frame is the same thing as gff phase
                // let fr = if let Some(fr) = r.frame() {
                //     Some(fr.to_string())
                // } else {
                //     None
                // };
    
                // parse essential attributes
                // we need gene id for all
                let transcript_id = if let Some(tid) = rec_attr_hm.get(Attributes::TranscriptID.as_ref()) {
                    Some(tid.to_owned())
                } else {
                    None
                };

                let gene_id = if let Some(gid) = rec_attr_hm.get(Attributes::GeneID.as_ref()) {
                    Some(gid.to_owned())
                } else {
                    None
                };

                let gene_name = if let Some(gn) = rec_attr_hm.get(Attributes::GeneName.as_ref()) {
                    Some(gn.to_owned())
                } else {
                    None
                };
                // let gene_id: String;
                // let gene_name: String;
                // if ogene_id.is_none() && ogene_name.is_none() {
                //     missing_gene_id += 1;
                //     missing_gene_name += 1;
                //     // if it is a gene 
                //     if is_gene {
                //         continue
                //     } else {
                //         gene_id = transcript_id.to_owned();
                //         gene_name = transcript_id.to_owned();
                //     }
                // } else {
                //     gene_id = if let Some(gid) = ogene_id {
                //         gid.to_owned()
                //     } else {
                //         missing_gene_id += 1;
                //         ogene_name.unwrap().to_owned()
                //     };
                //     gene_name = if let Some(gname) = ogene_name {
                //         gname.to_owned()
                //     } else {
                //         missing_gene_name += 1;
                //         ogene_id.unwrap().to_owned()
                //     };
                // }
                gr.seqid.push(seqid);
                gr.source.push(source);
                gr.feature_type.push(feature_type);
                gr.start.push(start);
                gr.end.push(end);
                gr.score.push(score);
                gr.strand.push(strand);
                gr.phase.push(phase);
                gr.gene_id.push(gene_id);
                gr.gene_name.push(gene_name);
                gr.transcript_id.push(transcript_id);
                rec_attr_hm.clear();
                recid += 1;
            }
        }

        let duration = start.elapsed();
        println!("Runtime: {:?}", duration);

        Ok(gr)
    }

    pub fn to_df(self) -> anyhow::Result<DataFrame> {
        // create dataframe!
        // for fields
        let mut df_vec = vec![
            Series::new("seqid", self.seqid),
            Series::new("source", self.source),
            Series::new("type", self.feature_type),
            Series::new("start", self.start),
            Series::new("end", self.end),
            Series::new("score", self.score),
            Series::new("strand", self.strand),
            Series::new("phase", self.phase),
        ];
        
        if let Some(attributes) = self.attributes {
            for (k,v) in attributes {
                df_vec.push(
                    Series::new(k.as_str(), v),
                )
            }
        }
        let gr = DataFrame::new(df_vec)?;
        Ok(gr)
    }
}
