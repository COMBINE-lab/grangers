use anyhow;
use noodles::{gff, gtf};
use polars::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::{collections::HashMap, path::Path};
use tracing::info;
use super::reader_utils::*;
use crate::grangers::grangers::Grangers;


#[derive(Copy, Clone)]
pub enum FileType {
    GTF,
    GFF,
}

impl FileType {
    pub fn get_essential(&self) -> &[&str] {
        match self {
            FileType::GTF => GTFESSENTIALATTRIBUTES.as_ref(),
            FileType::GFF => GFFESSENTIALATTRIBUTES.as_ref(),
        }
    }
}

#[derive(Copy, Clone)]
pub enum AttributeMode {
    Essential,
    Full,
}

impl AttributeMode {
    pub fn from(is_full: bool) -> AttributeMode {
        if is_full {
            AttributeMode::Full
        } else {
            AttributeMode::Essential
        }
    }

    pub fn is_full(&self) -> bool {
        match self {
            AttributeMode::Full => true,
            AttributeMode::Essential => false,
        }
    }
}

#[derive(Clone)]
struct Attributes {
    file_type: FileType,
    essential: HashMap<String, Vec<Option<String>>>,
    extra: Option<HashMap<String, Vec<Option<String>>>>,
    tally: usize,
}

impl Attributes {
    fn new(mode: AttributeMode, file_type: FileType) -> anyhow::Result<Attributes> {
        // create essential from an iterator
        let essential = HashMap::from_iter(
            file_type.get_essential()
                .iter()
                .map(|s| (s.to_string(), Vec::new())),
        );

        // if in full mode, create extra
        let extra = if mode.is_full() {
            Some(HashMap::new())
        } else {
            None
        };
        Ok(Attributes { file_type, essential, extra, tally: 0})
    }

    fn push(&mut self, hm: &mut HashMap<String, String>) {
        // parse essential attributes
        for &ea in self.file_type.get_essential() {
            if let Some(vec) = self.essential.get_mut(ea) {
                vec.push(hm.remove(ea))
            };
        }
        // the rest items are all extra attributes
        // parse them if we are in full modes
        if let Some(extra) = &mut self.extra {

            // append existing attributes
            extra.iter_mut().for_each(|(k, v)| {
                v.push(hm.remove(k));
            });

            // parse the rest attributes
            if !hm.is_empty() {
                // if there is any attribute left, create a new vector for it
                // the length should be the same as the essential attributes before insertion
                for (attr_name, attr_value) in hm {
                    extra.insert(attr_name.to_string(), {
                        let mut vec = vec![None; self.tally];
                        vec.push(Some(attr_value.to_string()));
                        vec
                    });
                }
            }
        }

        self.tally += 1;
    }
}

/// This is a wrapper of noodles Attribute struct.
/// As noodles define the Attribute struct twice, for gff and gtf,
/// this is a uniform wrapper

// pub struct Attributes;

#[derive(Copy, Clone)]
pub enum FeatureType {
    Gene,
    Transcript,
    Exon,
    Other,
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
#[derive(Clone)]
pub struct GStruct {
    seqid: Vec<String>,
    source: Vec<String>,
    feature_type: Vec<String>,
    start: Vec<u64>,
    end: Vec<u64>,
    score: Vec<Option<f32>>,
    strand: Vec<Option<String>>,
    phase: Vec<Option<String>>,
    attributes: Attributes,
    comments: Vec<String>,
    directives: Option<Vec<String>>,
}

// implement GTF reader
impl GStruct {
    fn _from_gtf<T: BufRead>(&mut self, rdr: &mut gtf::Reader<T>) -> anyhow::Result<()> {
        // initiate a reusable hashmap to take the attributes of each record
        let mut rec_attr_hm: HashMap<String, String> = HashMap::with_capacity(100);
        let mut n_comments = 0usize;
        let mut n_records = 0usize;

        // parse the file
        for l in rdr.lines() {
            let line = l?;
            match line {
                gtf::Line::Record(r) => {
                    n_records += 1;
                    // parse essential fields

                    GStruct::push(&mut self.seqid, r.reference_sequence_name().to_string());
                    GStruct::push(&mut self.source, r.source().to_string());
                    GStruct::push(&mut self.feature_type, r.ty().to_string());
                    GStruct::push(&mut self.start, r.start().get().to_owned() as u64);
                    GStruct::push(&mut self.end, r.end().get().to_owned() as u64);
                    GStruct::push(&mut self.score, r.score());
                    GStruct::push(&mut self.strand, if let Some(st) = r.strand() {
                        Some(st.as_ref().to_owned())
                    } else {
                        None
                    });

                    GStruct::push(&mut self.phase, if let Some(ph) = r.frame() {
                        Some(ph.to_string())
                    } else {
                        None
                    });

                    // parse attributes
                    rec_attr_hm.clear();
                    for attr in r.attributes().iter() {
                        rec_attr_hm.insert(attr.key().to_string(), attr.value().to_string());
                    }
                    self.attributes.push(&mut rec_attr_hm);
                },
                gtf::Line::Comment(c) => {
                    n_comments += 1;
                    self.comments.push(c);
                    continue;
                },
            }
        }
        info!(
            "Finished parsing the input file. Found {} comments and {} records.",
            n_comments, n_records
        );
        Ok(())
    }

    pub fn from_gtf<T: AsRef<Path>>(file_path: T, am: AttributeMode) -> anyhow::Result<GStruct> {
        // instantiate the struct
        let mut rdr: gtf::Reader<BufReader<File>> =
            gtf::Reader::new(BufReader::new(File::open(file_path)?));

        let mut gr = GStruct::new(am, FileType::GTF)?;
        gr._from_gtf(&mut rdr)?;
        Ok(gr)
    }
}

// implement GFF reader
impl GStruct {
    pub fn from_gff<T: AsRef<Path>>(file_path: T, am: AttributeMode) -> anyhow::Result<GStruct> {
        // instantiate the struct
        let mut rdr =
            gff::Reader::new(BufReader::new(File::open(file_path)?));

        let mut gr = GStruct::new(am, FileType::GFF)?;
        gr._from_gff(&mut rdr)?;
        Ok(gr)
    }
    fn _from_gff<T: BufRead>(&mut self, rdr: &mut gff::Reader<T>) -> anyhow::Result<()> {
        // initiate a reusable hashmap to take the attributes of each record
        let mut rec_attr_hm: HashMap<String, String> = HashMap::with_capacity(100);
        let mut n_comments = 0usize;
        let mut n_records = 0usize;
        let mut n_directives = 0usize;

        // parse the file
        for l in rdr.lines() {
            let line = l?;
            match line {
                gff::Line::Record(r) => {
                    n_records += 1;
                    GStruct::push(&mut self.seqid, r.reference_sequence_name().to_string());
                    GStruct::push(&mut self.source, r.source().to_string());
                    GStruct::push(&mut self.feature_type, r.ty().to_string());
                    GStruct::push(&mut self.start, r.start().get().to_owned() as u64);
                    GStruct::push(&mut self.end, r.end().get().to_owned() as u64);
                    GStruct::push(&mut self.score, r.score());
                    GStruct::push(&mut self.strand, match r.strand() {
                        gff::record::Strand::Forward | gff::record::Strand::Reverse => {
                            Some(r.strand().to_string())
                        }
                        _ => None,
                    });
                    GStruct::push(&mut self.phase, if let Some(ph) = r.phase() {
                        Some(ph.to_string())
                    } else {
                        None
                    });

                    // parse attributes
                    rec_attr_hm.clear();
                    for attr in r.attributes().iter() {
                        rec_attr_hm.insert(attr.key().to_string(), attr.value().to_string());
                    }
                    self.attributes.push(&mut rec_attr_hm);
                }
                gff::Line::Comment(c) => {
                    n_comments += 1;
                    self.comments.push(c);
                    continue;
                }
                gff::Line::Directive(d) => {
                    n_directives += 1;
                    // this must be Some
                    if let Some(directives) = &mut self.directives {
                        directives.push(d.to_string());
                    }
                    continue;
                }
            }
        }
        info!(
            "Finished parsing the input file. Found {} comments, {} directives, and {} records.",
            n_comments, n_directives, n_records
        );
        Ok(())
    }
}

// implenment general functions
impl GStruct {
    pub fn new(attribute_mode: AttributeMode, file_type: FileType) -> anyhow::Result<GStruct> {
        let gr = GStruct {
            seqid: Vec::new(),
            source: Vec::new(),
            feature_type: Vec::new(),
            start: Vec::new(),
            end: Vec::new(),
            score: Vec::new(),
            strand: Vec::new(),
            phase: Vec::new(),
            attributes: Attributes::new(attribute_mode,file_type)?,
            comments: Vec::new(),
            directives: match file_type {
                FileType::GTF => None,
                FileType::GFF => Some(Vec::new()),
            },
        };
        Ok(gr)
    }

    // TODO: might need a better generic type
    fn push<T: std::fmt::Debug+Clone>(vec: &mut Vec<T>, val: T) {
        vec.push(val);
    }

    pub fn to_grangers(self) -> anyhow::Result<Grangers> {
        // create dataframe!
        // we want to make some columns categorical because of this https://docs.rs/polars/latest/polars/docs/performance/index.html
        // fields
        let mut df_vec = vec![
            Series::new("seqid", self.seqid)
                .cast(&DataType::Categorical(None))
                .unwrap(),
            Series::new("source", self.source)
                .cast(&DataType::Categorical(None))
                .unwrap(),
            Series::new("type", self.feature_type)
                .cast(&DataType::Categorical(None))
                .unwrap(),
            Series::new("start", self.start),
            Series::new("end", self.end),
            Series::new("score", self.score),
            Series::new("strand", self.strand)
                .cast(&DataType::Categorical(None))
                .unwrap(),
            Series::new("phase", self.phase)
                .cast(&DataType::Categorical(None))
                .unwrap(),
        ];

        //for essential attributes
        // categorical as we want to do operations on them
        for (k, v) in self.attributes.essential {
            let s = Series::new(k.as_str(), v)
                .cast(&DataType::Categorical(None))
                .unwrap();
            df_vec.push(s);
        }

        // for extra attributes
        if let Some(attributes) = self.attributes.extra {
            for (k, v) in attributes {
                df_vec.push(Series::new(k.as_str(), v))
            }
        }
        let df = DataFrame::new(df_vec)?;

        let gr = Grangers ::new(df, self.comments, self.directives, self.attributes.file_type);
        Ok(gr)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const GTF_RECORD: &[u8] = b"##provider: GENCODE\nchr1\tHAVANA\tgene\t29554\t31109\t.\t+\t.\tgene_id \"ENSG00000243485\"; gene_version \"5\"; gene_type \"lncRNA\"; gene_name \"MIR1302-2HG\"; level 2; hgnc_id \"HGNC:52482\"; tag \"ncRNA_host\"; havana_gene \"OTTHUMG00000000959.2\";\nchr1\tHAVANA\ttranscript\t29554\t31097\t.\t+\t.\tgene_id \"ENSG00000243485\"; gene_version \"5\"; transcript_id \"ENST00000473358\"; transcript_version \"1\"; gene_type \"lncRNA\"; gene_name \"MIR1302-2HG\"; transcript_type \"lncRNA\"; transcript_name \"MIR1302-2HG-202\"; level 2; transcript_support_level \"5\"; hgnc_id \"HGNC:52482\"; tag \"not_best_in_genome_evidence\"; tag \"dotter_confirmed\"; tag \"basic\"; havana_gene \"OTTHUMG00000000959.2\"; havana_transcript \"OTTHUMT00000002840.1\";\nchr1\tHAVANA\texon\t29554\t30039\t.\t+\t.\tgene_id \"ENSG00000243485\"; gene_version \"5\"; transcript_id \"ENST00000473358\"; transcript_version \"1\"; gene_type \"lncRNA\"; gene_name \"MIR1302-2HG\"; transcript_type \"lncRNA\"; transcript_name \"MIR1302-2HG-202\"; exon_number 1; exon_id \"ENSE00001947070\"; exon_version \"1\"; level 2; transcript_support_level \"5\"; hgnc_id \"HGNC:52482\"; tag \"not_best_in_genome_evidence\"; tag \"dotter_confirmed\"; tag \"basic\"; havana_gene \"OTTHUMG00000000959.2\"; havana_transcript \"OTTHUMT00000002840.1\";\nchr1\tHAVANA\texon\t30564\t30667\t.\t+\t.\tgene_id \"ENSG00000243485\"; gene_version \"5\"; transcript_id \"ENST00000473358\"; transcript_version \"1\"; gene_type \"lncRNA\"; gene_name \"MIR1302-2HG\"; transcript_type \"lncRNA\"; transcript_name \"MIR1302-2HG-202\"; exon_number 2; exon_id \"ENSE00001922571\"; exon_version \"1\"; level 2; transcript_support_level \"5\"; hgnc_id \"HGNC:52482\"; tag \"not_best_in_genome_evidence\"; tag \"dotter_confirmed\"; tag \"basic\"; havana_gene \"OTTHUMG00000000959.2\"; havana_transcript \"OTTHUMT00000002840.1\";\nchr1\tHAVANA\ttranscript\t30267\t31109\t.\t+\t.\tgene_id \"ENSG00000243485\"; gene_version \"5\"; transcript_id \"ENST00000469289\"; transcript_version \"1\"; gene_type \"lncRNA\"; gene_name \"MIR1302-2HG\"; transcript_type \"lncRNA\"; transcript_name \"MIR1302-2HG-201\"; level 2; transcript_support_level \"5\"; hgnc_id \"HGNC:52482\"; tag \"not_best_in_genome_evidence\"; tag \"basic\"; havana_gene \"OTTHUMG00000000959.2\"; havana_transcript \"OTTHUMT00000002841.2\";";

    const GFF_RECORD: &[u8] = 
    b"##gff-version 3\n#description: evidence-based annotation of the human genome (GRCh38), version 43 (Ensembl 109)\n#provider: GENCODE\n#contact: gencode-help@ebi.ac.uk\n#format: gff3\n#date: 2022-11-29\n##sequence-region chr1 1 248956422\nchr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tID=ENSG00000290825.1;gene_id=ENSG00000290825.1;gene_type=lncRNA;gene_name=DDX11L2;level=2;tag=overlaps_pseudogene\nchr1\tHAVANA\ttranscript\t11869\t14409\t.\t+\t.\tID=ENST00000456328.2;Parent=ENSG00000290825.1;gene_id=ENSG00000290825.1;transcript_id=ENST00000456328.2;gene_type=lncRNA;gene_name=DDX11L2;transcript_type=lncRNA;transcript_name=DDX11L2-202;level=2;transcript_support_level=1;tag=basic,Ensembl_canonical;havana_transcript=OTTHUMT00000362751.1\nchr1\tHAVANA\texon\t11869\t12227\t.\t+\t.\tID=exon:ENST00000456328.2:1;Parent=ENST00000456328.2;gene_id=ENSG00000290825.1;transcript_id=ENST00000456328.2;gene_type=lncRNA;gene_name=DDX11L2;transcript_type=lncRNA;transcript_name=DDX11L2-202;exon_number=1;exon_id=ENSE00002234944.1;level=2;transcript_support_level=1;tag=basic,Ensembl_canonical;havana_transcript=OTTHUMT00000362751.1\nchr1\tHAVANA\texon\t12613\t12721\t.\t+\t.\tID=exon:ENST00000456328.2:2;Parent=ENST00000456328.2;gene_id=ENSG00000290825.1;transcript_id=ENST00000456328.2;gene_type=lncRNA;gene_name=DDX11L2;transcript_type=lncRNA;transcript_name=DDX11L2-202;exon_number=2;exon_id=ENSE00003582793.1;level=2;transcript_support_level=1;tag=basic,Ensembl_canonical;havana_transcript=OTTHUMT00000362751.1\nchr1\tHAVANA\texon\t13221\t14409\t.\t+\t.\tID=exon:ENST00000456328.2:3;Parent=ENST00000456328.2;gene_id=ENSG00000290825.1;transcript_id=ENST00000456328.2;gene_type=lncRNA;gene_name=DDX11L2;transcript_type=lncRNA;transcript_name=DDX11L2-202;exon_number=3;exon_id=ENSE00002312635.1;level=2;transcript_support_level=1;tag=basic,Ensembl_canonical;havana_transcript=OTTHUMT00000362751.1\n";

    #[test]
    fn test_from_gtf() {
        let mut rdr = gtf::Reader::new(GTF_RECORD);
        let mut gr = GStruct::new(AttributeMode::Full,FileType::GTF).unwrap();
        gr._from_gtf(&mut rdr).unwrap();
        // check values
        match gr {
            GStruct{seqid, source, feature_type, start, end, score, strand, phase, attributes, comments, directives} => {
                assert_eq!(seqid, vec![String::from("chr1");5]);
                assert_eq!(source, vec![String::from("HAVANA");5]);
                assert_eq!(feature_type, vec![String::from("gene"),String::from("transcript"),String::from("exon"),String::from("exon"),String::from("transcript")]);
                assert_eq!(start, vec![29554,29554,29554,30564,30267]);
                assert_eq!(end, vec![31109,31097,30039,30667,31109]);
                assert_eq!(score, vec![None;5]);
                assert_eq!(strand, vec![Some(String::from("+"));5]);
                assert_eq!(phase, vec![None;5]);
                match attributes {
                    Attributes{file_type,essential,extra, tally} => {
                        assert!(file_type.get_essential() == &["gene_id", "gene_name", "transcript_id"]);
                        assert!(essential.get("gene_id").unwrap().iter().map(|v| v.clone().unwrap().eq(&String::from("ENSG00000243485"))).collect::<Vec<bool>>().iter().all(|v| *v));

                        assert!(essential.get("gene_name").unwrap().iter().map(|v| v.clone().unwrap().eq(&String::from("MIR1302-2HG"))).collect::<Vec<bool>>().iter().all(|v| *v));

                        assert_eq!(essential.get("transcript_id").unwrap().iter().map(|v| if let Some(id) = v.clone() {id} else {String::from("none")}).collect::<Vec<String>>(), vec![String::from("none"), String::from("ENST00000473358"), String::from("ENST00000473358"), String::from("ENST00000473358"), String::from("ENST00000469289")]);
                        assert_eq!(extra.unwrap().get("gene_type").unwrap().iter().map(|v| if let Some(id) = v.clone() {id} else {String::from("none")}).collect::<Vec<String>>(), vec![String::from("lncRNA");5]);
                        assert_eq!(tally, 5);
                    },
                }
            },            
        }
    }

    #[test]
    fn test_from_gff() {
        let mut rdr = gff::Reader::new(GFF_RECORD);
        let mut gr = GStruct::new(AttributeMode::Full, FileType::GFF).unwrap();
        gr._from_gff(&mut rdr).unwrap();
        // check values
        match gr {
            GStruct{seqid, source, feature_type, start, end, score, strand, phase, attributes, comments, directives} => {
                assert_eq!(seqid, vec![String::from("chr1");5]);
                assert_eq!(source, vec![String::from("HAVANA");5]);
                assert_eq!(feature_type, vec![String::from("gene"),String::from("transcript"),String::from("exon"),String::from("exon"),String::from("exon")]);
                assert_eq!(start, vec![11869,11869,11869,12613,13221]);
                assert_eq!(end, vec![14409,14409,12227,12721,14409]);
                assert_eq!(score, vec![None;5]);
                assert_eq!(strand, vec![Some(String::from("+"));5]);
                assert_eq!(phase, vec![None;5]);
                match attributes {
                    Attributes{file_type,essential,extra,tally} => {
                        assert!(file_type.get_essential() == &["ID", "gene_id", "gene_name", "transcript_id"]);
                        assert!(essential.get("gene_id").unwrap().iter().map(|v| v.clone().unwrap().eq(&String::from("ENSG00000290825.1"))).collect::<Vec<bool>>().iter().all(|v| *v));

                        assert!(essential.get("gene_name").unwrap().iter().map(|v| v.clone().unwrap().eq(&String::from("DDX11L2"))).collect::<Vec<bool>>().iter().all(|v| *v));

                        assert_eq!(essential.get("transcript_id").unwrap().iter().map(|v| if let Some(id) = v.clone() {id} else {String::from("none")}).collect::<Vec<String>>(), vec![String::from("none"), String::from("ENST00000456328.2"), String::from("ENST00000456328.2"), String::from("ENST00000456328.2"), String::from("ENST00000456328.2")]);
                        assert_eq!(extra.unwrap().get("gene_type").unwrap().iter().map(|v| if let Some(id) = v.clone() {id} else {String::from("none")}).collect::<Vec<String>>(), vec![String::from("lncRNA");5]);
                        assert_eq!(tally, 5);
                    },
                }
            },            
        }
    }
}
