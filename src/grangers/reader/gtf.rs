use crate::grangers::grangers_utils::{is_gzipped, FileFormat};
use anyhow;
use flate2::bufread::GzDecoder;
use noodles::{gff, gtf};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::{collections::HashMap, path::Path};
use tracing::info;

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
pub struct Attributes {
    pub file_type: FileFormat,
    pub essential: HashMap<String, Vec<Option<String>>>,
    pub extra: Option<HashMap<String, Vec<Option<String>>>>,
    pub tally: usize,
}

impl Attributes {
    pub fn new(mode: AttributeMode, file_type: FileFormat) -> anyhow::Result<Attributes> {
        // create essential from an iterator
        let essential = HashMap::from_iter(
            file_type
                .get_essential()
                .iter()
                .map(|s| (s.to_string(), Vec::with_capacity(1_0000))),
        );

        // if in full mode, create extra
        let extra = if mode.is_full() {
            Some(HashMap::with_capacity(100))
        } else {
            None
        };
        Ok(Attributes {
            file_type,
            essential,
            extra,
            tally: 0,
        })
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
    pub seqid: Vec<String>,
    pub source: Vec<String>,
    pub feature_type: Vec<String>,
    pub start: Vec<i64>,
    pub end: Vec<i64>,
    pub score: Vec<Option<f32>>,
    pub strand: Vec<Option<String>>,
    pub phase: Vec<Option<String>>,
    pub attributes: Attributes,
    pub misc: Option<HashMap<String, Vec<String>>>,
}

// implement GTF reader
impl GStruct {
    pub fn from_gtf<T: AsRef<Path>>(file_path: T, am: AttributeMode) -> anyhow::Result<GStruct> {
        let mut gr = GStruct::new(am, FileFormat::GTF)?;
        if let Some(misc) = gr.misc.as_mut() {
            misc.insert(String::from("file_type"), vec![String::from("GTF")]);
        }

        let file = File::open(file_path)?;
        let mut inner_rdr = BufReader::new(file);
        // instantiate the struct
        if is_gzipped(&mut inner_rdr)? {
            info!("auto-detected gzipped file - reading via decompression");
            let mut rdr = gtf::Reader::new(BufReader::new(GzDecoder::new(inner_rdr)));
            gr._from_gtf(&mut rdr)?;
        } else {
            let mut rdr = gtf::Reader::new(inner_rdr);
            gr._from_gtf(&mut rdr)?;
        }

        Ok(gr)
    }

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
                    GStruct::push(&mut self.start, r.start().get().to_owned() as i64);
                    GStruct::push(&mut self.end, r.end().get().to_owned() as i64);
                    GStruct::push(&mut self.score, r.score());
                    GStruct::push(
                        &mut self.strand,
                        r.strand().map(|st| st.as_ref().to_owned()),
                    );

                    GStruct::push(&mut self.phase, r.frame().map(|ph| ph.to_string()));

                    // parse attributes
                    rec_attr_hm.clear();
                    for attr in r.attributes().iter() {
                        rec_attr_hm.insert(attr.key().to_string(), attr.value().to_string());
                    }
                    self.attributes.push(&mut rec_attr_hm);
                }
                gtf::Line::Comment(c) => {
                    n_comments += 1;
                    if let Some(misc) = self.misc.as_mut() {
                        misc.entry(String::from("comments"))
                            .and_modify(|v| v.push(c.clone()))
                            .or_insert(vec![c]);
                    }
                    continue;
                }
            }
        }
        info!(
            "Finished parsing the input file. Found {} comments and {} records.",
            n_comments, n_records
        );
        Ok(())
    }
}

// implement GFF reader
impl GStruct {
    pub fn from_gff<T: AsRef<Path>>(file_path: T, am: AttributeMode) -> anyhow::Result<GStruct> {
        let mut gr = GStruct::new(am, FileFormat::GFF)?;

        let file = File::open(file_path)?;
        let mut inner_rdr = BufReader::new(file);
        // instantiate the struct
        if is_gzipped(&mut inner_rdr)? {
            info!("auto-detected gzipped file - reading via decompression");
            let mut rdr = gff::Reader::new(BufReader::new(GzDecoder::new(inner_rdr)));
            gr._from_gff(&mut rdr)?;
        } else {
            let mut rdr = gff::Reader::new(inner_rdr);
            gr._from_gff(&mut rdr)?;
        }
        Ok(gr)
    }

    fn _from_gff<T: BufRead>(&mut self, rdr: &mut gff::Reader<T>) -> anyhow::Result<()> {
        // initiate a reusable hashmap to take the attributes of each record
        let mut rec_attr_hm: HashMap<String, String> = HashMap::with_capacity(100);
        let mut n_comments = 0usize;
        let mut n_records = 0usize;
        // parse the file
        for l in rdr.lines() {
            let line = l?;
            match line {
                gff::Line::Record(r) => {
                    n_records += 1;
                    GStruct::push(&mut self.seqid, r.reference_sequence_name().to_string());
                    GStruct::push(&mut self.source, r.source().to_string());
                    GStruct::push(&mut self.feature_type, r.ty().to_string());
                    GStruct::push(&mut self.start, r.start().get().to_owned() as i64);
                    GStruct::push(&mut self.end, r.end().get().to_owned() as i64);
                    GStruct::push(&mut self.score, r.score());
                    GStruct::push(
                        &mut self.strand,
                        match r.strand() {
                            gff::record::Strand::Forward | gff::record::Strand::Reverse => {
                                Some(r.strand().to_string())
                            }
                            _ => None,
                        },
                    );
                    GStruct::push(&mut self.phase, r.phase().map(|ph| ph.to_string()));

                    // parse attributes
                    rec_attr_hm.clear();
                    for attr in r.attributes().iter() {
                        rec_attr_hm.insert(attr.key().to_string(), attr.value().to_string());
                    }
                    self.attributes.push(&mut rec_attr_hm);
                }
                gff::Line::Comment(c) => {
                    n_comments += 1;
                    if let Some(misc) = self.misc.as_mut() {
                        misc.entry(String::from("comments"))
                            .and_modify(|v| v.push(c.clone()))
                            .or_insert(vec![c]);
                    }
                    continue;
                }
                gff::Line::Directive(d) => {
                    let dstring = d.to_string();
                    // this must be Some

                    if let Some(misc) = self.misc.as_mut() {
                        misc.entry(String::from("directives"))
                            .and_modify(|v| v.push(dstring.clone()))
                            .or_insert(vec![dstring]);
                    }
                    continue;
                }
            }
        }
        info!(
            "Finished parsing the input file. Found {} comments, and {} records.",
            n_comments, n_records
        );
        Ok(())
    }
}

// implenment general functions
impl GStruct {
    pub fn new(attribute_mode: AttributeMode, file_type: FileFormat) -> anyhow::Result<GStruct> {
        let gr = GStruct {
            seqid: Vec::with_capacity(1_0000),
            source: Vec::with_capacity(1_0000),
            feature_type: Vec::with_capacity(1_0000),
            start: Vec::with_capacity(1_0000),
            end: Vec::with_capacity(1_0000),
            score: Vec::with_capacity(1_0000),
            strand: Vec::with_capacity(1_0000),
            phase: Vec::with_capacity(1_0000),
            attributes: Attributes::new(attribute_mode, file_type)?,
            misc: Some(HashMap::new()),
        };
        Ok(gr)
    }

    // TODO: might need a better generic type
    fn push<T: std::fmt::Debug + Clone>(vec: &mut Vec<T>, val: T) {
        vec.push(val);
    }

    pub fn append<T: ToString + Clone>(vec: &mut Vec<T>, patch: &mut Vec<T>) {
        vec.append(patch);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const GTF_RECORD: &[u8] = b"##provider: GENCODE\nchr1\tHAVANA\tgene\t29554\t31109\t.\t+\t.\tgene_id \"ENSG00000243485\"; gene_version \"5\"; gene_type \"lncRNA\"; gene_name \"MIR1302-2HG\"; level 2; hgnc_id \"HGNC:52482\"; tag \"ncRNA_host\"; havana_gene \"OTTHUMG00000000959.2\";\nchr1\tHAVANA\ttranscript\t29554\t31097\t.\t+\t.\tgene_id \"ENSG00000243485\"; gene_version \"5\"; transcript_id \"ENST00000473358\"; transcript_version \"1\"; gene_type \"lncRNA\"; gene_name \"MIR1302-2HG\"; transcript_type \"lncRNA\"; transcript_name \"MIR1302-2HG-202\"; level 2; transcript_support_level \"5\"; hgnc_id \"HGNC:52482\"; tag \"not_best_in_genome_evidence\"; tag \"dotter_confirmed\"; tag \"basic\"; havana_gene \"OTTHUMG00000000959.2\"; havana_transcript \"OTTHUMT00000002840.1\";\nchr1\tHAVANA\texon\t29554\t30039\t.\t+\t.\tgene_id \"ENSG00000243485\"; gene_version \"5\"; transcript_id \"ENST00000473358\"; transcript_version \"1\"; gene_type \"lncRNA\"; gene_name \"MIR1302-2HG\"; transcript_type \"lncRNA\"; transcript_name \"MIR1302-2HG-202\"; exon_number 1; exon_id \"ENSE00001947070\"; exon_version \"1\"; level 2; transcript_support_level \"5\"; hgnc_id \"HGNC:52482\"; tag \"not_best_in_genome_evidence\"; tag \"dotter_confirmed\"; tag \"basic\"; havana_gene \"OTTHUMG00000000959.2\"; havana_transcript \"OTTHUMT00000002840.1\";\nchr1\tHAVANA\texon\t30564\t30667\t.\t+\t.\tgene_id \"ENSG00000243485\"; gene_version \"5\"; transcript_id \"ENST00000473358\"; transcript_version \"1\"; gene_type \"lncRNA\"; gene_name \"MIR1302-2HG\"; transcript_type \"lncRNA\"; transcript_name \"MIR1302-2HG-202\"; exon_number 2; exon_id \"ENSE00001922571\"; exon_version \"1\"; level 2; transcript_support_level \"5\"; hgnc_id \"HGNC:52482\"; tag \"not_best_in_genome_evidence\"; tag \"dotter_confirmed\"; tag \"basic\"; havana_gene \"OTTHUMG00000000959.2\"; havana_transcript \"OTTHUMT00000002840.1\";\nchr1\tHAVANA\ttranscript\t30267\t31109\t.\t+\t.\tgene_id \"ENSG00000243485\"; gene_version \"5\"; transcript_id \"ENST00000469289\"; transcript_version \"1\"; gene_type \"lncRNA\"; gene_name \"MIR1302-2HG\"; transcript_type \"lncRNA\"; transcript_name \"MIR1302-2HG-201\"; level 2; transcript_support_level \"5\"; hgnc_id \"HGNC:52482\"; tag \"not_best_in_genome_evidence\"; tag \"basic\"; havana_gene \"OTTHUMG00000000959.2\"; havana_transcript \"OTTHUMT00000002841.2\";";

    const GFF_RECORD: &[u8] = b"##gff-version 3\n#description: evidence-based annotation of the human genome (GRCh38), version 43 (Ensembl 109)\n#provider: GENCODE\n#contact: gencode-help@ebi.ac.uk\n#format: gff3\n#date: 2022-11-29\n##sequence-region chr1 1 248956422\nchr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tID=ENSG00000290825.1;gene_id=ENSG00000290825.1;gene_type=lncRNA;gene_name=DDX11L2;level=2;tag=overlaps_pseudogene\nchr1\tHAVANA\ttranscript\t11869\t14409\t.\t+\t.\tID=ENST00000456328.2;Parent=ENSG00000290825.1;gene_id=ENSG00000290825.1;transcript_id=ENST00000456328.2;gene_type=lncRNA;gene_name=DDX11L2;transcript_type=lncRNA;transcript_name=DDX11L2-202;level=2;transcript_support_level=1;tag=basic,Ensembl_canonical;havana_transcript=OTTHUMT00000362751.1\nchr1\tHAVANA\texon\t11869\t12227\t.\t+\t.\tID=exon:ENST00000456328.2:1;Parent=ENST00000456328.2;gene_id=ENSG00000290825.1;transcript_id=ENST00000456328.2;gene_type=lncRNA;gene_name=DDX11L2;transcript_type=lncRNA;transcript_name=DDX11L2-202;exon_number=1;exon_id=ENSE00002234944.1;level=2;transcript_support_level=1;tag=basic,Ensembl_canonical;havana_transcript=OTTHUMT00000362751.1\nchr1\tHAVANA\texon\t12613\t12721\t.\t+\t.\tID=exon:ENST00000456328.2:2;Parent=ENST00000456328.2;gene_id=ENSG00000290825.1;transcript_id=ENST00000456328.2;gene_type=lncRNA;gene_name=DDX11L2;transcript_type=lncRNA;transcript_name=DDX11L2-202;exon_number=2;exon_id=ENSE00003582793.1;level=2;transcript_support_level=1;tag=basic,Ensembl_canonical;havana_transcript=OTTHUMT00000362751.1\nchr1\tHAVANA\texon\t13221\t14409\t.\t+\t.\tID=exon:ENST00000456328.2:3;Parent=ENST00000456328.2;gene_id=ENSG00000290825.1;transcript_id=ENST00000456328.2;gene_type=lncRNA;gene_name=DDX11L2;transcript_type=lncRNA;transcript_name=DDX11L2-202;exon_number=3;exon_id=ENSE00002312635.1;level=2;transcript_support_level=1;tag=basic,Ensembl_canonical;havana_transcript=OTTHUMT00000362751.1\n";

    #[test]
    fn test_from_gtf() {
        let mut rdr = gtf::Reader::new(GTF_RECORD);
        let mut gr = GStruct::new(AttributeMode::Full, FileFormat::GTF).unwrap();
        gr._from_gtf(&mut rdr).unwrap();
        // check values
        match gr {
            GStruct {
                seqid,
                source,
                feature_type,
                start,
                end,
                score,
                strand,
                phase,
                attributes,
                misc: _,
            } => {
                assert_eq!(seqid, vec![String::from("chr1"); 5]);
                assert_eq!(source, vec![String::from("HAVANA"); 5]);
                assert_eq!(
                    feature_type,
                    vec![
                        String::from("gene"),
                        String::from("transcript"),
                        String::from("exon"),
                        String::from("exon"),
                        String::from("transcript")
                    ]
                );
                assert_eq!(start, vec![29554, 29554, 29554, 30564, 30267]);
                assert_eq!(end, vec![31109, 31097, 30039, 30667, 31109]);
                assert_eq!(score, vec![None; 5]);
                assert_eq!(strand, vec![Some(String::from("+")); 5]);
                assert_eq!(phase, vec![None; 5]);
                match attributes {
                    Attributes {
                        file_type,
                        essential,
                        extra,
                        tally,
                    } => {
                        assert!(
                            file_type == FileFormat::GTF
                        );
                        assert!(essential
                            .get("gene_id")
                            .unwrap()
                            .iter()
                            .map(|v| v.clone().unwrap().eq(&String::from("ENSG00000243485")))
                            .collect::<Vec<bool>>()
                            .iter()
                            .all(|v| *v));

                        assert!(essential
                            .get("gene_name")
                            .unwrap()
                            .iter()
                            .map(|v| v.clone().unwrap().eq(&String::from("MIR1302-2HG")))
                            .collect::<Vec<bool>>()
                            .iter()
                            .all(|v| *v));

                        assert_eq!(
                            essential
                                .get("transcript_id")
                                .unwrap()
                                .iter()
                                .map(|v| if let Some(id) = v.clone() {
                                    id
                                } else {
                                    String::from("none")
                                })
                                .collect::<Vec<String>>(),
                            vec![
                                String::from("none"),
                                String::from("ENST00000473358"),
                                String::from("ENST00000473358"),
                                String::from("ENST00000473358"),
                                String::from("ENST00000469289")
                            ]
                        );
                        assert_eq!(
                            extra
                                .unwrap()
                                .get("gene_type")
                                .unwrap()
                                .iter()
                                .map(|v| if let Some(id) = v.clone() {
                                    id
                                } else {
                                    String::from("none")
                                })
                                .collect::<Vec<String>>(),
                            vec![String::from("lncRNA"); 5]
                        );
                        assert_eq!(tally, 5);
                    }
                }
            }
        }
    }

    #[test]
    fn test_from_gff() {
        let mut rdr = gff::Reader::new(GFF_RECORD);
        let mut gr = GStruct::new(AttributeMode::Full, FileFormat::GFF).unwrap();
        gr._from_gff(&mut rdr).unwrap();
        // check values
        match gr {
            GStruct {
                seqid,
                source,
                feature_type,
                start,
                end,
                score,
                strand,
                phase,
                attributes,
                misc: _,
            } => {
                assert_eq!(seqid, vec![String::from("chr1"); 5]);
                assert_eq!(source, vec![String::from("HAVANA"); 5]);
                assert_eq!(
                    feature_type,
                    vec![
                        String::from("gene"),
                        String::from("transcript"),
                        String::from("exon"),
                        String::from("exon"),
                        String::from("exon")
                    ]
                );
                assert_eq!(start, vec![11869, 11869, 11869, 12613, 13221]);
                assert_eq!(end, vec![14409, 14409, 12227, 12721, 14409]);
                assert_eq!(score, vec![None; 5]);
                assert_eq!(strand, vec![Some(String::from("+")); 5]);
                assert_eq!(phase, vec![None; 5]);
                match attributes {
                    Attributes {
                        file_type,
                        essential,
                        extra,
                        tally,
                    } => {
                        assert!(
                            file_type == FileFormat::GFF

                        );
                        assert!(essential
                            .get("gene_id")
                            .unwrap()
                            .iter()
                            .map(|v| v.clone().unwrap().eq(&String::from("ENSG00000290825.1")))
                            .collect::<Vec<bool>>()
                            .iter()
                            .all(|v| *v));

                        assert!(essential
                            .get("gene_name")
                            .unwrap()
                            .iter()
                            .map(|v| v.clone().unwrap().eq(&String::from("DDX11L2")))
                            .collect::<Vec<bool>>()
                            .iter()
                            .all(|v| *v));

                        assert_eq!(
                            essential
                                .get("transcript_id")
                                .unwrap()
                                .iter()
                                .map(|v| if let Some(id) = v.clone() {
                                    id
                                } else {
                                    String::from("none")
                                })
                                .collect::<Vec<String>>(),
                            vec![
                                String::from("none"),
                                String::from("ENST00000456328.2"),
                                String::from("ENST00000456328.2"),
                                String::from("ENST00000456328.2"),
                                String::from("ENST00000456328.2")
                            ]
                        );
                        assert_eq!(
                            extra
                                .unwrap()
                                .get("gene_type")
                                .unwrap()
                                .iter()
                                .map(|v| if let Some(id) = v.clone() {
                                    id
                                } else {
                                    String::from("none")
                                })
                                .collect::<Vec<String>>(),
                            vec![String::from("lncRNA"); 5]
                        );
                        assert_eq!(tally, 5);
                    }
                }
            }
        }
    }
}
