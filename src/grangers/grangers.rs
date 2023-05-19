use anyhow::{Context,bail};
use std::ops::{Add, Sub, Mul};
use polars::{prelude::*,series::Series,lazy::prelude::*};
use std::fs;
use crate::grangers::reader;
use std::path::Path;
use super::reader::reader_utils::{setdiff, FIELDS};

use super::reader::fasta::SeqInfo;

// GTF files are 1-based with closed intervals.
/// The Grangers struct contains the following fields:
/// - df: the underlying polars dataframe
/// - comments: the comments in the GTF file
/// - chromsize: the chromosome size (set as none before calling `add_chromsize`)
/// - directives: the directives in the GFF file (set as none for GTF files)
#[derive(Clone)]
pub struct Grangers {
    file_type: reader::FileType,
    pub df: DataFrame,
    comments: Vec<String>,
    seqinfo: Option<SeqInfo>,
    directives: Option<Vec<String>>,
    attribute_names: Vec<String>,
}

// IO
impl Grangers {
    pub fn new(df: DataFrame, comments: Vec<String>, directives: Option<Vec<String>>, file_type: reader::FileType) -> Grangers {
        let attribute_names = setdiff(&df.get_column_names()[..], &FIELDS).iter().map(|s| s.to_string()).collect(); 
        Grangers {
            file_type,
            df,
            comments,
            seqinfo: None,
            directives,
            attribute_names,
        }
    }

    // Build the Grangers struct from a GTF file.\
    // Attributes except gene_id, gene_name and transcript_id 
    // will be ignored if only_essential is set to true.\
    // For a human GRh38 GTF file, the extra attributes 
    // can take more than 2GB of memory.
    pub fn from_gtf(file_path: &std::path::Path, only_essential: bool) -> anyhow::Result<Grangers> {
        let am = reader::AttributeMode::from(!only_essential);
        let gstruct = reader::GStruct::from_gtf(file_path, am)?;
        let gr = gstruct.into_grangers()?;
        Ok(gr)
    }

    // Build the Grangers struct from a GFF file.\
    // Attributes except ID, gene_id, gene_name and transcript_id 
    // will be ignored if only_essential is set to true.\
    // For a human GRh38 GTF file, the extra attributes 
    // can take more than 2GB of memory.
    pub fn from_gff(file_path: &std::path::Path, only_essential: bool) -> anyhow::Result<Grangers> {
        let am = reader::AttributeMode::from(!only_essential);
        let gstruct = reader::GStruct::from_gff(file_path, am)?;
        let gr = gstruct.into_grangers()?;
        Ok(gr)
    }

    // add seqinfo to the Grangers struct according to a fasta file
    pub fn add_seqinfo<T: AsRef<Path>>(&mut self, genome_file: T) -> anyhow::Result<()> {
        self.seqinfo = Some(SeqInfo::from_fasta(genome_file)?);
        Ok(())
    }

    pub fn write(&self, file_path: &std::path::Path) -> anyhow::Result<()> {
        fs::create_dir_all(file_path.parent()
            .with_context(|| format!("Could not get the parent directory of the given output file path {:?}",file_path.as_os_str()))?
        )?;
        match self.file_type {
            reader::FileType::GTF => {
                unimplemented!()
            }
            reader::FileType::GFF => {
                unimplemented!()
            }
        }
        // Ok(())
    }
}

// get struct fields
impl Grangers {
    // return the underlying polars dataframe
    pub fn df(&self) -> &DataFrame {
        &self.df
    }

    // return the comments
    pub fn comments(&self) -> &Vec<String> {
        &self.comments
    }

    // return the directives
    pub fn directives(&self) -> Option<&Vec<String>> {
        self.directives.as_ref()
    }

    pub fn seqinfo(&self) -> Option<&SeqInfo> {
        self.seqinfo.as_ref()
    }

    pub fn df_mut(&mut self) -> &mut DataFrame {
        &mut self.df
    }

    // return the comments
    pub fn comments_mut(&mut self) -> &mut Vec<String> {
        &mut self.comments
    }

    // return the directives
    pub fn directives_mut(&mut self) -> Option<&mut Vec<String>> {
        self.directives.as_mut()
    }

    pub fn seqinfo_mut(&mut self) -> Option<&mut SeqInfo> {
        self.seqinfo.as_mut()
    }
}

// get record fields
impl Grangers {
    fn column(&self, col_name: &str) -> anyhow::Result<&Series> {
        self.df
            .column(col_name)
            .with_context(|| format!("Could not get the column {} from the dataframe.", col_name))
    }
    pub fn seqid(&self) -> anyhow::Result<&Series> {
        self.column("seqid")
    }
    pub fn start(&self) -> anyhow::Result<&Series> {
        self.column("start")
    }
    pub fn end(&self) -> anyhow::Result<&Series> {
        self.column("end")
    }
    pub fn strand(&self) -> anyhow::Result<&Series> {
        self.column("strand")
    }
    pub fn score(&self) -> anyhow::Result<&Series> {
        self.column("score")
    }
    pub fn phase(&self) -> anyhow::Result<&Series> {
        self.column("phase")
    }
    pub fn feature_type(&self) -> anyhow::Result<&Series> {
        self.column("feature_type")
    }
    
    pub fn attribute(&self, name: &str) -> anyhow::Result<&Series> {
        if !self.attribute_names.contains(&name.to_string()) {
            bail!("Attribute {} do not exist.", name)
        }

        self.column(name)
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

    pub fn merge_by_gene(&mut self) -> anyhow::Result<()> {
        Ok(())
    }

    pub fn merge_by_transcript(&mut self) -> anyhow::Result<()> {
        Ok(())
    }

    /// flank the genomic features by a given width.
    /// The logic for this - (from the R/GenomicRanges & IRanges packages)

    /// - If `start` is `True` for a given range, the flanking occurs at the start, 
    /// otherwise the end. 
    /// - The `widths` of the flanks are given by the `width` parameter. 
    /// The widths can be negative, in which case the flanking region is 
    /// reversed so that it represents a prefix or suffix of the range. 

    /// Example:
    ///     gr.flank(3, True), where x indicates a range in gr and 
    ///     - indicates the resulting flanking region:\

    ///         ---xxxxxxx\

    ///     If start were FALSE, the range in gr becomes\

    ///         xxxxxxx---\
        
    ///     For negative width, i.e. gr.flank(x, -3, FALSE), 
    ///         where * indicates the overlap between x and the result:\

    ///         xxxx***\

    ///     If both is True, then, for all ranges in x, 
    ///         the flanking regions are extended into 
    ///         (or out of, if width is negative) the range, 
    ///         so that the result straddles the given endpoint 
    ///         and has twice the width given by width. 
    
    ///     This is illustrated below for gr.flank(3, both=TRUE):\

    ///         ---***xxxx\

    /// Args:
    ///     width (int): width to flank by.
    ///     start (bool, optional): only flank starts?. Defaults to True.
    ///     both (bool, optional): both starts and ends?. Defaults to False.
    ///     ignoreStrand (bool, optional): ignore strand?. Defaults to False.

    /// Returns:
    ///     GenomicRanges: a new `GenomicRanges` object with the flanked ranges.
    pub fn flank(
        gr: Grangers,
        width: i64,
        fo: Option<FlankOption>,
    ) -> anyhow::Result<Grangers> {
        let start;
        let both;
        let ignore_strand;
        if let Some(fo) = fo {
            start = fo.start;
            both = fo.both;
            ignore_strand = fo.ignore_strand;
        } else {
            start = true;
            both = false;
            ignore_strand = false;
        }
        let df = gr.df.lazy()
            .with_column(
                when(ignore_strand).then(
                    lit(true)
                ).otherwise(
                    col("strand")
                        .eq(lit("-"))
                        .neq(lit(start))
                )
                .alias("start_flags")
            )
            .with_column(
                    // when both is true
                    when(both)
                        .then(
                            // when start_flag is true
                            when(col("start_flags").eq(lit(true)))
                            .then(col("start") - lit(width).abs())
                            // when start_flag is false
                            .otherwise(col("end") - lit(width).abs() + lit(1)) 
                        )
                    // when both is false
                    .otherwise(
                        // if width >= 0:
                        when(width >= 0).then(
                            // tstart = all_starts[idx] - abs(width) if sf else all_ends[idx] + 1
    
                            when(col("start_flags").eq(lit(true)))
                                .then(
                                    col("start") - lit(width)
                                )
                                .otherwise(col("end") + lit(1))
                        ).otherwise(
                            // tstart = all_starts[idx] if sf else all_ends[idx] + abs(width) + 1
                            when(col("start_flags").eq(lit(true)))
                                .then(col("start"))
                                .otherwise(col("end") + lit(width) + lit(1)))
                        ).alias("start")
            )
            .select([
                // everything except end and start_flags
                all().exclude(["end", "start_flags"]),
                // new_ends.append(tstart + (width * (2 if both else 1) - 1))
                col("start").add(
                    (lit(width).abs().mul(when(lit(both)).then(lit(2)).otherwise(lit(1)))).sub(lit(1))
                ).alias("end")
            ])
            .collect()?;
        Ok(
            Grangers {
                df,
                file_type: gr.file_type,
                comments: gr.comments,
                seqinfo: gr.seqinfo,
                directives: gr.directives,
                attribute_names: gr.attribute_names,
            }
        )
    }

}

#[derive(Copy,Clone)]
pub struct FlankOption {
    start: bool,
    both: bool,
    ignore_strand: bool,
}

impl FlankOption {
    pub fn default() -> FlankOption {
        FlankOption {
            start: true,
            both: false,
            ignore_strand: false,
        }
    }

    pub fn new(start: bool, both: bool, ignore_strand: bool) -> FlankOption {
        FlankOption {
            start,
            both,
            ignore_strand,
        }
    }
}

#[cfg(test)]
mod tests {
    use polars::prelude::*;

    use super::*;
    use crate::gtf::{AttributeMode, FileType,GStruct,Attributes};
    #[test]
    fn test_graners() {
        // let df = df!(
        //     "seqid" => ["chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1"],
        //     "source" => ["HAVANA", "HAVANA", "HAVANA", "HAVANA", "HAVANA", "HAVANA", "HAVANA", "HAVANA", "HAVANA"],
        //     "feature_type" => ["gene", "transcript", "exon", "exon", "exon", "gene", "transcript", "exon", "exon"],
        //     "start" => [1, 1, 1, 21, 41, 101, 101, 101, 121],
        //     "end" => [50, 50, 10, 30, 50, 150, 150, 110,150],
        //     "score" => [100;9],
        //     "strand" => ["+", "+", "+", "+", "+", "-", "-", "-", "-"],
        //     "phase" => [0;9],
        //     "gene_id" =>["g1","g1","g1","g1","g1","g2","g2","g2","g2"],
        //     "gene_name" => ["g1","g1","g1","g1","g1","g2","g2","g2","g2"],
        //     "transcript_id" => [None,Some("t1"),Some("t1"),Some("t1"),None,Some("t2"),Some("t2"),Some("t2"),Some("t2")],
        // ).unwrap();
        // let comments = vec!["comment1".to_string(), "comment2".to_string()];
        // let directives = Some(vec!["directive1".to_string(), "directive2".to_string()]);
        // let file_type = reader::FileType::GTF;
        
        
        let mut gs = GStruct {
            seqid: vec![String::from("chr1");9],
            source: vec![String::from("HAVANA");9],
            feature_type: vec![String::from("gene"), String::from("transcript"), String::from("exon"), String::from("exon"), String::from("exon"), String::from("gene"), String::from("transcript"), String::from("exon"), String::from("exon")],
            start: vec![1, 1, 1, 21, 41, 101, 101, 101, 121],
            end: vec![50, 50, 10, 30, 50, 150, 150, 110,150],
            score: vec![Some(10.0);9],
            strand: vec![Some(String::from("+")), Some(String::from("+")), Some(String::from("+")), Some(String::from("+")), Some(String::from("+")), Some(String::from("-")), Some(String::from("-")), Some(String::from("-")), Some(String::from("-"))],
            phase: vec![Some(String::from("0"));9],
            attributes: Attributes::new(AttributeMode::Full,FileType::GTF).unwrap(),
            comments: vec!["comment1".to_string(), "comment2".to_string()],
            directives: Some(vec!["directive1".to_string(), "directive2".to_string()]),
        };
        let gsr = &mut gs;
        gsr.attributes.file_type = FileType::GTF;
        gsr.attributes.essential.insert(String::from("gene_id"), vec![Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g2")),Some(String::from("g2")),Some(String::from("g2")),Some(String::from("g2"))]);
        gsr.attributes.essential.insert(String::from("gene_name"), vec![Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g2")),Some(String::from("g2")),Some(String::from("g2")),Some(String::from("g2"))]);
        gsr.attributes.essential.insert(String::from("transcript_id"), vec![None,Some(String::from("t1")),Some(String::from(String::from("t1"))),Some(String::from("t1")),None,Some(String::from("t2")),Some(String::from("t2")),Some(String::from("t2")),Some(String::from("t2"))]);


        if let Some(extra) = &mut gsr.attributes.extra {
            extra.insert(String::from("gene_version"), vec![Some(String::from("1"));9]);
        }

        let mut gr = gs.into_grangers().unwrap();

        // test builder
        // assert_eq!(gr.df(), &df);
        // assert_eq!(gr.comments(), &comments);
        // assert_eq!(gr.directives(), directives.as_ref());
        // assert!(gr.file_type.is_gtf());
    }

    #[test]
    fn test_flank() {

        let mut gr = get_toy_gr().unwrap();

        // test flank with default parameters
        let fo = FlankOption {
            start: true,
            both: false,
            ignore_strand: false,
        };
        let gr1 = Grangers::flank(gr.clone(), 10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![91,91,91,111,131,251,251,211,251];
        let end: Vec<i64> = vec![100,100,100,120,140,260,260,220,260];
        assert_eq!(gr1.start().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), start);
        assert_eq!(gr1.end().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), end);

        let gr1 = Grangers::flank(gr.clone(), -10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![101,101,101,121,141,241,241,201,241];
        let end: Vec<i64> = vec![110,110,110,130,150,250,250,210,250];
        assert_eq!(gr1.start().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), start);
        assert_eq!(gr1.end().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), end);

        // test flank with default parameters and both=true
        let fo = FlankOption {
            start: true,
            both: true,
            ignore_strand: false,
        };
        let gr1 = Grangers::flank(gr.clone(), 10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![91,91,91,111,131,241,241,201,241];
        let end: Vec<i64> = vec![110,110,110,130,150,260,260,220,260];
        assert_eq!(gr1.start().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), start);
        assert_eq!(gr1.end().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), end);

        let gr1 = Grangers::flank(gr.clone(), -10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![91,91,91,111,131,241,241,201,241];
        let end: Vec<i64> = vec![110,110,110,130,150,260,260,220,260];
        assert_eq!(gr1.start().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), start);
        assert_eq!(gr1.end().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), end);

        // test flank with start = false
        let fo = FlankOption {
            start: false,
            both: false,
            ignore_strand: false,
        };
        let gr1 = Grangers::flank(gr.clone(), 10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![151,151,111,131,151,191,191,191,211];
        let end: Vec<i64> = vec![160,160,120,140,160,200,200,200,220];
        assert_eq!(gr1.start().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), start);
        assert_eq!(gr1.end().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), end);

        let gr1 = Grangers::flank(gr.clone(), -10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![141,141,101,121,141,201,201,201,221];
        let end: Vec<i64> = vec![150,150,110,130,150,210,210,210,230];
        assert_eq!(gr1.start().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), start);
        assert_eq!(gr1.end().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), end);

        // test flank with start = false and both=true
        let fo = FlankOption {
            start: false,
            both: true,
            ignore_strand: false,
        };
        let gr1 = Grangers::flank(gr.clone(), 10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![141,141,101,121,141,191,191,191,211];
        let end: Vec<i64> = vec![160,160,120,140,160,210,210,210,230];
        assert_eq!(gr1.start().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), start);
        assert_eq!(gr1.end().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), end);

        let gr1 = Grangers::flank(gr.clone(), -10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![141,141,101,121,141,191,191,191,211];
        let end: Vec<i64> = vec![160,160,120,140,160,210,210,210,230];
        assert_eq!(gr1.start().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), start);
        assert_eq!(gr1.end().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), end);

        // test flank with ignore_strand: true
        let fo = FlankOption {
            start: true,
            both: false,
            ignore_strand: true,
        };
        let gr1 = Grangers::flank(gr.clone(), 10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![91,91,91,111,131,191,191,191,211];
        let end: Vec<i64> = vec![100,100,100,120,140,200,200,200,220];
        assert_eq!(gr1.start().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), start);
        assert_eq!(gr1.end().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), end);

        let gr1 = Grangers::flank(gr.clone(), -10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![101,101,101,121,141,201,201,201,221];
        let end: Vec<i64> = vec![110,110,110,130,150,210,210,210,230];
        assert_eq!(gr1.start().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), start);
        assert_eq!(gr1.end().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), end);

        // test flank with ignore_strand: true and both=true
        let fo = FlankOption {
            start: true,
            both: true,
            ignore_strand: true,
        };
        let gr1 = Grangers::flank(gr.clone(), 10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![91,91,91,111,131,191,191,191,211];
        let end: Vec<i64> = vec![110,110,110,130,150,210,210,210,230];
        assert_eq!(gr1.start().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), start);
        assert_eq!(gr1.end().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), end);

        let gr1 = Grangers::flank(gr.clone(), -10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![91,91,91,111,131,191,191,191,211];
        let end: Vec<i64> = vec![110,110,110,130,150,210,210,210,230];
        assert_eq!(gr1.start().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), start);
        assert_eq!(gr1.end().unwrap().i64().unwrap().to_vec().into_iter().map(|x|x.unwrap()).collect::<Vec<i64>>(), end);

    }

    fn get_toy_gr() -> anyhow::Result<Grangers> {

    let mut gs = GStruct {
        seqid: vec![String::from("chr1");9],
        source: vec![String::from("HAVANA");9],
        feature_type: vec![String::from("gene"), String::from("transcript"), String::from("exon"), String::from("exon"), String::from("exon"), String::from("gene"), String::from("transcript"), String::from("exon"), String::from("exon")],
        start: vec![101, 101, 101, 121,141, 201, 201, 201, 221],
        end: vec![150, 150, 110, 130, 150, 250, 250, 210,250],
        score: vec![Some(10.0);9],
        strand: vec![Some(String::from("+")), Some(String::from("+")), Some(String::from("+")), Some(String::from("+")), Some(String::from("+")), Some(String::from("-")), Some(String::from("-")), Some(String::from("-")), Some(String::from("-"))],
        phase: vec![Some(String::from("0"));9],
        attributes: Attributes::new(AttributeMode::Full,FileType::GTF)?,
        comments: vec!["comment1".to_string(), "comment2".to_string()],
        directives: Some(vec!["directive1".to_string(), "directive2".to_string()]),
    };
    let gsr = &mut gs;
    gsr.attributes.file_type = FileType::GTF;
    gsr.attributes.essential.insert(String::from("gene_id"), vec![Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g2")),Some(String::from("g2")),Some(String::from("g2")),Some(String::from("g2"))]);
    gsr.attributes.essential.insert(String::from("gene_name"), vec![Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g1")),Some(String::from("g2")),Some(String::from("g2")),Some(String::from("g2")),Some(String::from("g2"))]);
    gsr.attributes.essential.insert(String::from("transcript_id"), vec![None,Some(String::from("t1")),Some(String::from(String::from("t1"))),Some(String::from("t1")),None,Some(String::from("t2")),Some(String::from("t2")),Some(String::from("t2")),Some(String::from("t2"))]);

    if let Some(extra) = &mut gsr.attributes.extra {
        extra.insert(String::from("gene_version"), vec![Some(String::from("1"));9]);
    }

    let mut gr = gs.into_grangers()?;
    Ok(gr)

    }
}