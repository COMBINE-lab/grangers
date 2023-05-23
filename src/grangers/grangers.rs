use super::reader::fasta::SeqInfo;
use super::reader::reader_utils::{setdiff, FIELDS};
use crate::grangers::reader;
use anyhow::{bail, Context};
use polars::frame::row::Row;
use polars::{lazy::prelude::*, prelude::*, series::Series};
use rust_lapper::{Interval, Lapper};
use std::collections::HashMap;
use std::fs;
use std::iter::FromIterator;
use std::path::Path;
use std::time::Instant;
use std::{
    collections::HashSet,
    ops::{Add, Mul, Sub},
};
use tracing::{warn, info};

const ESSENTIAL_FIELDS: [&str; 3] = ["seqnames", "start", "end"];

// GTF files are 1-based with closed intervals.
/// The Grangers struct contains the following fields:
/// - df: the underlying polars dataframe
/// - comments: the comments in the GTF file
/// - chromsize: the chromosome size (set as none before calling `add_chromsize`)
/// - directives: the directives in the GFF file (set as none for GTF files)
#[derive(Clone)]
pub struct Grangers {
    // file_type: Option<reader::FileType>,
    df: DataFrame,
    misc: Option<HashMap<String, Vec<String>>>,
    // comments: Option<Vec<String>>,
    seqinfo: Option<SeqInfo>,
    // directives: Option<Vec<String>>,
    // attribute_names: Option<Vec<String>>,
    lapper: Option<Lapper<u64, Vec<String>>>
}

// IO
impl Grangers {
    fn add_column<T: SeriesTrait>(&mut self, series: Series) -> anyhow::Result<()> {
        if let Err(_) = self.column("strand") {
            warn!("The dataframe does not contain the required column seqnames; Cannot proceed");
            self.df.with_column(series)?;
        }
        Ok(())
    }
    fn validate(df: &mut DataFrame) -> anyhow::Result<()> {
        // check essential fields
        for ef in ESSENTIAL_FIELDS {
            if let Err(_) = df.column(ef) {
                bail!(
                    "The dataframe does not contain the required column {}; Cannot proceed",
                    ef
                );
            }
        }

        // check strand
        if let Err(_) = df.column("strand") {
            warn!("The dataframe does not contain the required column seqnames; Cannot proceed");
            df.with_column(Series::new_null("predictions", df.height()))?;
        }

        Ok(())
    }
    /// Instantiate a new Grangers struct according to
    /// - a range dataframe that contain the ranges of the genomic features
    /// - an optional SeqInfo struct contains the reference information (usually chromosome)
    /// - an optional HashMap<String, Vec<String>>> contains additional information.
    /// The datafram should contain the ranges of the genomic features, and the SeqInfo struct should contain the chromosome sizes. The  
    pub fn new(
        mut df: DataFrame,
        seqinfo: Option<SeqInfo>,
        misc: Option<HashMap<String, Vec<String>>>,
        lapper: Option<Lapper<u64, Vec<String>>>,
    ) -> anyhow::Result<Grangers> {
        // validate the dataframe
        Grangers::validate(&mut df)?;

        Ok(Grangers { df, misc, seqinfo, lapper })
    }

    
    pub fn from_gstruct(gstruct: reader::GStruct) -> anyhow::Result<Grangers> {
        // create dataframe!
        // we want to make some columns categorical because of this https://docs.rs/polars/latest/polars/docs/performance/index.html
        // fields
        let mut df_vec = vec![
            Series::new("seqnames", gstruct.seqid),
            // .cast(&DataType::Categorical(None)).unwrap(),
            Series::new("source", gstruct.source),
            // .cast(&DataType::Categorical(None)).unwrap(),
            Series::new("feature_type", gstruct.feature_type),
            // .cast(&DataType::Categorical(None)).unwrap(),
            Series::new("start", gstruct.start),
            Series::new("end", gstruct.end),
            Series::new("score", gstruct.score),
            Series::new("strand", gstruct.strand),
            // .cast(&DataType::Categorical(None)).unwrap()
            Series::new("phase", gstruct.phase),
            // .cast(&DataType::Categorical(None)).unwrap(),
        ];

        //for essential attributes
        for (k, v) in gstruct.attributes.essential {
            let s = Series::new(k.as_str(), v);
            // .cast(&DataType::Categorical(None)).unwrap();
            df_vec.push(s);
        }

        // for extra attributes
        if let Some(attributes) = gstruct.attributes.extra {
            for (k, v) in attributes {
                df_vec.push(Series::new(k.as_str(), v))
            }
        }
        let df = DataFrame::new(df_vec)?;
        let gr = Grangers::new(df, None, gstruct.misc, None)?;
        Ok(gr)
    }


    // Build the Grangers struct from a GTF file.\
    // Attributes except gene_id, gene_name and transcript_id
    // will be ignored if only_essential is set to true.\
    // For a human GRh38 GTF file, the extra attributes
    // can take more than 2GB of memory.
    pub fn from_gtf(file_path: &std::path::Path, only_essential: bool) -> anyhow::Result<Grangers> {
        let am = reader::AttributeMode::from(!only_essential);
        let gstruct = reader::GStruct::from_gtf(file_path, am)?;
        let gr = Grangers::from_gstruct(gstruct)?;
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
        let gr = Grangers::from_gstruct(gstruct)?;
        Ok(gr)
    }

    // add seqinfo to the Grangers struct according to a fasta file
    pub fn add_seqinfo<T: AsRef<Path>>(&mut self, genome_file: T) -> anyhow::Result<()> {
        self.seqinfo = Some(SeqInfo::from_fasta(genome_file)?);
        Ok(())
    }

    pub fn write(&self, file_path: &std::path::Path) -> anyhow::Result<()> {
        fs::create_dir_all(file_path.parent().with_context(|| {
            format!(
                "Could not get the parent directory of the given output file path {:?}",
                file_path.as_os_str()
            )
        })?)?;
        // match self.file_type {
        //     reader::FileType::GTF => {
        //         unimplemented!()
        //     }
        //     reader::FileType::GFF => {
        //         unimplemented!()
        //     }
        // }
        unimplemented!()

        // Ok(())
    }
}

// get struct fields
impl Grangers {
    // return the underlying polars dataframe
    pub fn df(&self) -> &DataFrame {
        &self.df
    }

    pub fn seqinfo(&self) -> Option<&SeqInfo> {
        self.seqinfo.as_ref()
    }

    pub fn df_mut(&mut self) -> &mut DataFrame {
        &mut self.df
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
    pub fn lapper(&self) -> &Option<Lapper<u64, Vec<String>>> {
        &self.lapper
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

    /// The function consumes a Grangers object, flank the genomic features of it by a given width, and returns a new Grangers.
    /// The logic for flanking - (from the BiocPy/GenomicRanges)

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
    ///     gr: A `Grangers` object.
    ///     width (int): width to flank by.
    ///     flank_option (FlankOption, optional):
    ///     start (bool, optional): only flank starts?. Defaults to True.
    ///     both (bool, optional): both starts and ends?. Defaults to False.
    ///     ignoreStrand (bool, optional): ignore strand?. Defaults to False.

    /// Returns:
    ///     Grangers: a new `Grangers` object with the flanked ranges.
    pub fn flank(
        gr: Grangers,
        width: i64,
        option: Option<FlankOption>,
    ) -> anyhow::Result<Grangers> {
        let start;
        let both;
        let ignore_strand;
        if let Some(option) = option {
            start = option.start;
            both = option.both;
            ignore_strand = option.ignore_strand;
        } else {
            start = true;
            both = false;
            ignore_strand = false;
        }
        let df = gr
            .df
            .lazy()
            .with_column(
                when(ignore_strand)
                    .then(lit(true))
                    .otherwise(col("strand").eq(lit("-")).neq(lit(start)))
                    .alias("start_flags"),
            )
            .with_column(
                // when both is true
                when(both)
                    .then(
                        // when start_flag is true
                        when(col("start_flags").eq(lit(true)))
                            .then(col("start") - lit(width).abs())
                            // when start_flag is false
                            .otherwise(col("end") - lit(width).abs() + lit(1)),
                    )
                    // when both is false
                    .otherwise(
                        // if width >= 0:
                        when(width >= 0)
                            .then(
                                // tstart = all_starts[idx] - abs(width) if sf else all_ends[idx] + 1
                                when(col("start_flags").eq(lit(true)))
                                    .then(col("start") - lit(width))
                                    .otherwise(col("end") + lit(1)),
                            )
                            .otherwise(
                                // tstart = all_starts[idx] if sf else all_ends[idx] + abs(width) + 1
                                when(col("start_flags").eq(lit(true)))
                                    .then(col("start"))
                                    .otherwise(col("end") + lit(width) + lit(1)),
                            ),
                    )
                    .alias("start"),
            )
            .select([
                // everything except end and start_flags
                all().exclude(["end", "start_flags"]),
                // new_ends.append(tstart + (width * (2 if both else 1) - 1))
                col("start")
                    .add(
                        (lit(width)
                            .abs()
                            .mul(when(lit(both)).then(lit(2)).otherwise(lit(1))))
                        .sub(lit(1)),
                    )
                    .alias("end"),
            ])
            .collect()?;
        Ok(Grangers {
            df,
            seqinfo: gr.seqinfo,
            misc: gr.misc,
            lapper: gr.lapper,
        })
    }

    /// find the gaps between records in each group by a given column name used for finding groups. This function will call the `setdiff` function internally after setting the boundary of each group as the range from the smallest `start` value to the largest `end` value.
    /// To make the most sense, the `by` argument should be a column representing the name of high level features, for example, genes or transcripts.
    pub fn gaps(&self, by: &str) -> anyhow::Result<Grangers> {
        Ok(self.clone())
    }

    /// Find the set difference of genomic intervals with `other`.
    /// The `on` and `boundary_on` arguments are used as anchors to find the corresponding intervals in `self` and `boundary` respectively.
    /// For the boundary dataframe, all values in the `boundary_on` column should be unique, because it is used for defining the boundary of each feature.
    pub fn setdiff(
        &self,
        boundary: Grangers,
        on: &str,
        boundary_on: &str,
    ) -> anyhow::Result<Grangers> {
        // check if the `on` column exists in both objects
        if !self.df().get_column_names().contains(&on) {
            bail!(
                "The `on` column {} does not exist in the Grangers object.",
                on
            );
        }

        if boundary.df().get_column_names().contains(&boundary_on) {
            bail!(
                "The `boundary_on` column {} does not exist in the boundary object.",
                boundary_on
            );
        }

        // check if the boundary dataframe contains only one range for each feature
        if boundary.df().height() != boundary.df().column(boundary_on)?.unique()?.len() {
            bail!("The boundary dataframe contains more than one range for each feature.");
        }

        // check if the two `on` columns overlap
        let on_vec: Vec<String> = self
            .df()
            .column(on)?
            .utf8()?
            .into_iter()
            .map(|val| val.unwrap().to_string())
            .collect();

        let a = self
            .df()
            .select(["seqname", "start", "end", "strand", on])?
            .left_join(boundary.df(), [on], [boundary_on])?;

        Ok(self.clone())
    }

    /// this function find the boundary of each feature in the Grangers object by a given column name used for finding groups. The boundary of each group are defined as from the smallest `start` value to the largest `end` value.
    pub fn boundary(&self, by: &str) -> anyhow::Result<Grangers> {
        Ok(self.clone())
    }

    /// this function turns the seqinfo of the Grangers object into a boundary Grangers object.
    pub fn seqinfo_as_bounary(&self) -> anyhow::Result<Grangers> {
        Ok(self.clone())
    }

    /// merge the features by the given columns via the `by` argument to generate a new Grangers object.
    /// *** Argument:
    /// - `by`: a vector of string representing which group(s) to merge by. Each string should be a valid column name.
    /// - `ignore_strand`: whether to ignore the strand information when merging.
    /// - `slack`: the maximum distance between two features to be merged.
    pub fn merge<I, S>(
        &mut self,
        mo: MergeOption,
    ) -> anyhow::Result<Grangers> {
        // check if the `by` vector contains valid column namess
        for col_name in mo.by.iter() {
            if self.is_column(col_name.as_ref()) {
                bail!(
                    " `by` contains non-existing column - {}. Cannot proceed",
                    col_name
                );
            }
        }
        self.build_lapper(&mo.by)?;

        let mut lapper = if let Some(lapper) = self.lapper.clone() {
            lapper
        } else {
            bail!("Could not find the lapper that was just built; Please report this bug!")
        }; 

        lapper.merge_overlaps();

        if !lapper.overlaps_merged {
            info!("Did not find overlapping features. Nothing to merge.")
        }

        let lapper_vecs = LapperVecs::new(&lapper, &mo.by);
        let df = lapper_vecs.to_df()?;
        Grangers::new(df, self.seqinfo.clone(), self.misc.clone(), Some(lapper))
    }

    /// Instantiate a Lapper struct for interval search and overlap detection
    /// The Lapper struct works for unsigned integer only, so the start and end values will be converted to unsigned integer.
    /// The range is [start, end), i.e. inclusive start and exclusive end. Note that for GTF/GFF, the start is 1-based and the end is closed; for BED, the start is 0-based and open but the end is closed.
    pub fn build_lapper(
        &mut self,
        meta_cols: &Vec<String>,
    ) -> anyhow::Result<()> {
        // rust_lapper
        // only unsigned
        // [start, stop) Inclusive start, exclusive of stop
        type Iv = Interval<u64, Vec<String>>;

        let df = self.df();
        let mut meta_vec = vec![Vec::<String>::with_capacity(meta_cols.len()); df.height()];
        let mut iters = df
            .columns(meta_cols)?
            .iter()
            .map(|s| s.iter())
            .collect::<Vec<_>>();

        // build meta_vec
        for row in 0..df.height() {
            for iter in &mut iters {
                let value = iter
                    .next()
                    .expect("should have as many iterations as rows")
                    // .cast(&DataType::Utf8)? // this will cause a bug saying cannot cast non numeric to numeric
                    .to_string();
                // process value
                if let Some(v) = meta_vec.get_mut(row) {
                    v.push(value);
                } else {
                    bail!("meta_vec length is not equal to the dataframe height! Please Reort this bug!")
                }
            }
        }

        // build a lapper from the dataframe
        // check if negative start or end
        let min_start = df
            .column("start")?
            .cast(&DataType::Int64)?
            .i64()?
            .min()
            .with_context(|| {
                format!("Could not get the minimum start value from the dataframe.")
            })?;

        let mut pos_vec = vec![Vec::<String>::with_capacity(meta_cols.len()); df.height()];
        let mut iters = df
            .columns(["start", "end"])?
            .iter()
            .map(|s| s.iter())
            .collect::<Vec<_>>();

        let mut lapper_tree_vec = Vec::with_capacity(df.height());
        // build pos_vec
        for (rid, meta) in meta_vec.into_iter().enumerate() {
            let s: i64 = iters[0]
                .next()
                .expect("should have as many iterations as rows")
                .cast(&DataType::Int64)?
                .try_extract()?;
            let e: i64 = iters[1]
                .next()
                .expect("should have as many iterations as rows")
                .cast(&DataType::Int64)?
                .try_extract()?;

            // lappers take unsigned integer
            let (start, end) = if !min_start.is_negative() {
                (s as u64, e as u64)
            } else {
                ((s + min_start.abs()) as u64, (e + min_start.abs()) as u64)
            };

            lapper_tree_vec.push(Iv {
                start,
                stop: end,
                val: meta,
            });
        }

        // rust-lapper
        let start = std::time::Instant::now();
        self.lapper = Some(Lapper::new(lapper_tree_vec));
        let duration: std::time::Duration = start.elapsed();
        info!("build rust-lapper in {:?}", duration);
        Ok(())
    }

    /// check if the GRanges object has a column of a given name
    fn is_column(&self, col_name: &str) -> bool {
        let df = self.df();
        let mut valid = false;
        if !df.get_column_names().contains(&col_name) {
            valid = true;
        }
        valid
    }
}

struct LapperVecs {
    start: Vec<i64>,
    end: Vec<i64>,
    val: Vec<Vec<String>>,
    val_names: Vec<String>,
}

impl LapperVecs {
    fn new(lapper: &Lapper<u64, Vec<String>>, meta_cols: &Vec<String>) -> LapperVecs {
        // initialize the vectors
        let mut start: Vec<i64> = Vec::with_capacity(lapper.len());
        let mut end: Vec<i64> = Vec::with_capacity(lapper.len());
        let mut val: Vec<Vec<String>> = vec![Vec::<String>::with_capacity(lapper.len()); meta_cols.len()];

        // build the vectors
        for iv in lapper.iter() {
            start.push(iv.start as i64);
            end.push(iv.stop as i64);
            for (i, v) in iv.val.iter().enumerate() {
                val.get_mut(i).expect("lapper vec val is shorter than designed").push(v.to_string());
            }
        }

        
        LapperVecs {
            start,
            end,
            val,
            val_names: meta_cols.clone(),
        }
    }

    fn to_df(self) -> anyhow::Result<DataFrame> {
        let mut df_vec = vec![
            Series::new("start", self.start),
            Series::new("end", self.end),
        ];

        for (name, value) in self.val_names.into_iter().zip(self.val.into_iter()) {
            df_vec.push(Series::new(name.as_str(), value));
        }

        let df = DataFrame::new(df_vec)?;
        Ok(df)

    }
}

#[derive(Copy, Clone)]
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

pub struct MergeOption {
    pub by: Vec<String>,
    pub slack: usize,
}

impl MergeOption {
    pub fn default() -> MergeOption {
        MergeOption {
            by: vec![String::from("seqnames"), String::from("strand")],
            slack: 0,
        }
    }

    pub fn new(by: Vec<String>, ignore_strand: bool, slack: usize) -> MergeOption {

        // avoid duplicated columns
        let mut by_hash: HashSet<String> = by.into_iter().collect();
        
        if !ignore_strand {
            by_hash.insert(String::from("strand"));
        }


        // add chromosome name and strand if needed
        if by_hash.insert(String::from("seqnames")) {
            warn!("Added `seqnames` to the provided `by` vector as it is missing. Please instantiate it manually if you intend to not include it.")
        };

        MergeOption {
            by: by_hash.into_iter().collect(),
            slack,
        }
    }
}

#[cfg(test)]
mod tests {
    use polars::prelude::*;

    use super::*;
    use crate::gtf::{AttributeMode, Attributes, FileType, GStruct};
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
            seqid: vec![String::from("chr1"); 9],
            source: vec![String::from("HAVANA"); 9],
            feature_type: vec![
                String::from("gene"),
                String::from("transcript"),
                String::from("exon"),
                String::from("exon"),
                String::from("exon"),
                String::from("gene"),
                String::from("transcript"),
                String::from("exon"),
                String::from("exon"),
            ],
            start: vec![1, 1, 1, 21, 41, 101, 101, 101, 121],
            end: vec![50, 50, 10, 30, 50, 150, 150, 110, 150],
            score: vec![Some(10.0); 9],
            strand: vec![
                Some(String::from("+")),
                Some(String::from("+")),
                Some(String::from("+")),
                Some(String::from("+")),
                Some(String::from("+")),
                Some(String::from("-")),
                Some(String::from("-")),
                Some(String::from("-")),
                Some(String::from("-")),
            ],
            phase: vec![Some(String::from("0")); 9],
            attributes: Attributes::new(AttributeMode::Full, FileType::GTF).unwrap(),
            misc: Some(HashMap::new()),
        };
        let gsr = &mut gs;
        gsr.attributes.file_type = FileType::GTF;
        gsr.attributes.essential.insert(
            String::from("gene_id"),
            vec![
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g2")),
                Some(String::from("g2")),
                Some(String::from("g2")),
                Some(String::from("g2")),
            ],
        );
        gsr.attributes.essential.insert(
            String::from("gene_name"),
            vec![
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g2")),
                Some(String::from("g2")),
                Some(String::from("g2")),
                Some(String::from("g2")),
            ],
        );
        gsr.attributes.essential.insert(
            String::from("transcript_id"),
            vec![
                None,
                Some(String::from("t1")),
                Some(String::from(String::from("t1"))),
                Some(String::from("t1")),
                None,
                Some(String::from("t2")),
                Some(String::from("t2")),
                Some(String::from("t2")),
                Some(String::from("t2")),
            ],
        );

        if let Some(extra) = &mut gsr.attributes.extra {
            extra.insert(
                String::from("gene_version"),
                vec![Some(String::from("1")); 9],
            );
        }
        let gr = Grangers::from_gstruct(gs).unwrap();


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
        let start: Vec<i64> = vec![91, 91, 91, 111, 131, 251, 251, 211, 251];
        let end: Vec<i64> = vec![100, 100, 100, 120, 140, 260, 260, 220, 260];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        let gr1 = Grangers::flank(gr.clone(), -10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![101, 101, 101, 121, 141, 241, 241, 201, 241];
        let end: Vec<i64> = vec![110, 110, 110, 130, 150, 250, 250, 210, 250];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        // test flank with default parameters and both=true
        let fo = FlankOption {
            start: true,
            both: true,
            ignore_strand: false,
        };
        let gr1 = Grangers::flank(gr.clone(), 10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![91, 91, 91, 111, 131, 241, 241, 201, 241];
        let end: Vec<i64> = vec![110, 110, 110, 130, 150, 260, 260, 220, 260];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        let gr1 = Grangers::flank(gr.clone(), -10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![91, 91, 91, 111, 131, 241, 241, 201, 241];
        let end: Vec<i64> = vec![110, 110, 110, 130, 150, 260, 260, 220, 260];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        // test flank with start = false
        let fo = FlankOption {
            start: false,
            both: false,
            ignore_strand: false,
        };
        let gr1 = Grangers::flank(gr.clone(), 10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![151, 151, 111, 131, 151, 191, 191, 191, 211];
        let end: Vec<i64> = vec![160, 160, 120, 140, 160, 200, 200, 200, 220];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        let gr1 = Grangers::flank(gr.clone(), -10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![141, 141, 101, 121, 141, 201, 201, 201, 221];
        let end: Vec<i64> = vec![150, 150, 110, 130, 150, 210, 210, 210, 230];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        // test flank with start = false and both=true
        let fo = FlankOption {
            start: false,
            both: true,
            ignore_strand: false,
        };
        let gr1 = Grangers::flank(gr.clone(), 10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![141, 141, 101, 121, 141, 191, 191, 191, 211];
        let end: Vec<i64> = vec![160, 160, 120, 140, 160, 210, 210, 210, 230];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        let gr1 = Grangers::flank(gr.clone(), -10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![141, 141, 101, 121, 141, 191, 191, 191, 211];
        let end: Vec<i64> = vec![160, 160, 120, 140, 160, 210, 210, 210, 230];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        // test flank with ignore_strand: true
        let fo = FlankOption {
            start: true,
            both: false,
            ignore_strand: true,
        };
        let gr1 = Grangers::flank(gr.clone(), 10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![91, 91, 91, 111, 131, 191, 191, 191, 211];
        let end: Vec<i64> = vec![100, 100, 100, 120, 140, 200, 200, 200, 220];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        let gr1 = Grangers::flank(gr.clone(), -10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![101, 101, 101, 121, 141, 201, 201, 201, 221];
        let end: Vec<i64> = vec![110, 110, 110, 130, 150, 210, 210, 210, 230];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        // test flank with ignore_strand: true and both=true
        let fo = FlankOption {
            start: true,
            both: true,
            ignore_strand: true,
        };
        let gr1 = Grangers::flank(gr.clone(), 10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![91, 91, 91, 111, 131, 191, 191, 191, 211];
        let end: Vec<i64> = vec![110, 110, 110, 130, 150, 210, 210, 210, 230];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );

        let gr1 = Grangers::flank(gr.clone(), -10, Some(fo)).unwrap();
        let start: Vec<i64> = vec![91, 91, 91, 111, 131, 191, 191, 191, 211];
        let end: Vec<i64> = vec![110, 110, 110, 130, 150, 210, 210, 210, 230];
        assert_eq!(
            gr1.start()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            start
        );
        assert_eq!(
            gr1.end()
                .unwrap()
                .i64()
                .unwrap()
                .to_vec()
                .into_iter()
                .map(|x| x.unwrap())
                .collect::<Vec<i64>>(),
            end
        );
    }

    fn get_toy_gr() -> anyhow::Result<Grangers> {
        let mut gs = GStruct {
            seqid: vec![String::from("chr1"); 9],
            source: vec![String::from("HAVANA"); 9],
            feature_type: vec![
                String::from("gene"),
                String::from("transcript"),
                String::from("exon"),
                String::from("exon"),
                String::from("exon"),
                String::from("gene"),
                String::from("transcript"),
                String::from("exon"),
                String::from("exon"),
            ],
            start: vec![101, 101, 101, 121, 141, 201, 201, 201, 221],
            end: vec![150, 150, 110, 130, 150, 250, 250, 210, 250],
            score: vec![Some(10.0); 9],
            strand: vec![
                Some(String::from("+")),
                Some(String::from("+")),
                Some(String::from("+")),
                Some(String::from("+")),
                Some(String::from("+")),
                Some(String::from("-")),
                Some(String::from("-")),
                Some(String::from("-")),
                Some(String::from("-")),
            ],
            phase: vec![Some(String::from("0")); 9],
            attributes: Attributes::new(AttributeMode::Full, FileType::GTF)?,
            misc: Some(HashMap::new()),
        };
        let gsr = &mut gs;
        gsr.attributes.file_type = FileType::GTF;
        gsr.attributes.essential.insert(
            String::from("gene_id"),
            vec![
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g2")),
                Some(String::from("g2")),
                Some(String::from("g2")),
                Some(String::from("g2")),
            ],
        );
        gsr.attributes.essential.insert(
            String::from("gene_name"),
            vec![
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g1")),
                Some(String::from("g2")),
                Some(String::from("g2")),
                Some(String::from("g2")),
                Some(String::from("g2")),
            ],
        );
        gsr.attributes.essential.insert(
            String::from("transcript_id"),
            vec![
                None,
                Some(String::from("t1")),
                Some(String::from(String::from("t1"))),
                Some(String::from("t1")),
                None,
                Some(String::from("t2")),
                Some(String::from("t2")),
                Some(String::from("t2")),
                Some(String::from("t2")),
            ],
        );

        if let Some(extra) = &mut gsr.attributes.extra {
            extra.insert(
                String::from("gene_version"),
                vec![Some(String::from("1")); 9],
            );
        }

        let mut gr = Grangers::from_gstruct(gs)?;
        Ok(gr)
    }
}
