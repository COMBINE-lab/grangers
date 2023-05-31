use anyhow::bail;
use polars::prelude::DataFrame;
use std::collections::HashSet;
use tracing::warn;

#[derive(Copy, Clone)]
pub struct FlankOptions {
    pub start: bool,
    pub both: bool,
    pub ignore_strand: bool,
}

impl Default for FlankOptions {
    fn default() -> FlankOptions {
        FlankOptions {
            start: true,
            both: false,
            ignore_strand: false,
        }
    }
}

impl FlankOptions {
    pub fn new(start: bool, both: bool, ignore_strand: bool) -> FlankOptions {
        FlankOptions {
            start,
            both,
            ignore_strand,
        }
    }
}

/// Options used for the merge function
/// - by: a vector of string representing which column(s) to merge by. Each string should be a valid column name.
/// - slack: the minimum gap between two features to be merged. It should be a positive integer.
/// - output_count: whether to output the count of ranges of the merged features.
pub struct MergeOptions {
    pub by: Vec<String>,
    pub slack: i64,
    pub ignore_strand: bool,
}

impl Default for MergeOptions {
    fn default() -> MergeOptions {
        MergeOptions {
            by: vec![String::from("seqname"), String::from("strand")],
            slack: 1,
            ignore_strand: false,
        }
    }
}

impl MergeOptions {
    pub fn new<T: AsRef<str>>(
        by: &[T],
        ignore_strand: bool,
        slack: i64,
    ) -> anyhow::Result<MergeOptions> {
        // avoid duplicated columns
        let mut by_hash: HashSet<String> = by.into_iter().map(|n| n.as_ref().to_string()).collect();

        if slack < 1 {
            warn!("It usually doen't make sense to set a non-positive slack.")
        }

        if by_hash.take(&String::from("start")).is_some()
            | by_hash.take(&String::from("end")).is_some()
        {
            bail!("The provided `by` vector cannot contain the start or end column")
        };

        if ignore_strand {
            if by_hash.take(&String::from("strand")).is_some() {
                warn!("Remove `strand` from the provided `by` vector as the ignored_strand flag is set.")
            }
        } else {
            by_hash.insert(String::from("strand"));
        }

        // add chromosome name and strand if needed
        if by_hash.insert(String::from("seqname")) {
            warn!("Added `seqname` to the `by` vector as it is required.")
        };

        Ok(MergeOptions {
            by: by_hash.into_iter().collect(),
            slack,
            ignore_strand,
        })
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum ExtendOption {
    /// Extend the feature to the start
    Start,
    /// Extend the feature to the end
    End,
    /// Extend the feature to both sides
    Both,
}

pub struct GetSequenceOptions {}

/// Options used for dealing with out-of-boundary features
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum OOBOption {
    Truncate,
    Skip,
}

/// The name of columns used to store the information of features.
/// It is designed according to the GTF/GFF format.
/// https://useast.ensembl.org/info/website/upload/gff.html
/// The default values are the column names used in GTF/GFF files.
/// This will be used in almost all Grangers methods.
#[derive(Clone, PartialEq, Eq)]
pub struct FieldColumns {
    /// the name of the reference sequence column, usually it is called "seqname".\
    /// This corresponds to the first column in a GTF/GFF file.\
    /// This field is required, and the corresponding column should not contain missing values (nulls).\
    /// You can drop all rows with missing values in this column by calling `df.drop_nulls(Some(["seqname"]))`.
    pub seqname: String,
    /// the name of the source column, usually it is called "source".
    /// This is the second column in a GTF/GFF file. It records the source (HAVANA, ENSEMBL, etc) of the feature.
    pub source: Option<String>,
    /// The name of the feature type column whose values should be "gene", "transcript" and "exon", etc.
    /// This should be the third column in a GTF/GFF file.
    /// This field is required for transcriptome related methods, like `introns()`, `get_transcript_sequeneces()`, etc.
    pub feature_type: Option<String>,
    /// The column name of the start position column, usually it is called "start".
    /// This is the fourth column in a GTF/GFF file.
    /// This field is required, and the corresponding column should not contain missing values (nulls).
    pub start: String,
    /// The column name of the end position column, usually it is called "end".
    /// This is the fifth column in a GTF/GFF file.
    /// This field is required, and the corresponding column should not contain missing values (nulls).
    pub end: String,
    /// The column name of the score column, usually it is called "score".
    pub score: Option<String>,
    /// The column name of the strand column, usually it is called "strand".
    /// This is the seventh column in a GTF/GFF file.\
    /// This field is required, and the corresponding column should not contain missing values (nulls) for calling most Grangers methods.\
    /// If this field is missing, you can add one by calling `df.add_column("strand", vec!['.'; df.height()])`.\
    /// If it contains missing values, you can fill them by calling `df.fill_none("strand", '.')`.
    pub strand: String,
    /// The column name of the phase column, usually it is called "frame" or "phase".
    /// This is the eighth column in a GTF/GFF file.
    pub phase: Option<String>,
    /// The column name of the gene ID column, usually it is called "gene_id".
    pub gene_id: Option<String>,
    /// The column name of the transcript ID column, usually it is called "transcript_id".
    pub transcript_id: Option<String>,
    // The column name of the exon ID column, usually it is called "exon_id".
    pub exon_id: Option<String>,
    /// The column name of the exon number (order) column, usually it is called "exon_number".
    /// This column is used to sort the exons of a transcript.
    /// If this column is missing, the exons will be sorted by their start positions.
    pub exon_number: Option<String>,
}

impl FieldColumns {
    /// get a reference to the seqname field
    pub fn seqname(&self) -> &str {
        self.seqname.as_str()
    }
    /// get a reference to the source field
    pub fn source(&self) -> Option<&str> {
        self.source.as_deref()
    }
    /// get a reference to the feature_type field
    pub fn feature_type(&self) -> Option<&str> {
        self.feature_type.as_deref()
    }
    /// get a reference to the start field
    pub fn start(&self) -> &str {
        self.start.as_str()
    }
    /// get a reference to the end field
    pub fn end(&self) -> &str {
        self.end.as_str()
    }
    /// get a reference to the score field
    pub fn score(&self) -> Option<&str> {
        self.score.as_deref()
    }
    /// get a reference to the strand field
    pub fn strand(&self) -> &str {
        self.strand.as_str()
    }
    /// get a reference to the phase field
    pub fn phase(&self) -> Option<&str> {
        self.phase.as_deref()
    }
    /// get a reference to the gene_id field
    pub fn gene_id(&self) -> Option<&str> {
        self.gene_id.as_deref()
    }
    /// get a reference to the transcript_id field
    pub fn transcript_id(&self) -> Option<&str> {
        self.transcript_id.as_deref()
    }
    /// get a reference to the exon_id field
    pub fn exon_id(&self) -> Option<&str> {
        self.exon_id.as_deref()
    }
    /// get a reference to the exon_number field
    pub fn exon_number(&self) -> Option<&str> {
        self.exon_number.as_deref()
    }
}

impl Default for FieldColumns {
    fn default() -> Self {
        Self {
            seqname: "seqname".to_string(),
            source: Some("source".to_string()),
            feature_type: Some("feature_type".to_string()),
            start: "start".to_string(),
            end: "end".to_string(),
            score: Some("score".to_string()),
            strand: "strand".to_string(),
            phase: Some("phase".to_string()),
            gene_id: Some("gene_id".to_string()),
            transcript_id: Some("transcript_id".to_string()),
            exon_id: Some("exon_id".to_string()),
            exon_number: Some("exon_number".to_string()),
        }
    }
}

impl FieldColumns {
    /// create a new TxColumns struct.
    pub fn new<T: Into<String>>(
        seqname: T,
        source: Option<T>,
        feature_type: Option<T>,
        start: T,
        end: T,
        score: Option<T>,
        strand: T,
        phase: Option<T>,
        gene_id: Option<T>,
        transcript_id: Option<T>,
        exon_id: Option<T>,
        exon_number: Option<T>,
    ) -> Self {
        Self {
            seqname: seqname.into(),
            source: source.map(|v| v.into()),
            feature_type: feature_type.map(|v| v.into()),
            start: start.into(),
            end: end.into(),
            score: score.map(|v| v.into()),
            strand: strand.into(),
            phase: phase.map(|v| v.into()),
            gene_id: gene_id.map(|v| v.into()),
            transcript_id: transcript_id.map(|v| v.into()),
            exon_id: exon_id.map(|v| v.into()),
            exon_number: exon_number.map(|v| v.into()),
        }
    }

    /// validate the provided dataframe.
    /// If everything is Ok, return Ok(None)
    /// If some required fields are missing, return an error.
    /// if some optional fields are missing, return a FieldColumns struct with these missing fields set as None.
    /// If not, return a FieldColumns struct with all missing field as None.
    pub fn validate(&self, df: &DataFrame, complain: bool) -> anyhow::Result<Option<Self>> {
        let mut fc = self.clone();
        // check required fields
        if df.column(fc.seqname.as_str()).is_err()
            || df.column(fc.seqname.as_str())?.null_count() > 0
        {
            bail!(
                "The dataframe either does not contain the specified seqname column {} or this column contains missing value; Cannot proceed",
                fc.seqname
            );
        }
        if df.column(fc.start.as_str()).is_err() || df.column(fc.start.as_str())?.null_count() > 0 {
            bail!(
                "The dataframe does not contain the specified start column {} or this column contains missing value; Cannot proceed",
                fc.start
            );
        }
        if df.column(fc.end.as_str()).is_err() || df.column(fc.end.as_str())?.null_count() > 0 {
            bail!(
                "The dataframe does not contain the specified end column {} or this column contains missing value; Cannot proceed",
                fc.end
            );
        }
        if df.column(fc.strand.as_str()).is_err() {
            bail!(
                "The dataframe does not contain the specified strand column {}; Cannot proceed. You can add one by calling `df.add_column(\"strand\", vec!['.'; df.height()])`",
                fc.strand
            );
        }

        // check optional fields
        if let Some(s) = fc.source.to_owned() {
            if df.column(s.as_str()).is_err() {
                fc.source = None;
            }
            if complain {
                warn!("The provided source column {} is not found in the dataframe; It will be ignored", s)
            }
        }
        if let Some(s) = fc.feature_type.to_owned() {
            if df.column(s.as_str()).is_err() {
                fc.feature_type = None;
            }
            if complain {
                warn!("The provided feature_type column {} is not found in the dataframe; It will be ignored", s)
            }
        }
        if let Some(s) = fc.score.to_owned() {
            if df.column(s.as_str()).is_err() {
                fc.score = None;
            }
            if complain {
                warn!("The provided score column {} is not found in the dataframe; It will be ignored", s)
            }
        }
        if let Some(s) = fc.phase.to_owned() {
            if df.column(s.as_str()).is_err() {
                fc.phase = None;
            }
            if complain {
                warn!("The provided phase column {} is not found in the dataframe; It will be ignored", s)
            }
        }
        if let Some(s) = fc.gene_id.to_owned() {
            if df.column(s.as_str()).is_err() {
                fc.gene_id = None;
            }
            if complain {
                warn!("The provided gene_id column {} is not found in the dataframe; It will be ignored", s)
            }
        }
        if let Some(s) = fc.transcript_id.to_owned() {
            if df.column(s.as_str()).is_err() {
                fc.transcript_id = None;
            }
            if complain {
                warn!("The provided transcript_id column {} is not found in the dataframe; It will be ignored", s)
            }
        }
        if let Some(s) = fc.exon_id.to_owned() {
            if df.column(s.as_str()).is_err() {
                fc.exon_id = None;
            }
            if complain {
                warn!("The provided exon_id column {} is not found in the dataframe; It will be ignored", s)
            }
        }
        if let Some(s) = fc.exon_number.to_owned() {
            if df.column(s.as_str()).is_err() {
                fc.exon_number = None;
            }
            if complain {
                warn!("The provided exon_number column {} is not found in the dataframe; It will be ignored", s)
            }
        }
        // if fc is not the same as self, it means that some optional fields are missing.
        if &fc == self {
            Ok(None)
        } else {
            Ok(Some(fc))
        }
    }

    pub fn get_colname<T: AsRef<str>>(&self, field: T, complain: bool) -> Option<&str> {
        match field.as_ref() {
            "seqname" => Some(self.seqname.as_str()),
            "source" => self.source(),
            "feature_type" => self.feature_type(),
            "start" => Some(self.start.as_str()),
            "end" => Some(self.end.as_str()),
            "score" => self.score(),
            "strand" => Some(self.strand.as_str()),
            "phase" => self.phase(),
            "gene_id" => self.gene_id(),
            "transcript_id" => self.transcript_id(),
            "exon_id" => self.exon_id(),
            "exon_number" => self.exon_number(),
            _ => {
                if complain {
                    warn!(
                        "The provided field {} is not a valid field name; It will be ignored",
                        field.as_ref()
                    );
                }
                None
            }
        }
    }
}

#[derive(Clone, PartialEq, Eq)]
pub enum IntronsBy {
    /// find exons of each gene and (deduplicate if needed)
    Gene,
    /// find exons of each transcript (no deduplication)
    Transcript,
    /// find exons according to a custom column
    Other(String),
}

impl AsRef<str> for IntronsBy {
    fn as_ref(&self) -> &str {
        match self {
            IntronsBy::Gene => "gene_id",
            IntronsBy::Transcript => "transcript_id",
            IntronsBy::Other(s) => s.as_str(),
        }
    }
}
