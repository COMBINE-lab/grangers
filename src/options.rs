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
/// - ignore_strand: whether to ignore the strand information when merging.
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
        let mut by_hash: HashSet<String> = by.iter().map(|n| n.as_ref().to_string()).collect();

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
#[derive(Clone, PartialEq, Eq, Debug)]
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
    /// The column name of the gene_id column, usually it is called "gene_id".
    pub gene_id: Option<String>,
    /// The column name of the gene_name column, usually it is called "gene_id".
    pub gene_name: Option<String>,
    /// The column name of the transcript ID column, usually it is called "transcript_id".
    pub transcript_id: Option<String>,
    // The column name of the exon ID column, usually it is called "exon_id".
    // pub exon_id: Option<String>,
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

    /// get a reference to the gene_name field
    pub fn gene_name(&self) -> Option<&str> {
        self.gene_name.as_deref()
    }

    /// get a reference to the transcript_id field
    pub fn transcript_id(&self) -> Option<&str> {
        self.transcript_id.as_deref()
    }
    // get a reference to the exon_id field
    // pub fn exon_id(&self) -> Option<&str> {
    //     self.exon_id.as_deref()
    // }
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
            gene_name: Some("gene_name".to_string()),
            transcript_id: Some("transcript_id".to_string()),
            exon_number: Some("exon_number".to_string()),
        }
    }
}

impl FieldColumns {
    pub fn optional_fields(&self) -> [Option<&str>; 8] {
        [
            self.source(),
            self.feature_type(),
            self.score(),
            self.phase(),
            self.gene_id(),
            self.gene_name(),
            self.transcript_id(),
            self.exon_number(),
        ]
    }
    pub fn gtf_fields(&self) -> [&str; 8] {
        [
            self.seqname(),
            self.source().unwrap_or(""),
            self.feature_type().unwrap_or(""),
            self.start(),
            self.end(),
            self.score().unwrap_or(""),
            self.strand(),
            self.phase().unwrap_or(""),
        ]
    }

    
    pub fn gtf_attributes(&self) -> [Option<&str>; 4] {
        [
            self.gene_id(),
            self.gene_name(),
            self.transcript_id(),
            self.exon_number(),
        ]
    }

    

    pub fn essential_fields(&self) -> [&str; 4] {
        [self.seqname(), self.start(), self.end(), self.strand()]
    }
    /// Validate itself according to a provided dataframe.
    /// If everything is Ok, return true, else, return false.
    /// If fix is true, try to fix the fields by finding and using the columns with default field names(seqname, start, gene_id, etc) in the dataframe. if the essentail fields are invalid and cannot be fixed, return error
    pub fn is_valid(&self, df: &DataFrame, is_warn: bool, is_bail: bool) -> anyhow::Result<bool> {
        let mut is_valid = true;
        // check required fields
        if df.column(self.seqname()).is_err() {
            is_valid = false;
            if is_warn {
                warn!(
                    "The dataframe does not contain the specified seqname column {}; Cannot proceed. You can add one by calling `df.add_column(\"seqname\", vec!['.'; df.height()])`",
                    self.seqname()
                )
            }
        }
        if df.column(self.start()).is_err() {
            is_valid = false;
            if is_warn {
                warn!(
                    "The dataframe does not contain the specified start column {}; Cannot proceed. You can add one by calling `df.add_column(\"start\", vec!['.'; df.height()])`",
                    self.start()
                )
            }
        }
        if df.column(self.end()).is_err() {
            is_valid = false;
            if is_warn {
                warn!(
                    "The dataframe does not contain the specified end column {}; Cannot proceed. You can add one by calling `df.add_column(\"end\", vec!['.'; df.height()])`",
                    self.end()
                )
            }
        }
        if df.column(self.strand()).is_err() {
            is_valid = false;
            if is_warn {
                warn!(
                    "The dataframe does not contain the specified strand column {}; Cannot proceed. You can add one by calling `df.add_column(\"strand\", vec!['.'; df.height()])`",
                    self.strand()
                )
            }
        }
        // check additional fields
        if let Some(s) = self.source() {
            if df.column(s).is_err() {
                is_valid = false;
                if is_warn {
                    warn!("The provided source column {} is not found in the dataframe; It will be ignored", s)
                }
            }
        }
        if let Some(s) = self.feature_type() {
            if df.column(s).is_err() {
                is_valid = false;
                if is_warn {
                    warn!("The provided feature_type column {} is not found in the dataframe; It will be ignored", s)
                }
            }
        }
        if let Some(s) = self.score() {
            if df.column(s).is_err() {
                is_valid = false;
                if is_warn {
                    warn!("The provided score column {} is not found in the dataframe; It will be ignored", s)
                }
            }
        }
        if let Some(s) = self.phase() {
            if df.column(s).is_err() {
                is_valid = false;
                if is_warn {
                    warn!("The provided phase column {} is not found in the dataframe; It will be ignored", s)
                }
            }
        }

        if let Some(s) = self.gene_id() {
            if df.column(s).is_err() {
                is_valid = false;
                if is_warn {
                    warn!("The provided gene_id column {} is not found in the dataframe; It will be ignored", s)
                }
            }
        }

        if let Some(s) = self.gene_name() {
            if df.column(s).is_err() {
                is_valid = false;
                if is_warn {
                    warn!("The provided gene_name column {} is not found in the dataframe; It will be ignored", s)
                }
            }
        }

        if let Some(s) = self.transcript_id() {
            if df.column(s).is_err() {
                is_valid = false;
                if is_warn {
                    warn!("The provided transcript_id column {} is not found in the dataframe; It will be ignored", s)
                }
            }
        }

        if let Some(s) = self.exon_number() {
            if df.column(s).is_err() {
                is_valid = false;
                if is_warn {
                    warn!("The provided exon_number column {} is not found in the dataframe; It will be ignored", s)
                }
            }
        }

        if !is_valid & is_bail {
            bail!(
                "The FieldColumns is not valid; Please try fix it by calling FieldColumns::fix()."
            )
        }

        if !is_valid & is_warn {
            warn!(
                "The FieldColumns is not valid; Please try fix it by calling FieldColumns::fix()."
            )
        }

        Ok(is_valid)
    }

    pub fn fix(&mut self, df: &DataFrame, is_warn: bool) -> anyhow::Result<()> {
        // try fix required fields
        if df.column(self.seqname()).is_err() {
            if is_warn {
                warn!(
                    "cannot find the specified seqname column {} in the dataframe; try to fix",
                    self.seqname()
                );
            }
            if df.column("seqname").is_ok() {
                self.seqname = "seqname".to_string();
            } else {
                bail!("The dataframe does not contain the specified seqname column {} or a column named \"seqname\"; Cannot fix.", self.seqname());
            }
        }
        if df.column(self.start()).is_err() {
            if is_warn {
                warn!(
                    "cannot find the specified start column {} in the dataframe; try to fix",
                    self.start()
                );
            }
            if df.column("start").is_ok() {
                self.start = "start".to_string();
            } else {
                bail!("The dataframe does not contain the specified start column {} or a column named \"start\"; Cannot fix.", self.start());
            }
        }
        if df.column(self.end()).is_err() {
            if is_warn {
                warn!(
                    "cannot find the specified end column {} in the dataframe; try to fix",
                    self.end()
                );
            }
            if df.column("end").is_ok() {
                self.end = "end".to_string();
            } else {
                bail!("The dataframe does not contain the specified end column {} or a column named \"end\"; Cannot fix.", self.end());
            }
        }
        if df.column(self.strand()).is_err() {
            if is_warn {
                warn!(
                    "cannot find the specified strand column {} in the dataframe; try to fix",
                    self.strand()
                );
            }
            if df.column("strand").is_ok() {
                self.strand = "strand".to_string();
            } else {
                bail!("The dataframe does not contain the specified strand column {} or a column named \"strand\"; Cannot fix. If this is desired, you can add a dummy strand column by calling `df.add_column(\"strand\", vec!['.'; df.height()])`", self.strand());
            }
        }

        // try fix optional fields
        if let Some(s) = self.source() {
            if df.column(s).is_err() {
                if is_warn {
                    warn!(
                        "cannot find the specified source column {} in the dataframe; try to fix",
                        s
                    );
                }
                self.source = if df.column("source").is_ok() {
                    Some("source".to_string())
                } else {
                    None
                }
            }
        }
        if let Some(s) = self.feature_type() {
            if df.column(s).is_err() {
                if is_warn {
                    warn!("cannot find the specified feature_type column {} in the dataframe; try to fix", s);
                }
                self.feature_type = if df.column("feature_type").is_ok() {
                    Some("feature_type".to_string())
                } else {
                    None
                }
            }
        }
        if let Some(s) = self.score() {
            if df.column(s).is_err() {
                if is_warn {
                    warn!(
                        "cannot find the specified score column {} in the dataframe; try to fix",
                        s
                    );
                }
                self.score = if df.column("score").is_ok() {
                    Some("score".to_string())
                } else {
                    None
                }
            }
        }
        if let Some(s) = self.phase() {
            if df.column(s).is_err() {
                if is_warn {
                    warn!(
                        "cannot find the specified phase column {} in the dataframe; try to fix",
                        s
                    );
                }
                self.phase = if df.column("phase").is_ok() {
                    Some("phase".to_string())
                } else {
                    None
                }
            }
        }
        if let Some(s) = self.gene_id() {
            if df.column(s).is_err() {
                if is_warn {
                    warn!(
                        "cannot find the specified gene_id column {} in the dataframe; try to fix",
                        s
                    );
                }
                self.gene_id = if df.column("gene_id").is_ok() {
                    Some("gene_id".to_string())
                } else {
                    None
                }
            }
        }

        if let Some(s) = self.gene_name() {
            if df.column(s).is_err() {
                if is_warn {
                    warn!(
                        "cannot find the specified gene_name column {} in the dataframe; try to fix",
                        s
                    );
                }
                self.gene_name = if df.column("gene_name").is_ok() {
                    Some("gene_name".to_string())
                } else {
                    None
                }
            }
        }

        if let Some(s) = self.transcript_id() {
            if df.column(s).is_err() {
                if is_warn {
                    warn!("cannot find the specified transcript_id column {} in the dataframe; try to fix", s);
                }
                self.transcript_id = if df.column("transcript_id").is_ok() {
                    Some("transcript_id".to_string())
                } else {
                    None
                }
            }
        }

        if let Some(s) = self.exon_number() {
            if df.column(s).is_err() {
                if is_warn {
                    warn!("cannot find the specified exon_number column {} in the dataframe; try to fix", s);
                }
                self.exon_number = if df.column("exon_number").is_ok() {
                    Some("exon_number".to_string())
                } else {
                    None
                }
            }
        }

        Ok(())
    }

    pub fn update<T: AsRef<str>>(&mut self, field: T, value: T) -> anyhow::Result<()> {
        let value = value.as_ref().to_string();
        match field.as_ref() {
            "seqname" => self.seqname = value,
            "source" => self.source = Some(value),
            "feature_type" => self.feature_type = Some(value),
            "start" => self.start = value,
            "end" => self.end = value,
            "score" => self.score = Some(value),
            "strand" => self.strand = value,
            "phase" => self.phase = Some(value),
            "gene_id" => self.gene_id = Some(value),
            "gene_name" => self.gene_name = Some(value),
            "transcript_id" => self.transcript_id = Some(value),
            "exon_number" => self.exon_number = Some(value),
            _ => bail!("invalid field name: {}", field.as_ref()),
        }

        Ok(())
    }

    pub fn field<T: AsRef<str>>(&self, field: T) -> Option<&str> {
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
            "gene_name" => self.gene_name(),
            "transcript_id" => self.transcript_id(),
            "exon_number" => self.exon_number(),
            _ => None,
        }
    }
    pub fn field_checked<T: AsRef<str>>(
        &self,
        field: T,
        is_bail: bool,
    ) -> anyhow::Result<Option<&str>> {
        match field.as_ref() {
            "seqname" => Ok(Some(self.seqname.as_str())),
            "source" => Ok(self.source()),
            "feature_type" => Ok(self.feature_type()),
            "start" => Ok(Some(self.start.as_str())),
            "end" => Ok(Some(self.end.as_str())),
            "score" => Ok(self.score()),
            "strand" => Ok(Some(self.strand.as_str())),
            "phase" => Ok(self.phase()),
            "gene_id" => Ok(self.gene_id()),
            "gene_name" => Ok(self.gene_name()),
            "transcript_id" => Ok(self.transcript_id()),
            "exon_number" => Ok(self.exon_number()),
            _ => {
                if is_bail {
                    bail!(
                        "The provided field {} is not a valid field name; Cannot proceed",
                        field.as_ref()
                    );
                }

                warn!(
                    "The provided field {} is not a valid field name; It will be ignored",
                    field.as_ref()
                );
                Ok(None)
            }
        }
    }
}
