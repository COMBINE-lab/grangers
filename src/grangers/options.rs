use std::collections::HashSet;
use tracing::{warn};
use anyhow::{bail};

#[derive(Clone)]
pub struct IntronsOptions {
    /// The representitive name of exon feature in the type_col
    pub exon_name: String,
    /// The column name of the feature_type column.
    pub type_col: String,
}

impl Default for IntronsOptions {
    fn default() -> Self {
        Self {
            exon_name: "exon".to_string(),
            type_col: "feature_type".to_string(),
        }
    }
}

impl IntronsOptions {
    pub fn new<T: ToString>(exon_name: T, type_col: T) -> Self {
        Self {
            exon_name: exon_name.to_string(),
            type_col: type_col.to_string(),
        }
    }
}


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
            by: vec![String::from("seqnames"), String::from("strand")],
            slack: 1,
            ignore_strand: false,
        }
    }
}

impl MergeOptions {
    pub fn new<T: ToString>(
        by: Vec<T>,
        ignore_strand: bool,
        slack: i64,
    ) -> anyhow::Result<MergeOptions> {
        // avoid duplicated columns
        let mut by_hash: HashSet<String> = by.into_iter().map(|n| n.to_string()).collect();

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
        if by_hash.insert(String::from("seqnames")) {
            warn!("Added `seqnames` to the `by` vector as it is required.")
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

pub struct GetSequenceOptions {

}

/// Options used for dealing with out-of-boundary features
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum OOBOption {
    Truncate,
    Skip,
}

