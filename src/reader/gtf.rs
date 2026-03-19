use crate::grangers_utils::{is_gzipped, FileFormat};
use anyhow::{self, Context};
use flate2::bufread::MultiGzDecoder;
use noodles::{gff, gtf};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::{collections::HashMap, path::Path};
use tracing::{info, warn};

#[derive(Copy, Clone)]
/// Represents the modes available for selecting attributes during data processing.
///
/// This enum is used to specify how attributes should be handled, particularly when dealing with
/// genomic data or similar structured datasets. Depending on the mode selected, different sets of
/// attributes may be included in the output or analysis.
///
/// # Variants
///
/// * `Essential` - In this mode, only a core set of essential attributes are included or processed.
///   This is typically used to focus on the most critical data elements and reduce complexity or processing time.
///
/// * `Full` - This mode includes all available attributes for each data entry, ensuring comprehensive
///   coverage and detail. It is used when complete data representation is necessary for the analysis or output.
///
/// # Examples
///
/// Usage of `AttributeMode` can depend on the context, such as configuring data extraction or processing functions:
///
/// ```rust
/// let mode = AttributeMode::Essential;
/// ```
///
/// This might instruct a data processing function to only consider the most important attributes,
/// streamlining the analysis for efficiency or clarity.
pub enum AttributeMode {
    Essential,
    Full,
}

impl AttributeMode {
    /// Creates an `AttributeMode` instance based on a boolean flag.
    ///
    /// This method allows for convenient instantiation of `AttributeMode` based on a simple boolean value,
    /// where `true` maps to `AttributeMode::Full` and `false` maps to `AttributeMode::Essential`.
    ///
    /// # Arguments
    ///
    /// * `is_full`: A boolean value indicating whether the full attribute mode should be used.
    ///
    /// # Returns
    ///
    /// Returns an `AttributeMode` variant corresponding to the boolean input: `Full` if `is_full` is `true`,
    /// otherwise `Essential`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// let mode_full = AttributeMode::from(true);
    /// assert_eq!(mode_full, AttributeMode::Full);
    ///
    /// let mode_essential = AttributeMode::from(false);
    /// assert_eq!(mode_essential, AttributeMode::Essential);
    /// ```
    pub fn from(is_full: bool) -> AttributeMode {
        if is_full {
            AttributeMode::Full
        } else {
            AttributeMode::Essential
        }
    }
    /// Checks if the `AttributeMode` is `Full`.
    ///
    /// This method allows for easy checking of whether an `AttributeMode` instance represents
    /// the full set of attributes (`Full`) or just the essential subset (`Essential`).
    ///
    /// # Returns
    ///
    /// Returns `true` if the mode is `Full`, otherwise `false`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// let mode = AttributeMode::Full;
    /// assert!(mode.is_full());
    ///
    /// let mode = AttributeMode::Essential;
    /// assert!(!mode.is_full());
    /// ```
    pub fn is_full(&self) -> bool {
        match self {
            AttributeMode::Full => true,
            AttributeMode::Essential => false,
        }
    }
}

#[derive(Clone)]
/// Stores attributes related to genomic data processing, categorizing them into essential and extra attributes.
///
/// This struct is particularly useful for managing genomic feature attributes extracted from files,
/// allowing for differentiated handling based on attribute importance and processing mode.
///
/// # Fields
///
/// * `file_type`: An instance of `FileFormat` specifying the format of the source genomic data file.
/// * `essential`: A `HashMap` storing essential attributes. Keys are attribute names and values
///   are vectors of `Option<String>` representing the attribute values for each genomic feature.
/// * `extra`: An optional `HashMap` storing additional, non-essential attributes when operating in full mode.
/// * `tally`: A counter indicating the number of genomic features processed, used to ensure attribute vectors
///   are correctly sized.
///
pub struct Attributes {
    pub file_type: FileFormat,
    pub essential: HashMap<String, Vec<Option<String>>>,
    pub extra: Option<HashMap<String, Vec<Option<String>>>>,
    pub tally: usize,
}

impl Attributes {
    /// Constructs a new `Attributes` instance based on the provided attribute mode and file format.
    ///
    /// ### Arguments
    ///
    /// * `mode`: An instance of `AttributeMode` determining whether to include extra attributes.
    /// * `file_type`: The format of the genomic data file, influencing which attributes are considered essential.
    ///
    /// ### Returns
    ///
    /// Returns an `anyhow::Result<Attributes>` containing the new instance if successful, or an error if creation fails.
    ///
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
    /// Adds a new set of attributes for a genomic feature to the `Attributes` instance.
    ///
    /// ### Arguments
    ///
    /// * `hm`: A mutable reference to a `HashMap<String, String>` containing the attribute names and values to add.
    ///
    /// # Examples
    ///
    /// Creating a new `Attributes` instance in essential mode for a hypothetical file format:
    ///
    /// ```rust
    /// let attr_mode = AttributeMode::Essential;
    /// let file_format = FileFormat::new_custom_format(); // assuming this is a valid method
    /// let attributes = Attributes::new(attr_mode, file_format)?;
    /// ```
    ///
    /// Adding attributes for a genomic feature:
    ///
    /// ```rust
    /// let mut feature_attributes = HashMap::new();
    /// feature_attributes.insert("gene".to_string(), "BRCA1".to_string());
    /// attributes.push(&mut feature_attributes);
    /// ```
    ///
    /// # Note
    ///
    /// When adding attributes with the `push` method, the attributes are categorized into essential and extra based
    /// on the file format's definition of essential attributes and the current mode of the `Attributes` instance.
    /// Extra attributes are only stored if the instance is in full mode.
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
/// Represents the type of a genomic feature.
///
/// This enumeration categorizes different types of genomic features commonly found in bioinformatics analyses,
/// such as genes, transcripts, and exons. It provides a structured way to refer to these different feature types,
/// facilitating data processing and annotation tasks.
///
/// # Variants
///
/// * `Gene` - Represents a gene, a fundamental unit of heredity and a primary sequence element in genomic studies.
/// * `Transcript` - Represents a transcript, which is the RNA copy of a gene used in the process of gene expression.
/// * `Exon` - Represents an exon, a segment of a DNA or RNA molecule containing information coding for a protein or peptide sequence.
/// * `Other` - Represents any other type of genomic feature not covered by the specific categories listed above.
///
pub enum FeatureType {
    Gene,
    Transcript,
    Exon,
    Other,
}

impl std::str::FromStr for FeatureType {
    type Err = anyhow::Error;

    /// Parses a string slice into a `FeatureType`.
    ///
    /// This method provides a mechanism to convert textual feature types into the respective `FeatureType` enumeration variants.
    ///
    /// ### Arguments
    ///
    /// * `s`: A string slice representing the name of a genomic feature type.
    ///
    /// ### Returns
    ///
    /// Returns a `Result<FeatureType, anyhow::Error>`:
    /// * `Ok(FeatureType)` for a recognized feature type or `FeatureType::Other` for unrecognized strings.
    /// * `Err(anyhow::Error)` is theoretically possible but not currently implemented since all errors default to `Other`.
    /// ### Examples
    ///
    /// Converting a string to a `FeatureType`:
    ///
    /// ```rust
    /// use std::str::FromStr;
    /// let gene_type = FeatureType::from_str("gene").unwrap();
    /// assert_eq!(gene_type, FeatureType::Gene);
    ///
    /// let unknown_type = FeatureType::from_str("nonexistent").unwrap();
    /// assert_eq!(unknown_type, FeatureType::Other);
    /// ```
    ///
    /// ### Errors
    ///
    /// The `from_str` method returns an `anyhow::Result<FeatureType>`:
    /// * `Ok(FeatureType)` if the string successfully maps to a `FeatureType`.
    /// * `Err(anyhow::Error)` if there is an unexpected error during parsing, though in current implementation,
    ///   it will always return `Ok(FeatureType)` since unrecognized types default to `FeatureType::Other`.
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
/// Represents a generic structure for genomic features and annotations.
///
/// This struct is used to store information typically found in genomic data formats such as GFF, GTF,
/// or custom annotation files. It provides a comprehensive representation of genomic features, including
/// identifiers, sources, types, positions, scores, strands, phases, and associated attributes.
///
/// `GStruct` can be intialized from a GTF or GFF file, or constructed manually by providing the fields to the `new` method.
///
/// # Fields
///
/// * `seqid`: A vector of `String` representing sequence identifiers, such as chromosome names or contig IDs.
/// * `source`: A vector of `String` indicating the sources of the genomic features, such as the database or
///   algorithm that generated the annotation.
/// * `feature_type`: A vector of `String` describing the types of genomic features, such as 'gene', 'exon', or 'CDS'.
/// * `start`: A vector of `i64` indicating the start positions of the genomic features.
/// * `end`: A vector of `i64` indicating the end positions of the genomic features.
/// * `score`: A vector of `Option<f32>` representing the scores associated with the genomic features, which can
///   be null if no score is provided.
/// * `strand`: A vector of `Option<String>` indicating the strands of the genomic features, typically '+' or '-',
///   but can be null if the strand is not specified.
/// * `phase`: A vector of `Option<String>` representing the phase of the genomic features, important for features
///   like CDS; can be null if not applicable.
/// * `attributes`: An `Attributes` instance storing additional information associated with each genomic feature,
///   structured as essential and extra attributes based on the data processing mode.
/// * `misc`: An optional `HashMap<String, Vec<String>>` for storing miscellaneous information that does not fit
///   into the structured fields above.
///
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
    /// Constructs a `GStruct` instance from a GTF (Gene Transfer Format) file.
    ///
    /// This function reads genomic feature information from a GTF file and initializes a `GStruct`
    /// instance with the data extracted. It supports both plain text and gzipped GTF files, automatically
    /// detecting the file format. Based on the specified `AttributeMode`, it categorizes attributes
    /// into essential and extra.
    ///
    /// # Type Parameters
    ///
    /// * `T`: A type that can be referenced as a file path, implementing the `AsRef<Path>` trait.
    ///
    /// # Arguments
    ///
    /// * `file_path`: The file path to the GTF file to be read. Can be either plain text or gzipped.
    /// * `am`: The `AttributeMode` determining how to handle additional attributes found within the GTF file.
    ///
    /// # Returns
    ///
    /// Returns `anyhow::Result<GStruct>`:
    /// * `Ok(GStruct)`: A `GStruct` instance populated with data from the GTF file if successful.
    /// * `Err(anyhow::Error)`: An error if there is a problem opening the file, reading from it, or parsing its content.
    ///
    /// # Examples
    ///
    /// Reading genomic features from a GTF file and creating a `GStruct`:
    ///
    /// ```rust
    /// use std::path::Path;
    ///
    /// let gstruct = GStruct::from_gtf(Path::new("path/to/data.gtf"), AttributeMode::Essential)?;
    /// ```
    ///
    /// # Note
    ///
    /// The function also inserts a 'file_type' attribute into the `misc` field of the resulting `GStruct`,
    /// indicating that the data was sourced from a GTF file. This can be useful for downstream processing or
    /// metadata tracking.
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
            let mut rdr = gtf::io::Reader::new(BufReader::new(MultiGzDecoder::new(inner_rdr)));
            gr._from_gtf(&mut rdr)?;
        } else {
            let mut rdr = gtf::io::Reader::new(inner_rdr);
            gr._from_gtf(&mut rdr)?;
        }

        Ok(gr)
    }

    fn _from_gtf<T: BufRead>(&mut self, rdr: &mut gtf::io::Reader<T>) -> anyhow::Result<()> {
        // initiate a reusable hashmap to take the attributes of each record
        let mut rec_attr_hm: HashMap<String, String> = HashMap::with_capacity(100);
        let mut n_comments = 0usize;
        let mut n_records = 0usize;

        // parse the file
        for l in rdr.lines() {
            let line = l?;
            match line.kind()    {
                gtf::line::Kind::Record => {
                    let r = line.as_record().with_context(|| format!("Failed parsing a record line: {:#?}", line.as_record()))?
                    .with_context(|| format!("Failed parsing a record line: {:#?}", line.as_record()))?;
                    n_records += 1;
                    // parse essential fields

                    GStruct::push(&mut self.seqid, r.reference_sequence_name().to_string());
                    GStruct::push(&mut self.source, r.source().to_string());
                    GStruct::push(&mut self.feature_type, r.ty().to_string());
                    GStruct::push(&mut self.start, r.start()?.get() as i64);
                    GStruct::push(&mut self.end, r.end()?.get() as i64);
                    GStruct::push(&mut self.score, r.score().transpose()?);
                    let strand_opt = r.strand().map(|st| {
                        let st_debug = format!("{:?}", st);
                        match st_debug.as_str() {
                            "Reverse" => String::from("-"),
                            _ => String::from("+"),
                        }
                    })?;
                    GStruct::push(&mut self.strand, Some(strand_opt));

                    if let Some(ph_result) = r.phase() {
                        let ph = ph_result?;
                        let ph_debug = format!("{:?}", ph);
                        let phase_str = match ph_debug.as_str() {
                            "Zero" => Some(String::from("0")),
                            "One" => Some(String::from("1")),
                            "Two" => Some(String::from("2")),
                            _ => Some(String::from("0")),
                        };
                        GStruct::push(&mut self.phase, phase_str);
                    } else {
                        GStruct::push(&mut self.phase, None);
                    }

                    // parse attributes
                    rec_attr_hm.clear();
                    // In noodles-gtf, attributes().iter() returns an iterator of (key, value) tuples
                    let attributes = r.attributes()?;
                    for att in gff::feature::record::Attributes::iter(&attributes) {
                        let (key, value) = att?;
                        match value {
                            gff::feature::record::attributes::field::Value::String(s) => {
                                rec_attr_hm.insert(key.to_string(), s.to_string());
                            }
                            gff::feature::record::attributes::field::Value::Array(a) => {
                                rec_attr_hm.insert(
                                    key.to_string(), 
                                    a.iter()
                                        .map(|result| result.map(|s| s.to_string()))
                                        .collect::<Result<Vec<_>, _>>()?
                                        .join(",")
                                );
                            }
                        }
                    }
                    self.attributes.push(&mut rec_attr_hm);
                }
                gtf::line::Kind::Comment => {
                    if let Some(c) = line.as_comment() {
                        n_comments += 1;
                        if let Some(misc) = self.misc.as_mut() {
                            misc.entry(String::from("comments"))
                                .and_modify(|v| v.push(c.to_string()))
                                .or_insert(vec![c.to_string()]);
                        }
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
    /// Constructs a `GStruct` instance from a GFF (Generic Feature Format) file.
    ///
    /// This function reads genomic feature information from a GFF file and populates a `GStruct`
    /// instance with the data extracted. It supports both plain text and gzipped GFF files, automatically
    /// detecting and handling the file format accordingly. Attributes within the file are handled based on
    /// the specified `AttributeMode`, categorizing them into essential and additional attributes.
    ///
    /// # Type Parameters
    ///
    /// * `T`: A type that can be referenced as a file path, implementing the `AsRef<Path>` trait.
    ///
    /// # Arguments
    ///
    /// * `file_path`: The path to the GFF file to be read. The file can be in plain text or gzipped format.
    /// * `am`: The `AttributeMode` determining how additional attributes found within the GFF file should be handled.
    ///
    /// # Returns
    ///
    /// Returns `anyhow::Result<GStruct>`:
    /// * `Ok(GStruct)`: A `GStruct` instance populated with data from the GFF file if successful.
    /// * `Err(anyhow::Error)`: An error if there is an issue with opening the file, reading from it, or parsing its content.
    ///
    /// # Examples
    ///
    /// Reading genomic features from a GFF file and initializing a `GStruct`:
    ///
    /// ```rust
    /// use std::path::Path;
    ///
    /// let gstruct = GStruct::from_gff(Path::new("path/to/data.gff"), AttributeMode::Essential)?;
    /// ```
    ///
    /// # Note
    ///
    /// The function initializes the `GStruct` instance by setting up appropriate structures to store
    /// essential and, depending on the `AttributeMode`, extra attributes. It ensures the data from the
    /// GFF file is correctly interpreted and stored for downstream analysis or processing.
    pub fn from_gff<T: AsRef<Path>>(file_path: T, am: AttributeMode) -> anyhow::Result<GStruct> {
        let mut gr = GStruct::new(am, FileFormat::GFF)?;

        let file = File::open(file_path)?;
        let mut inner_rdr = BufReader::new(file);
        // instantiate the struct
        if is_gzipped(&mut inner_rdr)? {
            info!("auto-detected gzipped file - reading via decompression");
            let mut rdr = gff::io::Reader::new(BufReader::new(MultiGzDecoder::new(inner_rdr)));
            gr._from_gff(&mut rdr)?;
        } else {
            let mut rdr = gff::io::Reader::new(inner_rdr);
            gr._from_gff(&mut rdr)?;
        }
        Ok(gr)
    }

    fn _from_gff<T: BufRead>(&mut self, rdr: &mut gff::io::Reader<T>) -> anyhow::Result<()> {
        // initiate a reusable hashmap to take the attributes of each record
        let mut rec_attr_hm: HashMap<String, String> = HashMap::with_capacity(100);
        let mut n_comments = 0usize;
        let mut n_records = 0usize;
        let mut n_strand_none = 0usize;
        let mut n_strand_unknown = 0usize;

        // parse the file
        for l in rdr.lines() {
            let line = l?;
            match line.kind() {
                gff::line::Kind::Record => {
                    let r = line
                        .as_record()
                        .with_context(|| format!("Failed parseing a record line: {:#?}", line))??;
                    n_records += 1;
                    GStruct::push(&mut self.seqid, r.reference_sequence_name().to_string());
                    GStruct::push(&mut self.source, r.source().to_string());
                    GStruct::push(&mut self.feature_type, r.ty().to_string());
                    GStruct::push(&mut self.start, r.start()?.get() as i64);
                    GStruct::push(&mut self.end, r.end()?.get() as i64);

                    if let Some(s) = r.score() {
                        GStruct::push(&mut self.score, Some(s?));
                    } else {
                        GStruct::push(&mut self.score, None);
                    }

                    let strand_str = {
                        let strand_val = r.strand()?;
                        let strand_debug = format!("{:?}", strand_val);
                        match strand_debug.as_str() {
                            "Reverse" => Some(String::from("-")),
                            "Forward" => Some(String::from("+")),
                            "None" => {
                                n_strand_none += 1;
                                Some(String::from("+"))
                            }
                            "Unknown" => {
                                n_strand_unknown += 1;
                                Some(String::from("+"))
                            }
                            _ => Some(String::from("+")),
                        }
                    };
                    GStruct::push(&mut self.strand, strand_str);

                    if let Some(p) = r.phase() {
                        let phase_val = p?;
                        let phase_debug = format!("{:?}", phase_val);
                        let phase_result = match phase_debug.as_str() {
                            "Zero" => Some(String::from("0")),
                            "One" => Some(String::from("1")),
                            "Two" => Some(String::from("2")),
                            _ => None,
                        };
                        GStruct::push(&mut self.phase, phase_result);
                    } else {
                        GStruct::push(&mut self.phase, None);
                    }

                    // parse attributes
                    rec_attr_hm.clear();
                    for attr in r.attributes().iter() {
                        let (attrk, attrv) = attr?;

                        match attrv {
                            gff::record::attributes::field::Value::String(val) => {
                                rec_attr_hm.insert(attrk.to_string(), val.clone().to_string());
                            }
                            gff::record::attributes::field::Value::Array(a) => {
                                let mut arr = Vec::new();
                                for s in a.iter() {
                                    arr.push(s.to_string());
                                }

                                rec_attr_hm.insert(attrk.to_string(), arr.join(","));

                                // anyhow::bail!("Currently, having multiple values associated with a single GFF attributed is not supported.");
                            }
                        }
                    }
                    self.attributes.push(&mut rec_attr_hm);
                }
                gff::line::Kind::Comment => {
                    let c = line
                        .as_comment()
                        .with_context(|| format!("failed parsing a comment line: {:#?}", line))?;
                    n_comments += 1;
                    if let Some(misc) = self.misc.as_mut() {
                        misc.entry(String::from("comments"))
                            .and_modify(|v| v.push(c.to_string()))
                            .or_insert(vec![c.to_string()]);
                    }
                    continue;
                }
                gff::line::Kind::Directive => {
                    let d = line
                        .as_directive()
                        .with_context(|| format!("failed parsing a directive line: {:#?}", line))?;
                    // we create a string containing the key and value fields separated by space for the directive
                    let dstring = format!("{} {}", d.key(), d.value().map(|v| v.to_string()).unwrap_or_else(|| String::from("")));

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

        if n_strand_none > 0 {
            warn!(
                "{} records have no strand information, set to '+'",
                n_strand_none
            );
        }

        if n_strand_unknown > 0 {
            warn!(
                "{} records have unknown strand information, set to '+'",
                n_strand_unknown
            );
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
    /// Constructs a new instance of `GStruct` for storing and managing genomic feature information.
    ///
    /// This method initializes `GStruct` with pre-allocated storage for various genomic feature
    /// attributes and sets up the attributes according to the specified attribute handling mode and file type.
    ///
    /// # Arguments
    ///
    /// * `attribute_mode`: An `AttributeMode` determining how additional attributes found within the data should be handled.
    /// * `file_type`: The `FileFormat` indicating the format of the genomic data file being processed.
    ///
    /// # Returns
    ///
    /// Returns `anyhow::Result<GStruct>`:
    /// * `Ok(GStruct)`: If the `GStruct` instance was successfully created.
    /// * `Err(anyhow::Error)`: If there was an error initializing the `Attributes`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// let gstruct = GStruct::new(AttributeMode::Essential, FileFormat::GFF)?;
    /// ```
    ///
    /// This function is primarily used during the initial setup phase of genomic data processing pipelines.
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

    /// Adds a value to the end of a vector.
    ///
    /// A generic method provided to push a debuggable and cloneable value onto a vector.
    ///
    /// # Type Parameters
    ///
    /// * `T`: The type of the elements in the vector. Must implement `std::fmt::Debug` and `Clone` traits.
    ///
    /// # Arguments
    ///
    /// * `vec`: A mutable reference to a vector of type `T`.
    /// * `val`: The value to be added to the vector.
    ///
    /// # Examples
    ///
    /// ```rust
    /// let mut vec = Vec::new();
    /// GStruct::push(&mut vec, "Example");
    /// assert_eq!(vec, ["Example"]);
    /// ```
    ///
    /// This method is used to add individual elements to the various attribute vectors of a `GStruct` instance.
    // TODO: might need a better generic type
    fn push<T: std::fmt::Debug + Clone>(vec: &mut Vec<T>, val: T) {
        vec.push(val);
    }

    /// Appends elements from one vector to another.
    ///
    /// A generic method for appending all elements from one vector to another.
    ///
    /// # Type Parameters
    ///
    /// * `T`: The type of the elements in the vectors. Must implement `ToString` and `Clone` traits.
    ///
    /// # Arguments
    ///
    /// * `vec`: A mutable reference to the main vector to which elements will be appended.
    /// * `patch`: A mutable reference to the vector containing elements to append to the main vector.
    ///
    /// # Examples
    ///
    /// ```rust
    /// let mut main_vec = vec!["First".to_string()];
    /// let mut patch_vec = vec!["Second".to_string(), "Third".to_string()];
    /// GStruct::append(&mut main_vec, &mut patch_vec);
    /// assert_eq!(main_vec, ["First", "Second", "Third"]);
    /// ```
    ///
    /// This method is used for combining vectors, typically used for merging data from different sources
    /// or batches into a single `GStruct` instance's attribute vectors.
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
        let mut rdr = gtf::io::Reader::new(GTF_RECORD);
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
                        assert!(file_type == FileFormat::GTF);
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
        let mut rdr = gff::io::Reader::new(GFF_RECORD);
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
                        assert!(file_type == FileFormat::GFF);
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
