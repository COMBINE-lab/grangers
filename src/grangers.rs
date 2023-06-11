mod grangers;
mod grangers_utils;
pub mod reader;
pub use grangers::*;
pub use grangers_utils::FileFormat;
pub use reader::GStruct;

// This function traverses the GTF file:
// - If everything is good, count the # of transcripts for each gene and number of exons for each txp to instantiate the struct
// - If there is any invalid record, report an error and write a clean_gtf.gtf file.
//     We care only about records with type gene, transcript and exon.
//     1. Each transcript/exon record must have a valid `transcript_id` attribute.
//         - If this is not satisfied, it returns an error.
//         - Only the records with a valid transcript_id will be written to the clean_gtf.gtf.
//     2. For gene_id and gene_name attributes,
//         - If these two attributes do not exist in any record, An error will be returned.
//             At the same time, in the clean_gtf.gtf, the two fields will be imputed using the transcript_id fields.
//             However, this means each transcript will have a count because there is no transcript to gene mapping.
//         - If one of these two fields is completely missing, a warning will be generated, and the missing field will be imputed using the other.
//         - if some records miss gene_id and/or gene_name, a warning will be printed,
//             and the missing values will be imputed by the following rules:
//                 For records miss gene_id or gene_name, impute the missing one using the other one;
//                 If both are missing, impute them using transcript_id, which cannot be missing.
//    3. If there is no "transcript" or "gene" feature record, a warning will be printed.
//         Moreover, those missing records will be imputed using the "exon" feature records:
//             The Start and End site of the gene/transcript will be imputed as the bounds of their corresponding exons.
//    4. If the boundaries defined in the transcripts'/genes' feature records do not match those implied by their exons' feature records,
//         report a warning but still use transcripts'/genes' feature records to extract unspliced sequences.
//         To be specific, if some but not all transcripts/genes have their corresponding transcripts'/genes' feature records,
//         or the Start and/or End site defined in the transcript/gene feature records do not match the corresponding exons' bounds,
//         then the existing transcripts'/genes' feature records will be used to extract unspliced transcripts.
//         At the same time, in the clean_gtf.gtf, all genes/transcripts that appeared in the exon feature records
//         will have their corresponding transcripts'/genes' feature records, in which the boundaries match the corresponding exons' bounds.
