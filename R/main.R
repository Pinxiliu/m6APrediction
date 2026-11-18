# Functions reused from Practical 5. Add documentation in your new package.
#' Convert DNA sequences into categorical nucleotide features
#'
#' Converts each DNA sequence into a data frame of nucleotide positions.
#'
#' @param dna_strings A character vector containing DNA sequences of equal length.
#' @return A data frame containing categorical nucleotide features.
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A sites for multiple DNA sequences
#'
#' Predicts m6A probabilities and classification results for multiple input sequences.
#'
#' @param ml_fit A trained random forest model.
#' @param feature_df A data frame containing the required sequence features.
#' @param positive_threshold Probability threshold for classifying positive m6A sites.
#' @return A data frame containing the original features and prediction results.
#' @import randomForest
#' @importFrom stats predict
#' @export
#' @examples
#' model_path <- system.file("extdata", "rf_fit.rds", package = "m6APrediction")
#' data_path <- system.file("extdata", "m6A_input_example.csv", package = "m6APrediction")
#' ml_fit <- readRDS(model_path)
#' feature_df <- read.csv(data_path)
#' prediction_multiple(ml_fit, head(feature_df, 3))
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}

#' Predict the m6A status of a single DNA sequence
#'
#' Predicts the m6A probability and classification result for one input sequence.
#'
#' @param ml_fit A trained random forest model.
#' @param gc_content GC content of the sequence.
#' @param RNA_type RNA transcript type.
#' @param RNA_region Region of the RNA transcript.
#' @param exon_length Length of the exon.
#' @param distance_to_junction Distance to the nearest exon junction.
#' @param evolutionary_conservation Evolutionary conservation score.
#' @param DNA_5mer Five-nucleotide DNA sequence.
#' @param positive_threshold Probability threshold for classifying positive m6A sites.
#' @return A named character vector containing the predicted probability and m6A status.
#' @export
#' @examples
#' model_path <- system.file("extdata", "rf_fit.rds", package = "m6APrediction")
#' ml_fit <- readRDS(model_path)
#' prediction_single(
#'   ml_fit = ml_fit,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "GGACA",
#'   positive_threshold = 0.5
#' )
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}
