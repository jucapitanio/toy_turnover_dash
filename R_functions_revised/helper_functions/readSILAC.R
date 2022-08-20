#' This is a helper function used by the \code{importIP2} function to read in tab-delimited SILAC proteomic data from the IP2 pipeline and transform the data to save protein IDs on peptide lines. This function returns only raw, unfiltered, peptide-level data (singleton and non-singleton peptides), and the SILAC analysis.
#' 
#' @param filepath path to the TMT data file to read in
#' @return an unfiltered dataset containing only peptide-level protein sequences

.readSILAC <- function(filepath) {
  
  data <- read_csv(filepath)
  
  # wrangle the singleton peptide data to be consistent with the rest of the df
  data <- data %>% 
    rename(prof_score = singleton, uniprot = Proteins, description = `Protein Description`) %>% 
    mutate(singleton = if_else(grepl(pattern = "singleton", x = Sequence),T,F)) %>% 
    mutate(Sequence = gsub("singleton","",Sequence))
  
  return(data)
}