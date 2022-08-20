#' This is a helper function used by the \code{importIP2} function to read in tab-delimited TMT proteomic data from the IP2 pipeline and transform the data to save protein IDs on peptide lines. This function returns only raw, unfiltered, peptide-level data, and the associated spectrum counts across each time point.
#' 
#' @param filepath path to the TMT data file to read in
#' @return an unfiltered dataset containing only peptide-level protein sequences

.readTMT <- function(filepath) {
  
  if (endsWith(filepath,"xls")) {
    data <- read_delim(filepath,"\t", escape_double = FALSE, trim_ws = TRUE)
  } else if (endsWith(filepath, "csv")) {
    data <- read_delim(filepath,",")
  } else {
    print("urecognized file format")
  }
  
  # throw out IP2's auto-normalized or empty columns, separate protein column into uniprot id and description.
  data <- data %>%
    select(-starts_with("norm")) %>% 
    select_if(function(x) any(!is.na(x))) %>% 
    separate(col = protein, into = c("uniprot", "description"), sep = " ",
             remove = F, extra = "merge") %>% 
    rename(Unique = unique, Sequence = sequence) %>% 
    #mutate(Unique = if_else(Unique == "U",T,F))
    mutate(Unique = if_else(is.na(Unique),F,T))
  
  return(data)
}