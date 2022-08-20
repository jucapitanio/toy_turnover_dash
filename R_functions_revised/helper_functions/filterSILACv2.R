#' This is a helper function for \code{filterIP2} which will filter SILAC data and throw out some data points.
#'
#' This function does a few things to the data:
#'
#' - throw out all peptides with REGRESSION_FACTOR equal to -1
#'
#' - throw out singleton peptide data with SINGLETON_SCORE below 0.8, UNLESS there is no regular peptide data for that sequence. Any remaining duplicate singleton/regular peptides are averaged into a single entry (summarized).
#'
#' - throw out rows containing NA values
#'
#' - discards columns collected by IP2 irrelevant for normalization, such as REGRESSION_FACTOR and SINGLETON_SCORE, REV_SLOPE_RATIO, and PROBABILITY_SCORE
#'
#' @param data a dataframe containing SILAC data from the IP2 pipeline
#' @return a dataframe of TMT data with some datapoints thrown out

.filterSILAC <- function(data) {
  # Remove poor datapoints:
  data <- data %>%
    filter(!(Ratio == 0 & det_factor == 0)) %>% 
    filter(prof_score >= 0.8)
  
  groupings <- c("Unique", "Sequence")
  if (!is.null(data$Localization)) {
    groupings <- c(groupings, "Localization")
  }
  
  # Separate data and annotation:
  
  data_num <- data %>% select(c(groupings, "samInt", "refInt", "areaRatio")) %>% 
    mutate(areaRatio = as.numeric(areaRatio))
  data_anno <- data %>% select(c(groupings, "description"))
  
  # Data wrangling per peptide:
  
  data_num <- data_num %>% ungroup() %>% 
    group_by_at(vars(groupings)) %>% 
    summarise_all(median, na.rm=T) %>% 
    filter(!is.na(samInt) | !is.na(refInt)) %>% ungroup() %>% 
    mutate(PctHeavy = 100 * (1/(1+areaRatio)))
  
  # Fix the annotation data:
  data_anno <- data_anno %>% 
    mutate(description = gsub(pattern = " ]",replacement = "",x = gsub("[","",description,fixed = T),fixed = T)) %>%
    separate_rows(description, sep = "\\s,\\s") %>% 
    mutate(description = gsub(pattern = "^\\s+",replacement = "",x = gsub("\\s+$","",description))) %>% 
    separate(col = description, into = c("uniprot","prot_name"), sep = " ", remove = F, extra = "merge") %>% 
    separate(col = prot_name, into = c(NA,"gene_sym"),sep = " GN=",remove = F) %>% 
    separate(col = gene_sym,into = c("gene_sym",NA),sep = " PE=",remove = T) %>% ungroup() %>% unique()
  
  # Merge the data and the anno back together:
  data <- right_join(data_anno, data_num, by = c(groupings))

  return(data)
} # end of function