#' This is a helper function for \code{filterIP2} which will transform some of the TMT data in order to throw out some unreliable data points. The function assumes that the timepoint columns in the data have been properly renamed and rearranged in chronological order.
#' 
#' This function does a few things to the data:
#' 
#' - for identical peptide sequences, intensities of the same TMT channel derived from different evidence entries are summed
#' 
#' - for light TMT, throw out sequences where the first timepoint reads zero
#' 
#' - for heavy TMT, throw out sequences where the last timepoint reads zero
#' 
#' - throw out sequences that, at a given timepoint, reach zero intensity in-between two nonzero intensities
#' 
#' Since light-TMT is a decay function, throw out a light-TMT sequence if the current timepoint reads 0 and a future timepoint reads nonzero. Likewise, since heavy-TMT is a growth function, throw out a heavy-TMT sequence if the current timepoint reads 0 and a past timepoint reads nonzero.
#' 
#' @param data a dataframe containing TMT data from the IP2 pipeline, where timepoints are prefixed by "H_"
#' @param type "light" or "heavy" TMT data
#' @return a dataframe of TMT data with some datapoints thrown out

.filterTMT <- function(data, type) {
  # error checking; assumes that data is formatted correctly
  if (type != "light" && type != "heavy") {
    stop("Argument 'type' must be one of \"light\" or \"heavy\"")
  }
  
  # hardcoded timepoint column prefix, used to select tp columns
  tp_prefix = "H_"
  
  # rename timecourse columns to something concrete
  names_map <- data %>%
    select(starts_with(tp_prefix)) %>%
    colnames(.) %>% as.data.frame(.)
  colnames(names_map) <- "old_names"
  names_map <- names_map %>% 
    separate(old_names, into = c(NA, "time_pt"),sep = "_",remove = F) %>% 
    arrange(as.numeric(time_pt)) %>% 
    mutate(placeholder = paste0("T", seq(1, length(old_names))))
  
  colnames(data) <- 
    plyr::mapvalues(colnames(data), names_map$old_names, names_map$placeholder)
  
  maxtime <- paste0("T", nrow(names_map))
  timepoints <- seq(nrow(names_map)-1)
  
  # summarize identical sequences
  groupings <- c("Unique", "Sequence")
  if (!is.null(data$Localization)) {
    groupings <- c(groupings, "Localization")
  }
  data <- data %>%
    select(groupings, starts_with("T")) %>% 
    ungroup() %>%
    group_by_at(vars(groupings)) %>%
    summarize_all(sum) %>%
    mutate(remove = "no")
  
  # throw out some data
  if (type == "light") {
    data <- data %>%
      filter(T1 != 0)
    
    # look for zero intensities in the middle of the time course, throw out
    for (row in seq(nrow(data))) {
      for (tp in timepoints) {
        curr_tp <- paste0("T", tp)
        next_tp <- paste0("T", tp + 1)
        
        if (data[row, curr_tp] == 0 &&
            data[row, next_tp] != 0) {
          data$remove[row] <- "yes"
        }
      }
    }
    
  } else {
    # type == "heavy" must be true
    data <- data %>%
      filter(get(maxtime) != 0)
    
    # look for zero intensities in the middle of the time course, throw out
    for (row in seq(nrow(data))) {
      for (tp in timepoints) {
        curr_tp <- paste0("T", tp)
        next_tp <- paste0("T", tp + 1)
        
        if (data[row, curr_tp] != 0 &&
            data[row, next_tp] == 0) {
          data$remove[row] <- "yes"
        }
      }
    }
  }
  
  data <- data %>%
    filter(remove != "yes") %>%
    select(-remove, -`TMT purity`, -`TMT S/N`) %>%
    `colnames<-`(plyr::mapvalues(
      colnames(.), names_map$placeholder, as.character(names_map$old_names)
    ))

  return(data)
} # end of function