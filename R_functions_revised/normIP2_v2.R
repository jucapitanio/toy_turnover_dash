#' The normalization function can perform the following:
#' 
#' (Row normalization) Scale TMT-quant light and heavy tables so you know which TMT tag contributed what proportion of the signal, scaling to 100%. Then multiply that by the total peak signal from the SILAC data for that peptide.
#' 
#' (Sum normalization) Ensure each timepoint column sums to the same value across the timecourse.
#' 
#' By default, the function will apply both normalizations to the data, but you may forgo either normalization if you choose to. Skipping both row and sum normalization will return the combined light/heavy TMT filtered data without modification.
#' 
#' The function will return a single, large dataframe containing both light and heavy TMT data, discarding the SILAC data
#' 
#' @param data the list of filtered data returned by \code{filterIP2}
#' @param row_norm a boolean value indicating whether the function will use SILAC peak areas to row-normalize TMT data
#' @param sum_norm a boolean value indicating whether the function will sum-normalize each timepoint column 
#' @param outputdir: (optional) directory to write out imported datasets as .csv
#' @param ... additional parameters will be passed on to \code{readr::write_csv}
#' @return a cleaned and filtered list of three dataframes: "TMT_light", "TMT_heavy", and "SILAC"

normIP2 <- function(data, row_norm = TRUE, sum_norm = TRUE, outputdir = NULL, ...) {
  
  # hardcoded timepoint column prefix, used to select tp columns
  tp_prefix = "H_"
  
  status_norm <- ""
  
  # perform row normalization using SILAC peak areas
  if (row_norm) {
    for (type in c("TMT_light", "TMT_heavy")) {
      data[[type]] <- data[[type]] %>% 
        # get the proportion at each timepoint
        ungroup() %>%
        mutate(row_signal = rowSums(select(., starts_with(tp_prefix)))) %>% 
        gather(key, value, starts_with(tp_prefix)) %>% 
        mutate(value = value / row_signal) %>% 
        select(-row_signal) %>% ungroup()
        
        # normalize by scaling to SILAC peak area
        
      data[[type]] <- inner_join(data[["SILAC"]], data[[type]],by = c("Unique", "Sequence", "Localization")) %>% 
        {
          if (type == "TMT_light") {
            ungroup(.) %>% 
              mutate(value = value * samInt) %>% 
              spread(key, value)
            
          } else {  # type must be "TMT_heavy"
            ungroup(.) %>% 
              mutate(value = value * refInt) %>% 
              spread(key, value)
          }
        } %>% select(-areaRatio, -samInt, -refInt)
     }
    
    status_norm <- "row"
  }
  
  # combine light and heavy TMT datasets
  data[["TMT_light"]]$Type <- "TMT_light"
  data[["TMT_heavy"]]$Type <- "TMT_heavy"
  combinedTMT <- bind_rows(data[["TMT_light"]], data[["TMT_heavy"]])
  
  # perform sum normalization for each timepoint, per localization
  if (sum_norm) {
    
    # decide later if this should be a helper function or if it stays here.
    .sumnormTMT <- function(combinedTMT, tp_prefix) {
      light_cc <- combinedTMT %>% filter(Type == "TMT_light")
      heavy_cc <- combinedTMT %>% filter(Type == "TMT_heavy")
      Nsum_each <- colSums(light_cc %>% select(starts_with(tp_prefix))) + colSums(heavy_cc %>% select(starts_with(tp_prefix)))
      Nsums <- median(Nsum_each) / Nsum_each
      NS_light <- sweep(light_cc %>% select(starts_with(tp_prefix)), 2, Nsums, FUN="*")
      NS_heavy <- sweep(heavy_cc %>% select(starts_with(tp_prefix)), 2, Nsums, FUN="*")
      NS_light <- bind_cols(light_cc %>% select(-starts_with(tp_prefix)), NS_light)
      NS_heavy <- bind_cols(heavy_cc %>% select(-starts_with(tp_prefix)), NS_heavy)
      NS <- bind_rows(NS_light,NS_heavy)
      return(NS)
    }
    
    # if localization exists, normalize per localization
    if (!is.null(combinedTMT$Localization)) {
      combinedTMT <- combinedTMT %>% 
        group_by(Localization) %>% 
        group_modify(~ .sumnormTMT(.x,tp_prefix))
    # if not, just normalize entire table
    } else {
      combinedTMT <- .sumnormTMT(combinedTMT,tp_prefix)
    }
      
    status_norm <- paste0(status_norm, "sum")
  }
  
  # write out the combined, normalized TMT dataset to file
  if (!missing(outputdir)) {
    
    if (!row_norm && !sum_norm) {
      status_norm <- "not"
    }
    
    write_csv(combinedTMT, paste0(outputdir, "/TMT_combined_", status_norm, "_normalized.csv"), ...)
  }
  
  return(combinedTMT)
}