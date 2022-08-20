#' Filters out known contaminants in additon to any additional filters supplied to the optional argument. These filters will be applied to all three dataframes: "TMT_light", "TMT_heavy", and "SILAC" data. The following contaminants are always filtered for: "Reverse_sp", "contaminant_" in the `uniprot` column, and "Bos taurus", "keratin" in the `description` column.
#' 
#' Additionally, this function will throw out the following data points:
#' 
#' For TMT data, light peptides are discarded if the first timepoint is 0 and heavy peptides are discarded if the last timepoint is 0. Further, a peptide sequence is discarded when a a peptide reaches zero intensity in-between two nonzero intensities at a given timepoint (an impossible measurement).
#' 
#' For SILAC data, a peptide sequence is discarded if its regression factor is equal to -1 (i.e. the IP2 calculation failed for this sequence). Additionally, the function will throw out singleton peptide data with regression scores below 0.8, UNLESS there is no regular peptide data for that sequence; any remaining duplicate singleton/regular peptides are averaged.
#' 
#' PRECONDITION: the TMT-SILAC data has been imported via the \code{importIP2} function, and the dataframes in the returned list have not been renamed. Additionally, this function expects that the TMT timecourse columns have been renamed to start with the prefix "H_" and rearranged in chronological order.
#' 
#' @param data the list of raw data returned by \code{importIP2}
#' @param filter_uniprot (optional) a character vector of exact uniprots to filter out of the uniprot column
#' @param filter_description (optional) a character vector of keywords to filter out of the description column
#' @param outputdir: (optional) directory to write out imported datasets as .csv
#' @param ... additional parameters will be passed on to \code{readr::write_csv}
#' @return a cleaned and filtered list of three dataframes: "TMT_light", "TMT_heavy", and "SILAC"

filterIP2 <- function(data, filter_uniprot = NULL, filter_description = NULL, outputdir = NULL, ...) {
  # error checking code
  if (is.null(data[["TMT_light"]]) || is.null(data[["TMT_heavy"]]) || is.null(data[["SILAC"]])) {
    stop("Argument 'data' is missing one of \"TMT_light\", \"TMT_heavy\", or \"SILAC\".")
    
  } else if (data[["TMT_light"]] %>% select(starts_with("H_")) %>% is_empty() ||
             data[["TMT_heavy"]] %>% select(starts_with("H_")) %>% is_empty()) {
    stop("Timepoint columns of \"TMT_light\", \"TMT_heavy\" dataframes in 'data' have not been renamed to begin with prefix \"H_\"")
  }
  
  # filter out contaminants
  data <- lapply(data, function(df) {
    df %>% 
      filter(
        !(grepl(pattern = "no description", x = df$description, ignore.case = T) |
            grepl(pattern = "contaminant_", x = df$uniprot, ignore.case = T) |
            grepl(pattern = "Reverse_sp", x = df$uniprot, ignore.case = T) |
            grepl(pattern = "Bos taurus", x=df$description, ignore.case = T) |
            grepl(pattern = "keratin", x = df$description, ignore.case = T))
      )
  })
  
  if (!missing(filter_uniprot)) {
    data <- lapply(data, function(df) {
      df %>%
        filter( !(toupper(uniprot) %in% toupper(filter_uniprot)) )
    })
  }

  if (!missing(filter_description)) {
    regex <- paste(filter_description, collapse = "|")
    data <- lapply(data, function(df) {
      df %>% 
        filter(!(grepl(pattern = regex, x = df$description, ignore.case = T)))
    })
  }
  
  # sanity check: throw error if filters discard all data
  if (any(c(dim(data[[1]]), dim(data[[2]]), dim(data[[3]])) == 0)) {
    stop("Specified filters discarded all observations. Please provide a different set of filters.")
  }
  
  # throw out some data points
  data[["TMT_light"]] <- .filterTMT(data[["TMT_light"]], "light")
  data[["TMT_heavy"]] <- .filterTMT(data[["TMT_heavy"]], "heavy")
  suppressWarnings(data[["SILAC"]] <- .filterSILAC(data[["SILAC"]]))

  # write out each dataset to a different file
  if (!missing(outputdir)) {
    for(df in c("TMT_light", "TMT_heavy", "SILAC")) {
      write_csv(data[[df]], paste0(outputdir, "/", df, "_filtered.csv"), ...)
      save(data, file = paste0(outputdir, "/", "all_data_filtered.csv"))
    }
  }
  
  return(data)
} # end of function