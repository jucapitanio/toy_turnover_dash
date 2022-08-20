#' Filters out data columns for any missing TMT tags that were quantified by IP2 using the "Isobaric N-plex labeling (more than 6 reporter ions)" mode. This IP2 quantification method will need to be run with 10 tags, despite how many you have just because the kit used allows for 10 tags. This is very specific for rare cases when the 10-plex kit is used with 7-9 timepoints, so this function will rarely be used. You should also choose the number of columns you expect to keep, for example, if you had 8 samples quantified with the 10-plex kit, you should choose to remove 2 columns and keep 8. These filters will be applied to the 2 TMT dataframes: "TMT_light", "TMT_heavy".
#
#' PRECONDITION: the TMT-SILAC data has been imported via the \code{importIP2} function, and the dataframes in the returned list have not been renamed. Additionally, this function expects that the TMT timecourse columns have been renamed to start with the prefix "H." and rearranged in chronological order.
#' 
#' @param data the list of raw data returned by \code{importIP2}
#' @param keep_colnum an integer indicating the number of columns that should be removed
#' @param outputdir: (optional) directory to write out filtered datasets as .csv
#' @param ... additional parameters will be passed on to \code{readr::write_csv}
#' @return a cleaned and filtered list of three dataframes: "TMT_light", "TMT_heavy", and "SILAC"


rmTMTtagIP2 <- function(raw_data, Times, outputdir=NULL, ...) {
    
    # error checking code
    if (is.null(raw_data[["TMT_light"]]) || is.null(raw_data[["TMT_heavy"]]) || is.null(raw_data[["SILAC"]])) {
        stop("Argument 'raw_data' is missing one of \"TMT_light\", \"TMT_heavy\", or \"SILAC\".")
        
    } 
    raw_data[["TMT_light"]]$`TMT purity` <- as.numeric(raw_data[["TMT_light"]]$`TMT purity`)
    raw_data[["TMT_heavy"]]$`TMT purity` <- as.numeric(raw_data[["TMT_heavy"]]$`TMT purity`)
    TMT_all <- bind_rows(raw_data[["TMT_light"]], raw_data[["TMT_heavy"]], .id = "lightheavy")
    
    discard_blanks <- function(df, key, new_colnames) {
        # identify which 2 columns are blank in this run
        
        blank_cols <- df %>% 
            select(contains("m/z")) %>% 
            colSums(., na.rm = TRUE) %>% 
            sort() %>% 
            names() %>% 
            .[c(1:(10-length(Times)))]
        
        # discard blank columns and rename for consistency
        filtered <- df %>% 
            select(-blank_cols) %>% 
            rename_at(vars(contains("m/z")), funs(paste0("H_", Times)))
        
        return(filtered)
    }
    
    if (dim(TMT_all %>% select(contains("m/z")))[2] == length(Times)) {
        
        TMT_all_fixed <- TMT_all %>% 
            rename_at(vars(contains("m/z")), funs(paste0("H_", Times))) %>% 
            group_by(lightheavy) %>% 
            group_split(keep=FALSE)
        
    } else {
    
    TMT_all_fixed <- TMT_all %>% 
        group_by(lightheavy, Localization) %>% 
        group_modify(~ discard_blanks(.x, .y)) %>% 
        group_by(lightheavy) %>% 
        group_split(keep=FALSE)
    
    }
    
    # sanity check: throw error if filters discard all raw_data
    if (any(c(dim(TMT_all_fixed[[1]]), dim(TMT_all_fixed[[2]])) == 0)) {
        stop("Specified filters discarded all observations. Please check function parameters to ensure accurate use.")
    }
    
    # put the reorganized time-course back with raw_data
    raw_data[["TMT_light"]] = TMT_all_fixed[[1]]
    raw_data[["TMT_heavy"]] = TMT_all_fixed[[2]]
    
    # write out each raw_dataset to a different file
    if (!missing(outputdir)) {
        for(df in c("TMT_light", "TMT_heavy")) {
            write_csv(raw_data[[df]], paste0(outputdir, "/", df, "_TMTtag_filtered.csv"), ...)
        }
    }
    
    return(raw_data)
}
    