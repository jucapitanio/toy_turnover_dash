#' Reads in TMT-SILAC protein datasets from IP2 and filters the datasets for the peptide-level data only.
#' 
#' This function expects TMT-light, TMT-heavy, and SILAC data from the same experiment, and expects the same number of files for each data type.
#' 
#' @param files a character vector of file paths pointing to IP2 data (tab-delimited)
#' @param datatype a character vector specifying \code{“light”} or \code{“heavy”} or \code{“silac”} corresponding to the files specified
#' @param localization (optional) a character vector specifying cellular localization (e.g. a cellular fraction) that each dataset comes from. Will add a “localization” column to the combined data.
#' @param outputdir: (optional) directory to write out imported datasets as .csv
#' @param ... additional parameters will be passed on to \code{readr::write_csv}
#' @return a list of the imported data tables, summarized by \code{datatype}

importIP2 <- function(files, datatype, localization=NULL, outputdir=NULL, ...) {

  #for troubleshooting:
  # files = design$files
  # datatype = design$datatype
  # localization = design$localization
  #
  # check for invalid required arguments
  if (length(files) != length(datatype)) {
    stop("Number of files to import does not match the number of datatype")
  } else {
    
    table <- table(datatype)
    if (dim(table) != 3 ||
        !all(c("light", "heavy", "silac") %in% datatype)) {
      stop("Please supply 'light', 'heavy', and 'silac' data")
    } else if (!all(sapply(table, function(x) {
      x == table[3]
    }))) {
      stop("Each datatype must have the same number of files")
    }
    rm(table)
  }
  
  # check for invalid optional argument
  if (!missing(localization) &&
      !all(sapply(list(files, datatype, localization), function(x) {
        length(x) == length(localization)
      }))) {
    stop("Number of localizations does not match number of files or datatype")
  }
  
  #go through all files and import via corresponding helper functions
  imported_data <- (mapply(
    path = files,
    type = datatype,
    FUN = function(path, type) {
      if (type == "light" || type == "heavy") {
        return(.readTMT(path))
        
      } else if (type == "silac") {
        return(.readSILAC(path))
      }
    },
    SIMPLIFY = F
  ))
  
  # add localization column to the datasets
  if (!missing(localization)) {
    imported_data <- mapply(
      df = imported_data,
      loc = localization,
      
      FUN = function(df, loc) {
        df %>% 
          mutate(Localization = loc)
      }
    )
    
    # merge different localizations into a single dataframe, by datatype
    TMT_light = data_frame()
    TMT_heavy = data_frame()
    SILAC = data_frame()
    
    for (i in 1:length(files)) {
      if (datatype[i] == "light") {
        TMT_light <-
          rbind(imported_data[[files[i]]], TMT_light)
        
      } else if (datatype[i] == "heavy") {
        TMT_heavy <-
          rbind(imported_data[[files[i]]], TMT_heavy)
        
      } else if (datatype[i] == "silac") {
        SILAC <-
          rbind(imported_data[[files[i]]], SILAC)
      }
    }
    
    imported_data <- list(
      TMT_light = TMT_light,
      TMT_heavy = TMT_heavy,
      SILAC = SILAC
    )
  }
  
  if (!missing(outputdir)) {
    # write out each dataset to a different file
    for(df in c("TMT_light", "TMT_heavy", "SILAC")) {
      write_csv(imported_data[[df]], paste0(outputdir, "/", df, "_raw.csv"), ...)
    }
  }
  
  return(imported_data)
}