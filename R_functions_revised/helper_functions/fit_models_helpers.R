.get_param_NLS <- function(NLSport, param, groupings) {
    
    All_param <- data.frame(matrix(NA, nrow = dim(NLSport)[1], ncol = 1))
    colnames(All_param) <- "fit"
    
    for (i in 1:dim(NLSport)[1]) {
        
        All_param$fit[i] <- ifelse(!(is.null(coef(NLSport$fit[[i]])[param])), coef(NLSport$fit[[i]])[param], NA) 
        
    }
    colnames(All_param) <- param
    return(All_param)
}

.get_errors_Port <- function(NLSport, all_norm_fc, groupings) {
    
    All_Err <- data.frame(matrix(NA, nrow = dim(NLSport)[1], ncol = 1))
    
    colnames(All_Err) <- "Rsq"
    
    test <- str_c(paste0(groupings," == ", "NLSport$",groupings,"[i]"), collapse =" & ")
    
    for (i in 1:dim(NLSport)[1]) {
        
        data_NLS <- all_norm_fc %>% ungroup() %>% filter(eval(str2expression(test))) %>% select(time, ratio)
        nls_model_NLSPort <- NLSport$fit[[i]]
        
        if (!(is.null(nls_model_NLSPort))) {
            All_Err[i,"Rsq"] <- modelr::rsquare(model = nls_model_NLSPort, data = data_NLS)
        } else {
            All_Err[i,"Rsq"] <- c(NA)
        }
        
        rm(nls_model_NLSPort, data_NLS)
        svMisc::progress(i,dim(NLSport)[1])
    }
    
    return(All_Err)
}

ErrorsModelr <- function(modlist) {
    
    Rsquared = modelr::rsquare(model = modlist,
                               data = modlist$model %>% 
                                   rename("ratio" = "log(ratio + 0.01)") %>%
                                   mutate(ratio = exp(ratio)))
    
    Errdf <- cbind(Rsquared)
    return(Errdf)
    
}