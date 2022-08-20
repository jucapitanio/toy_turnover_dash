# This is a function I made to estimate the half-life of the proteins based on the percentage of heavy proteins present across the entire summed timecourse.
# It applies a grid-search like idea so that we can use the data from MS1 and not only that from MS3 for our estimates.

guess_Thalf_pctHeavy <- function(Times, all_norm_fc) {
    
    ### Calculate a guestimate Thalf using an assymptote with no offset
    
    # Let's do a loop to calculate the pcnt_heavy for half-lives from 0.1h to 100000h using the formula from a simple exponential model with no offset:
    Thalf <- c(seq(0.1,0.9,0.1),seq(1,100000,0.5))
    pct_heavy <- c()
    Times <- as.numeric(Times)
    
    for (i in 1:length(Thalf)) {
        alpha <- log(2)/Thalf[i]
        light <- 0 + (100 - 0) * exp(-alpha * Times)
        heavy <- 100 + (0 - 100) * exp(-alpha * Times)
        pct_heavy[i] <- sum(heavy) / (sum(light)+sum(heavy))
    }
    
    df <- as.data.frame(matrix(nrow = length(Thalf),ncol = 1))
    df$Thalf <- as.numeric(Thalf)
    df$pcnt_heavy = as.numeric(pct_heavy)
    # Reordering needed to use with match_closest below.
    df <- df %>% select(-V1) %>% arrange(pcnt_heavy)
    
    
    # figure out the Thalf value given the pcnt heavy in the entire timecourse?
    source("R_functions_revised/helper_functions/match_closest.R")
    
    match_thalf <- function(x_val){df[match_closest(x_val/100, df$pcnt_heavy), "Thalf"]}
    match_thalf <- Vectorize(match_thalf)
    
    all_norm_fc <- all_norm_fc %>% 
        mutate(guess_Thalf = match_thalf(PctHeavy)) %>%
        mutate(guess_alpha = alpha <- log(2)/guess_Thalf)
    
    return(all_norm_fc)
    
}
