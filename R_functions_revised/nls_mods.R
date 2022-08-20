# Create a function to test all the different methods of robust nls fitting.
# Also vary a fixed vs estimated start point for time 0. 
# Use the exponential plateau function (no estimate of the offset, set to 0/1)

# Collect the predicted values for all models next to the actual values - save 
# table for later plotting if necessary

# Create a new table to collect the following for all the test models:
# Calculate RMSE and correlation between pred and ratios to be used to determine
# best model.
# Collect the values of alpha (and R0 if not fixed), calculate the half-lives
# 
# Run this function 4 times, for each peptide per cell compartment, for the 
# cloud of peptides per protein per cell compartment and the 2 above for all
# compartments together.
# 
# Return a list with the 2 tables (data and preds, eval criteria and estimates)

nls_mods <- function(all_norm_fc, groupings, outputdir, savemod = F, Times) {
    
    # Import and reformat data:
    all_norm_fc_l <- all_norm_fc %>% 
        pivot_longer(cols = starts_with("H_"), names_to = "time", names_prefix = "H_", values_to = "ratio") %>% 
        mutate(time = as.numeric(time)) %>% 
        select(one_of(c(groupings, "time", "ratio","guess_alpha")))
    
    # Fit the different NLS models (note, CM is too slow to be useful and never converges so i removed it)
    NLS_light <- all_norm_fc_l %>% filter(Type == "TMT_light" & !is.na(ratio)) %>% 
        ungroup() %>% group_by_at(groupings) %>% 
        do(fit_port = tryCatch(nls(ratio ~ 0 + (1 - 0) * exp(-alpha * time), 
                                   data = ., 
                                   start = list(alpha = first(.$guess_alpha)),
                                   lower = list(alpha = 0),
                                   upper = list(alpha = 10),
                                   algorithm="port", 
                                   control = nls.control(minFactor=1/4096, maxiter=200)), 
                               error=function(e) NULL),
           fit_robMport = tryCatch(nlrob(ratio ~ 0 + (1 - 0) * exp(-alpha * time), 
                                         data = ., 
                                         start = list(alpha = first(.$guess_alpha)),
                                         lower = list(alpha = 0),
                                         upper = list(alpha = 10),
                                         algorithm="port", 
                                         control = nls.control(minFactor=1/4096, maxiter=200)), 
                                   error=function(e) NULL))
    
    NLS_heavy <- all_norm_fc_l %>% filter(Type == "TMT_heavy" & !is.na(ratio)) %>% 
        ungroup() %>% group_by_at(groupings) %>% 
        do(fit_port = tryCatch(nls(ratio ~ 1 + (0 - 1) * exp(-alpha * time), 
                                   data = ., 
                                   start = list(alpha = first(.$guess_alpha)),
                                   lower = list(alpha = 0),
                                   upper = list(alpha = 10),
                                   algorithm="port", 
                                   control = nls.control(minFactor=1/4096, maxiter=200)), 
                               error=function(e) NULL),
           fit_robMport = tryCatch(nlrob(ratio ~ 1 + (0 - 1) * exp(-alpha * time), 
                                         data = ., 
                                         start = list(alpha = first(.$guess_alpha)),
                                         lower = list(alpha = 0),
                                         upper = list(alpha = 10),
                                         algorithm="port", 
                                         control = nls.control(minFactor=1/4096, maxiter=200)), 
                                   error=function(e) NULL))
    
    NLS_mods <- bind_rows(NLS_light, NLS_heavy)
    rm(NLS_light, NLS_heavy)
    
    # Let's try this with a loop instead of the apply functions. Then we'll use join to put it together.
    # First to collect the alpha values.
   
    NLS_mods$alpha_robMport <- rep(NA,length(NLS_mods$fit_robMport))
    
    for (i in 1:length(NLS_mods$fit_robMport)) {
        if (!is.null(NLS_mods$fit_robMport[i][[1]])) {
            NLS_mods$alpha_robMport[i] <- coef(NLS_mods$fit_robMport[i][[1]])["alpha"]
        } else {
            NLS_mods$alpha_robMport[i] <- NA
        }
        
    }
    
    NLS_mods$alpha_port <- rep(NA,length(NLS_mods$fit_port))
    
    for (i in 1:length(NLS_mods$fit_port)) {
        if (!is.null(NLS_mods$fit_port[i][[1]])) {
            NLS_mods$alpha_port[i] <- coef(NLS_mods$fit_port[i][[1]])["alpha"]
        } else {
            NLS_mods$alpha_port[i] <- NA
        }
        
    }
    
    # Now to collect the predicted values:
    
    predict_port <- as.data.frame(matrix(nrow = length(NLS_mods$fit_port), ncol = length(Times)))
    colnames(predict_port) <- paste0("pred_port_", Times)
    NLS_mods <- bind_cols(NLS_mods, predict_port)
    rm(predict_port)
    
    for (i in 1:length(NLS_mods$fit_port)) {
        if (!is.null(NLS_mods$fit_port[i][[1]])) {
            NLS_mods[i,paste0("pred_port_", Times)] <- t(predict(NLS_mods$fit_port[i][[1]]))
        } else {
            NLS_mods[i,paste0("pred_port_", Times)] <- NA
        }
        
    }
    
    predict_robMport <- as.data.frame(matrix(nrow = length(NLS_mods$fit_robMport), ncol = length(Times)))
    colnames(predict_robMport) <- paste0("pred_robMport_", Times)
    NLS_mods <- bind_cols(NLS_mods, predict_robMport)
    rm(predict_robMport)
    
    for (i in 1:length(NLS_mods$fit_robMport)) {
        if (!is.null(NLS_mods$fit_robMport[i][[1]])) {
            NLS_mods[i,paste0("pred_robMport_", Times)] <- t(predict(NLS_mods$fit_robMport[i][[1]]))
        } else {
            NLS_mods[i,paste0("pred_robMport_", Times)] <- NA
        }
        
    }
    
    # If saving models:
    if (!missing(outputdir) & savemod) {
        dir.create(outputdir)
        save(NLS_mods,file = paste0(outputdir,"/NLS_models_results.RData"))
    }
    
    # Remove the models from the table and separate one with the alphas and another with the predictions.
    NLS_alphas <- NLS_mods %>% select(groupings, contains("alpha"))
    # Calculate half-lives:
    NLS_alphas <- NLS_alphas %>% 
        mutate(Thalf_port = log(2)/alpha_port,
               Thalf_robMport = log(2)/alpha_robMport)
    
    
    NLS_preds_port <- NLS_mods %>% select(groupings, contains("pred_port")) %>% 
        pivot_longer(cols = starts_with("pred_port"), names_prefix = "pred_port_", names_to = "time", values_to = "pred_ratio_port") %>% mutate(time = as.numeric(time))
    NLS_preds_robMport <- NLS_mods %>% select(groupings, contains("pred_robMport")) %>% 
        pivot_longer(cols = starts_with("pred_robMport"), names_prefix = "pred_robMport_", names_to = "time", values_to = "pred_ratio_robMport") %>% mutate(time = as.numeric(time))
    
    NLS_preds <- full_join(NLS_preds_port, NLS_preds_robMport, by=c(groupings,"time"))
    rm(NLS_preds_port, NLS_preds_robMport)
    
    # Join predictions and actual data, then calculate cor and rmse for each peptide:
    NLS_preds <- full_join(NLS_preds, all_norm_fc_l %>% select(-guess_alpha), by=c(groupings,"time"))
    
    NLS_cor_er <- NLS_preds %>% 
        group_by_at(groupings) %>% 
        mutate(cor_port = cor(ratio,jitter(pred_ratio_port,amount = 0.001)), 
               rmse_port = sqrt(mean((ratio-pred_ratio_port)^2)),
               cor_robMport = cor(ratio,jitter(pred_ratio_robMport,amount = 0.001)), 
               rmse_robMport = sqrt(mean((ratio-pred_ratio_robMport)^2))) %>% 
        select(-time,-ratio,-starts_with("pred_")) %>% unique()
    # Move cor and rmse (unique) per peptide to the same table as the alpha values
    NLS_alphas_cor_er <- full_join(NLS_alphas,NLS_cor_er, by = groupings)
    NLS_alphas_cor_er <- full_join(NLS_alphas_cor_er,all_norm_fc %>% select(-starts_with("H_")), by = groupings)
    rm(NLS_alphas,NLS_cor_er)

    # Save csv
    if (!missing(outputdir)) {
        dir.create(outputdir)
        write_csv(NLS_preds, path = paste0(outputdir,"/NLS_models_predictions.csv"))
        write_csv(NLS_alphas_cor_er, path = paste0(outputdir,"/NLS_models_alpha_err_half.csv"))
    }
    
    All_mods <- list(NLS_preds = NLS_preds, NLS_results = NLS_alphas_cor_er)
    
    return(All_mods)
    
}
