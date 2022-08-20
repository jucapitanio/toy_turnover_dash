FC_calc <- function(all_normalized_tmt, Times) {
    
    # Normalize the light labelled peptides:
    NS_light <- all_normalized_tmt %>% filter(Type == "TMT_light") %>% ungroup()
    NS_lightFC <- sweep(as.matrix(NS_light %>% select(starts_with("H_"))), 1, as.matrix(NS_light %>% select(H_0)), "/")
    NS_lightFC <- bind_cols(NS_light %>% select(-starts_with("H_")), as.data.frame(NS_lightFC)[,paste0("H_", Times)]) %>% drop_na(starts_with("H_"))
    
    # Normalize the heavy labelled peptides:
    NS_heavy <- all_normalized_tmt %>% filter(Type == "TMT_heavy") %>% ungroup() %>% select(-PctHeavy)
    NS_all <- inner_join(NS_light %>% select(-Type), NS_heavy, by= c("Unique", "Sequence", "uniprot", "description", "Localization","prot_name","gene_sym"), suffix = c(".light", ".heavy"))
    temp_rat <- as.data.frame(lapply(Times, function(x) NS_all[,paste0("H_",x,".heavy")] / (NS_all[,paste0("H_",x,".heavy")] + NS_all[, paste0("H_",x,".light")])))
    colnames(temp_rat) <- paste0("H_",Times)
    NS_heavyFC <- bind_cols(NS_all %>% select(-starts_with("H_")), temp_rat) %>% drop_na(starts_with("H_"))
    
     # combine heavy and light:
    all_FC <- bind_rows(NS_lightFC, NS_heavyFC)
    
    # save data used to fit models:
    #if (exists("outputdir")) {
    #    write_csv(all_FC, paste0(outputdir, "/", "final data for model fit.csv"))
    #}
    
    return(all_FC)
}