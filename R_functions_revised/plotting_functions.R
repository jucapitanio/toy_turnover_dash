pcnt_label_plot <- function(normalized_tmt, Localization_facet = F) {
  
  if (Localization_facet) {
    Norm <- normalized_tmt %>% ungroup() %>% group_by(Localization, Type) %>% 
      summarise_at(.vars = vars(starts_with("H_")),.funs = sum)
  } else {
    Norm <- normalized_tmt %>% ungroup() %>% group_by(Type) %>% 
      summarise_at(.vars = vars(starts_with("H_")),.funs = sum)
  }
  
  Norm <- Norm %>% pivot_longer(cols = starts_with("H_"), names_to = "Times", names_prefix = "H_", values_to = "intensity")
  
  ggplot(data = Norm, aes(fill = Type, y = intensity, x = reorder(Times, as.numeric(Times)))) +
    geom_bar(stat = "identity") +
    xlab("Time (hours)") +
    ylab("Intensities (sum)") +
    {if(Localization_facet)facet_wrap(.~ Localization, scales = "free_y",nrow = 1)} +
    theme(legend.position = "bottom", legend.title = element_blank())
  
}


proteome_mean_plot <- function(all_norm_fc, Localization_facet = F) {
  if (Localization_facet) {
    Prot_means <- all_norm_fc %>%
      group_by(Localization, Type) %>%
      summarise_at(.vars = vars(starts_with("H_")), mean, na.rm = T) %>%
      gather(key = "Hours", value = "ratio", -Localization, -Type) %>%
      mutate(Hours = as.numeric(str_extract(Hours, "\\-*\\d+\\.*\\d*")))

    Prot_CIs <- all_norm_fc %>%
      group_by(Localization, Type) %>%
      summarise_at(.vars = vars(starts_with("H_")), function(x) {
        DescTools::MeanCI(x, conf.level = 0.95, sides = "two.sided", na.rm = T)["mean"] - DescTools::MeanCI(x, conf.level = 0.95, sides = "two.sided", na.rm = T)["lwr.ci"]
      }) %>%
      gather(key = "Hours", value = "ratio", -Localization, -Type) %>%
      mutate(Hours = as.numeric(str_extract(Hours, "\\-*\\d+\\.*\\d*")))

    Proteome <- full_join(Prot_means, Prot_CIs, by = c("Localization", "Type", "Hours"), suffix = c(".mean", ".CI95"))

    ggplot(Proteome, aes(x = Hours, y = ratio.mean, colour = Type, group = Type)) +
      geom_errorbar(aes(ymin = ratio.mean - ratio.CI95, ymax = ratio.mean + ratio.CI95), colour = "black", width = .2, position = position_dodge(0.1)) +
      geom_line(position = position_dodge(0.1)) +
      geom_point(position = position_dodge(0.1), size = .1, shape = 19) + # 21 is filled circle
      xlab("Time (Hours)") +
      ylab("Proportion of Labeled peptides") +
      scale_colour_hue(name = NULL, l = 40) +
      theme_bw() +
      theme(legend.position = "none") +
      facet_wrap(. ~ Localization, nrow = 1)
  } else {
    Prot_means <- all_norm_fc %>%
      group_by(Type) %>%
      summarise_at(.vars = vars(starts_with("H_")), mean, na.rm = T) %>%
      gather(key = "Hours", value = "ratio", -Type) %>%
      mutate(Hours = as.numeric(str_extract(Hours, "\\-*\\d+\\.*\\d*")))

    Prot_CIs <- all_norm_fc %>%
      group_by(Type) %>%
      summarise_at(.vars = vars(starts_with("H_")), function(x) {
        DescTools::MeanCI(x, conf.level = 0.95, sides = "two.sided", na.rm = T)["mean"] - DescTools::MeanCI(x, conf.level = 0.95, sides = "two.sided", na.rm = T)["lwr.ci"]
      }) %>%
      gather(key = "Hours", value = "ratio", -Type) %>%
      mutate(Hours = as.numeric(str_extract(Hours, "\\-*\\d+\\.*\\d*")))

    Proteome <- full_join(Prot_means, Prot_CIs, by = c("Type", "Hours"), suffix = c(".mean", ".CI95"))

    ggplot(Proteome, aes(x = Hours, y = ratio.mean, colour = Type, group = Type)) +
      geom_errorbar(aes(ymin = ratio.mean - ratio.CI95, ymax = ratio.mean + ratio.CI95), colour = "black", width = .2, position = position_dodge(0.1)) +
      #geom_smooth(method="nls", formula=y~exp(x*b), method.args=list(start=c(b=0.5), control=nls.control(maxiter=2000)), se=FALSE) +
      geom_line(position = position_dodge(0.1)) +
      geom_point(position = position_dodge(0.1), size = .1, shape = 19) + # 21 is filled circle
      xlab("Time (Hours)") +
      ylab("Proportion of Labeled peptides") +
      scale_colour_hue(name = NULL, l = 40) +
      theme_bw() +
      theme(legend.position = "none")
  }
}

proteome_mean_plot2 <- function(all_norm_fc, Localization_facet = F) {
  if (Localization_facet) {
    Prot_means <- all_norm_fc %>%
      group_by(Localization, Type) %>%
      summarise_at(.vars = vars(starts_with("H_")), mean, na.rm = T) %>%
      gather(key = "Hours", value = "ratio", -Localization, -Type) %>%
      mutate(Hours = as.numeric(str_extract(Hours, "\\-*\\d+\\.*\\d*")))
    
    Prot_CIs <- all_norm_fc %>%
      group_by(Localization, Type) %>%
      summarise_at(.vars = vars(starts_with("H_")), function(x) {
        DescTools::MeanCI(x, conf.level = 0.95, sides = "two.sided", na.rm = T)["mean"] - DescTools::MeanCI(x, conf.level = 0.95, sides = "two.sided", na.rm = T)["lwr.ci"]
      }) %>%
      gather(key = "Hours", value = "ratio", -Localization, -Type) %>%
      mutate(Hours = as.numeric(str_extract(Hours, "\\-*\\d+\\.*\\d*")))
    
    Proteome <- full_join(Prot_means, Prot_CIs, by = c("Localization", "Type", "Hours"), suffix = c(".mean", ".CI95"))
    
    ggplot(Proteome, aes(x = Hours, y = ratio.mean, colour = Type, group = Type)) +
      geom_errorbar(aes(ymin = ratio.mean - ratio.CI95, ymax = ratio.mean + ratio.CI95), colour = "black", width = .2, position = position_dodge(0.1)) +
      geom_smooth(method="nls", formula=y ~ 0 + (1 - 0) * exp(-alpha * x), method.args=list(start=c(alpha=0.05), control=nls.control(maxiter=2000)), se=FALSE,data = Proteome %>% filter(Type == "TMT_light")) +
      geom_smooth(method="nls", formula=y ~ 1 + (0 - 1) * exp(-alpha * x), method.args=list(start=c(alpha=0.05), control=nls.control(maxiter=2000)), se=FALSE,data = Proteome %>% filter(Type == "TMT_heavy")) +
      #geom_line(position = position_dodge(0.1)) +
      geom_point() + # 21 is filled circle
      xlab("Time (Hours)") +
      ylab("Proportion of Labeled peptides") +
      scale_colour_hue(name = NULL, l = 40) +
      geom_point() +
      theme_bw() +
      theme(legend.position = "none") +
      facet_wrap(. ~ Localization, nrow = 1)
  } else {
    Prot_means <- all_norm_fc %>%
      group_by(Type) %>%
      summarise_at(.vars = vars(starts_with("H_")), mean, na.rm = T) %>%
      gather(key = "Hours", value = "ratio", -Type) %>%
      mutate(Hours = as.numeric(str_extract(Hours, "\\-*\\d+\\.*\\d*")))
    
    Prot_CIs <- all_norm_fc %>%
      group_by(Type) %>%
      summarise_at(.vars = vars(starts_with("H_")), function(x) {
        DescTools::MeanCI(x, conf.level = 0.95, sides = "two.sided", na.rm = T)["mean"] - DescTools::MeanCI(x, conf.level = 0.95, sides = "two.sided", na.rm = T)["lwr.ci"]
      }) %>%
      gather(key = "Hours", value = "ratio", -Type) %>%
      mutate(Hours = as.numeric(str_extract(Hours, "\\-*\\d+\\.*\\d*")))
    
    Proteome <- full_join(Prot_means, Prot_CIs, by = c("Type", "Hours"), suffix = c(".mean", ".CI95"))
    
    ggplot(Proteome, aes(x = Hours, y = ratio.mean, colour = Type, group = Type)) +
      geom_errorbar(aes(ymin = ratio.mean - ratio.CI95, ymax = ratio.mean + ratio.CI95), colour = "black", width = .2, position = position_dodge(0.1)) +
      geom_smooth(method="nls", formula=y ~ 0 + (1 - 0) * exp(-alpha * x), method.args=list(start=c(alpha=0.05), control=nls.control(maxiter=2000)), se=FALSE,data = Proteome %>% filter(Type == "TMT_light")) +
      geom_smooth(method="nls", formula=y ~ 1 + (0 - 1) * exp(-alpha * x), method.args=list(start=c(alpha=0.05), control=nls.control(maxiter=2000)), se=FALSE,data = Proteome %>% filter(Type == "TMT_heavy")) +
      #geom_line(position = position_dodge(0.1)) +
      geom_point() + # 21 is filled circle
      xlab("Time (Hours)") +
      ylab("Proportion of Labeled peptides") +
      scale_colour_hue(name = NULL, l = 40) +
      theme_bw() +
      theme(legend.position = "none")
  }
}
