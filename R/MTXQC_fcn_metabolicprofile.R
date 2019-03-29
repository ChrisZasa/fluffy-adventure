
metabolic_profile_clean = function(dataframe, groups, value_mean) {
  #' Determination of the metabolic profile
  #' 
  #' The metabolic profile summarises the quantities 
  #' of different derivates of the same metabolite.
  #' Specify the grouping of your statistcs in the vector
  #' group
  #' 
  #' @param dataframe containing information about quantities associated with metabolite names and corresponding derivates.
  #'
  #' 
  #' @param groups conditions that should be used for determining statistics
  #' 
  #' @param value_mean quantitative information that should be used for calculation of statistics
  #'
  #' @return dataframe containing statistics of metabolic profile regarding defined conditions
  #' 
  #'   
  
  #df = quant_profile
  #groups = parvec
  #value_mean = "Q_ss"
  
  
  ### 1 - determine quant_profile for Derivate
  
  groups_1 = c("Letter_Derivate", groups)
  
  quant_derivate = ddply(df, groups_1, .fun = function(xx) {
    c(n_Deriv = length(xx[, value_mean]),
      m_Deriv = mean(xx[, value_mean], na.rm = TRUE),
      sd_Deriv = sd(xx[, value_mean], na.rm = TRUE))
  })
  
  temp = df[,c("Letter_Derivate", "Pathway", groups, "sd_Q_ss")]
  q_merge = merge(quant_derivate, temp)
  
  ### 2 - Replace missing sd values   
  q_merge$sd_man = ifelse(is.na(q_merge$sd_Deriv), 
                          q_merge$sd_Q_ss, q_merge$sd_Deriv)
  
  # clean data set & export
  profile_sel = unique(q_merge[,c('Letter_Derivate','Pathway', 
                                  groups ,'n_Deriv','m_Deriv','sd_man')])
  
  write.csv(profile_sel, paste0(path_setup, set_output, set_pp, 
                                "Metabolite_profile_summary_", 
                                params$quant, ".csv"), row.names = F)
  
  
  ### 3 - Summarise for each condition
  group_2 = c(groups, "Pathway")
  
  sum_class = ddply(profile_sel, group_2, summarise, 
                    sum_c = sum(m_Deriv))
  
  profile_plot = merge(profile_sel, sum_class)
  
  #fraction of met. within pathway
  profile_plot$frac_pw = profile_plot$m_Deriv / profile_plot$sum_c * 100
  
  #sum all metabolite per condition
  sum_global = ddply(profile_plot, groups, summarise, 
                     global_all = sum(m_Deriv))
  
  profile_plot = merge(profile_plot, sum_global)
  profile_plot$frac_global = profile_plot$m_Deriv / profile_plot$global_all * 100
  
  #discard pathway fraction <1%
  profile_plot_sel = subset(profile_plot, profile_plot$frac_pw >= 1.00)
  
  if (length(groups) == 1) {
    print( 
      ggplot(profile_plot_sel, aes(frac_pw, frac_global,
                                   ymax = 1.3 * frac_global, size = m_Deriv, 
                                   label = Letter_Derivate)) +
        geom_point(aes(fill = Pathway), shape = 21) +
        scale_size_area(max_size = 20) +
        theme_bw() +
        geom_text_repel(aes(frac_pw, frac_global, label = Letter_Derivate),
                        size = 3,  
                        box.padding = unit(0.8, 'lines'), 
                        point.padding = unit(1.6, 'lines'), 
                        segment.color = '#cccccc',
                        segment.size = 0.5,  
                        force = 2,
                        max.iter = 3e3) +
        scale_x_log10() +
        facet_grid(as.formula(paste("~",groups[1]))) +
        scale_fill_mtxqc(values = "hot") +
        #ggtitle(paste0('Metabolic profile of: ', var, " (", params$quant, ", ", params$analysis)) +
        ggtitle(paste('Metabolic profile based on (ls): ', params$quant)) +
        xlab('Fraction of metabolite within its pathway in (%)') +
        ylab('Fraction of metabolite per total metabolite content in (%)') +
        theme(strip.background = element_rect(fill = "white"))
    )
  }
  
  
  
  if (length(groups) == 2) {
    print( 
      ggplot(profile_plot_sel, aes(frac_pw, frac_global,
                                   ymax = 1.3 * frac_global, size = m_Deriv, 
                                   label = Letter_Derivate)) +
        geom_point(aes(fill = Pathway), shape = 21) +
        scale_size_area(max_size = 20) +
        theme_bw() +
        geom_text_repel(aes(frac_pw, frac_global, label = Letter_Derivate),
                        size = 3,  
                        box.padding = unit(0.8, 'lines'), 
                        point.padding = unit(1.6, 'lines'), 
                        segment.color = '#cccccc',
                        segment.size = 0.5,  
                        force = 2,
                        max.iter = 3e3) +
        scale_x_log10() +
        facet_grid(as.formula(paste(groups[1],"~",groups[2]))) +
        scale_fill_mtxqc(palette = "hot") +
        #ggtitle(paste0('Metabolic profile of: ', var, " (", params$quant, ", ", params$analysis)) +
        ggtitle(paste('Metabolic profile based on (ls): ', params$quant)) +
        xlab('Fraction of metabolite within its pathway in (%)') +
        ylab('Fraction of metabolite per total metabolite content in (%)') +
        theme(strip.background = element_rect(fill = "white"))
    )
  }
  
  # if (length(parvec) = 3) {
  #   
  #   for (var in unique(profile_plot[parvec[1]])) {
  #     print( 
  #       ggplot( subset(profile_plot, parvec[1] %in% var),
  #                 #subset(profile_plot, profile_plot[parvec[1]] == "sample"),
  #                     aes(frac_pw, frac_global,
  #                         ymax = 1.3 * frac_global, size = m_Deriv,
  #                         label = Letter_Derivate)) +
  #         geom_point(aes(fill = Pathway), shape = 21) +
  #         scale_size_area(max_size = 20) +
  #         theme_bw() +
  #         geom_text_repel(aes(frac_pw, frac_global, label = Letter_Derivate),
  #           size = 3,
  #           box.padding = unit(0.8, 'lines'),
  #           point.padding = unit(1.6, 'lines'),
  #           segment.color = '#cccccc',
  #           segment.size = 0.5,
  #           force = 2,
  #           max.iter = 3e3) +
  #         scale_x_log10() +
  #         facet_grid(as.formula(paste(parvec[1],"~",parvec[2]))) +
  #         scale_fill_manual(values = color_pathway) +
  #         #ggtitle(paste0('Metabolic profile of: ', var, " (", params$quant, ", ", params$analysis)) +
  #         ggtitle(paste('Metabolic profile based on (ls): ', params$quant)) +
  #         xlab('Fraction of metabolite within its pathway in (%)') +
  #         ylab('Fraction of metabolite per total metabolite content in (%)') +
  #         theme(strip.background = element_rect(fill = "white"))
  #     )
  #   }
  # }
  
  return(profile_plot_sel)
  message("Metabolic profile calculated and exported. Plot only printed in case of 2 variable condition.")
}



metabolic_profile = function(df, groups, value_mean) {
  #' Determination of the metabolic profile
  #' The metabolic profile summarises the quantities 
  #' of different derivates of the same metabolite.
  #' Specify the grouping of your statistcs in the vector
  #' group
  
  #df = quant_profile
  #groups = parvec
  #value_mean = "Q_ss"
  
  
  ### 1 - determine quant_profile for Derivate
  
  groups_1 = c("Letter_Derivate", groups)
  
  quant_derivate = ddply(df, groups_1, .fun = function(xx) {
    c(n_Deriv = length(xx[, value_mean]),
      m_Deriv = mean(xx[, value_mean], na.rm = TRUE),
      sd_Deriv = sd(xx[, value_mean], na.rm = TRUE))
  })
  
  temp = df[,c("Letter_Derivate", "Pathway", groups, "sd_Q_ss")]
  q_merge = merge(quant_derivate, temp)
  
  ### 2 - Replace missing sd values   
  q_merge$sd_man = ifelse(is.na(q_merge$sd_Deriv), 
                          q_merge$sd_Q_ss, q_merge$sd_Deriv)
  
  # clean data set & export
  profile_sel = unique(q_merge[,c('Letter_Derivate','Pathway', 
                                  groups ,'n_Deriv','m_Deriv','sd_man')])
  
  write.csv(profile_sel, paste0(path_setup, set_output, set_pp, 
                                "Metabolite_profile_summary_", 
                                params$quant, ".csv"), row.names = F)
  
  
  ### 3 - Summarise for each condition
  group_2 = c(groups, "Pathway")
  
  sum_class = ddply(profile_sel, group_2, summarise, 
                    sum_c = sum(m_Deriv))
  
  profile_plot = merge(profile_sel, sum_class)
  
  #fraction of met. within pathway
  profile_plot$frac_pw = profile_plot$m_Deriv / profile_plot$sum_c * 100
  
  #sum all metabolite per condition
  sum_global = ddply(profile_plot, groups, summarise, 
                     global_all = sum(m_Deriv))
  
  profile_plot = merge(profile_plot, sum_global)
  profile_plot$frac_global = profile_plot$m_Deriv / profile_plot$global_all * 100
  
  #discard pathway fraction <1%
  profile_plot_sel = subset(profile_plot, profile_plot$frac_pw >= 1.00)
  
  if (length(groups) == 1) {
    print( 
      ggplot(profile_plot_sel, aes(frac_pw, frac_global,
                                   ymax = 1.3 * frac_global, size = m_Deriv, 
                                   label = Letter_Derivate)) +
        geom_point(aes(fill = Pathway), shape = 21) +
        scale_size_area(max_size = 20) +
        theme_bw() +
        geom_text_repel(aes(frac_pw, frac_global, label = Letter_Derivate),
                        size = 3,  
                        box.padding = unit(0.8, 'lines'), 
                        point.padding = unit(1.6, 'lines'), 
                        segment.color = '#cccccc',
                        segment.size = 0.5,  
                        force = 2,
                        max.iter = 3e3) +
        scale_x_log10() +
        facet_grid(as.formula(paste("~",groups[1]))) +
        scale_fill_mtxqc(palette = mypal) +
        #ggtitle(paste0('Metabolic profile of: ', var, " (", params$quant, ", ", params$analysis)) +
        ggtitle(paste('Metabolic profile based on (ls): ', params$quant)) +
        xlab('Fraction of metabolite within its pathway in (%)') +
        ylab('Fraction of metabolite per total metabolite content in (%)') +
        theme(strip.background = element_rect(fill = "white"))
    )
  }
  
  
  if (length(groups) == 2) {
    print( 
      ggplot(profile_plot_sel, aes(frac_pw, frac_global,
                                   ymax = 1.3 * frac_global, size = m_Deriv, 
                                   label = Letter_Derivate)) +
        geom_point(aes(fill = Pathway), shape = 21) +
        scale_size_area(max_size = 20) +
        theme_bw() +
        geom_text_repel(aes(frac_pw, frac_global, label = Letter_Derivate),
                        size = 3,  
                        box.padding = unit(0.8, 'lines'), 
                        point.padding = unit(1.6, 'lines'), 
                        segment.color = '#cccccc',
                        segment.size = 0.5,  
                        force = 2,
                        max.iter = 3e3) +
        scale_x_log10() +
        facet_grid(as.formula(paste(groups[1],"~",groups[2]))) +
        scale_fill_mtxqc(palette = mypal) +
        #ggtitle(paste0('Metabolic profile of: ', var, " (", params$quant, ", ", params$analysis)) +
        ggtitle(paste('Metabolic profile based on (ls): ', params$quant)) +
        xlab('Fraction of metabolite within its pathway in (%)') +
        ylab('Fraction of metabolite per total metabolite content in (%)') +
        theme(strip.background = element_rect(fill = "white"))
    )
  }
  
  # if (length(parvec) = 3) {
  #   
  #   for (var in unique(profile_plot[parvec[1]])) {
  #     print( 
  #       ggplot( subset(profile_plot, parvec[1] %in% var),
  #                 #subset(profile_plot, profile_plot[parvec[1]] == "sample"),
  #                     aes(frac_pw, frac_global,
  #                         ymax = 1.3 * frac_global, size = m_Deriv,
  #                         label = Letter_Derivate)) +
  #         geom_point(aes(fill = Pathway), shape = 21) +
  #         scale_size_area(max_size = 20) +
  #         theme_bw() +
  #         geom_text_repel(aes(frac_pw, frac_global, label = Letter_Derivate),
  #           size = 3,
  #           box.padding = unit(0.8, 'lines'),
  #           point.padding = unit(1.6, 'lines'),
  #           segment.color = '#cccccc',
  #           segment.size = 0.5,
  #           force = 2,
  #           max.iter = 3e3) +
  #         scale_x_log10() +
  #         facet_grid(as.formula(paste(parvec[1],"~",parvec[2]))) +
  #         scale_fill_manual(values = color_pathway) +
  #         #ggtitle(paste0('Metabolic profile of: ', var, " (", params$quant, ", ", params$analysis)) +
  #         ggtitle(paste('Metabolic profile based on (ls): ', params$quant)) +
  #         xlab('Fraction of metabolite within its pathway in (%)') +
  #         ylab('Fraction of metabolite per total metabolite content in (%)') +
  #         theme(strip.background = element_rect(fill = "white"))
  #     )
  #   }
  # }
  
  return(profile_plot_sel)
  message("Metabolic profile calculated and exported. Plot only printed in case of 2 variable condition.")
}


