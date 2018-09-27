
### Metabolic profile splitted for more than two parameter


#```{r metabolic_profile, echo=FALSE, fig.cap="Metabolic profiles according to defined groupings.", fig.width= 8, fig.height=6}
#MOD!

for (var in unique(metabolic_profile[parvec[1]])) {
  print(
    metabolic_profile %>%
      group_by(get(params$par1) == var) %>%
      ggplot(aes(frac_pw, frac_global, ymax = 1.3 * frac_global, size = m_Deriv,
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
      #facet_grid(Condition ~ Days) +
      facet_grid(as.formula(paste(parvec[1], "~", parvec[2]))) +
      #facet_grid(as.formula(paste("~", parvec[2]))) +
      scale_fill_manual(values = color_pathway) +
      ggtitle(paste0('Metabolic profile of: ', var, " (", params$quant, ", ", params$analysis, ")")) +
      #ggtitle(paste('Metabolic profile based on (ls): ', params$quant)) +
      xlab('Fraction of metabolite within its pathway in (%)') +
      ylab('Fraction of metabolite per total metabolite content in (%)') +
      theme(strip.background = element_rect(fill = "white"))
  )
}
#```  