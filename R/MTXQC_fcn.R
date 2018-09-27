# Function exporting yaml parameters (expSetup)
### Stupid Thing... don't overwrite!!!
export_MTXQCsetup = function(x, subfolder, filedef, setting) { 
  
  #Import all required file names
  nick_peakareas_file <- as.character(filedef[which(filedef$AssociatedFile == "sample_area"), "Filename"])
  nick_mz_file <- as.character(filedef[which(filedef$AssociatedFile == "mz_73"), "Filename"])
  nick_mid_file <- as.character(filedef[which(filedef$AssociatedFile == "pSIRM_se"), "Filename"])
  nick_inc_file <- as.character(filedef[which(filedef$AssociatedFile == "inc"), "Filename"])
  

  if (setting == "ExpSetup") {
    params_df = as.data.frame(do.call(rbind, x))
    params_df$arg = rownames(params_df)
    row.names(params_df) = NULL
  
    colnames(params_df) = c("Value", "Parameter")
    write.csv(params_df, paste0(subfolder, "MTXQC_params.csv"), row.names = F)
    message("MTXQC_params.csv written.")
    
    #check if not a labeling experiment
    if (x$data == "qMTX") {
      check_inc = FALSE
    } else {
      check_inc = TRUE
    }
    
    if (x$inputformat == "maui") {
        
        params_maui <- data.frame(Parameter = c("spath", "matrix", "mz", "mid", "inc_data", "inputformat",
          "intstd","alkanes", "peakchroma", "mqt", "inc"),
          Value = c(subfolder, nick_peakareas_file, nick_mz_file, nick_mid_file, nick_inc_file, "maui", 
            TRUE, TRUE, TRUE, TRUE, check_inc))
        
        write.csv(params_maui, paste0(subfolder, "Maui_params.csv"), row.names = FALSE)
        message("Maui_params.csv generated.")
    
      } else {
        message("Proceed with MTXQC_metmax in order to generate required input files.")
      }
  
  } else {
    #in case of MTXQC_Metmax
    nick_peakareas_file = x$matrix
    nick_mz_file = x$mz
    nick_mid_file = x$mid
    
    check_intst = x$intstd
    check_alkanes = x$alkanes
    check_peakchroma = x$peakchroma
    check_mqt = x$mqt
    check_inc = x$inc
    
    params_df <- data.frame(Parameter = c("spath", "matrix", "mz", "mid", "inc_data", "inputformat",
                                     "intstd","alkanes", "peakchroma", "mqt", "inc"),
                       Value = c(subfolder, nick_peakareas_file, nick_mz_file, nick_mid_file, nick_inc_file, "metmax", 
                         check_intst, check_alkanes, check_peakchroma, check_mqt, check_inc))
    
    write.csv(params_df, paste0(subfolder, setting, "_params.csv"), row.names = FALSE)
    message("Metmax_params.csv exported.")
  }
  return(params_df)
}


# Function checking for multiple data frames (all MTXQC input files)
# if they are not empty

# create a list of available data frame
inp_f <- names(which(unlist(eapply(.GlobalEnv,is.data.frame))))

check_inputfiles = function(my_list){

  #ckeck if nrow is not NULL
  my_list <- lapply(my_list, function(x) {

    
      if (empty(my_list) == TRUE) {
           message("No annotation file for additional calibration curves loaded! Please define one")
          
          #     knit_exit()
      } else {
        
        message("A")
      }
  
    return(x)
    }
  )
}



# Functions used in frame of the MTXQC to perform data transformation, absolute
# quantification and evaluation of sample peak areas

melt_call = function(dataset, col_names_id, new_col_names){
  #'alternative melting of data frames into the long format
  #'new columns are named as specified in new_col_names
  #'
  
  data_melt = melt(dataset, id = col_names_id, na.rm = T)
  colnames(data_melt) = c(col_names_id, new_col_names)
  return(data_melt)
}


qcurve_top5_rsquare = function(df, path, diff_set){
  #'Determination of calibration curves based on the
  #'ManualQuantTable. 
  #'Two versions are currently implemented:
  #'  (1) - considering different Batch_Ids (diff_set == "yes")
  #'  (2) - data only one setup (else option)
  #'
  
  
  if (diff_set == "yes") {
    
    #here we get the r-squared for each linear regression curve
    df = ddply(df, c("Metabolite_short", "Batch_Id", "Origin"), 
      transform, adj_r_squared = summary(lm(Concentration ~ ChromIntensities))$adj.r.squared)
    
    #here we get the y-intercept for each linear regression curve
    df = ddply(df, c("Metabolite_short", "Batch_Id", "Origin"), 
      transform, intercept = coefficients(lm(Concentration ~ ChromIntensities))[1])
    
    #here we get the slope for each linear regression curve
    df = ddply(df, c("Metabolite_short", "Batch_Id", "Origin"), 
      transform, slope = coefficients(lm(Concentration ~ ChromIntensities))[2])
    
    #max and min value
    df = ddply(df, c("Metabolite_short", "Batch_Id", "Origin"), 
      transform, max = max(Concentration), min = min(Concentration))
    
    #let's write these data into a file
    write.table(df, file = paste0(path, "output/quant/top5_QMQcurveInfo.csv"), row.names = F)
    message('top5_QMQcurveInfo.csv generated!')
  } else {
    
    #here we get the r-squared for each linear regression curve
    df = ddply(df, c("Metabolite_short", "Batch_Id"), transform, adj_r_squared = summary(lm(Concentration ~ ChromIntensities))$adj.r.squared)
    
    #here we get the y-intercept for each linear regression curve
    df = ddply(df, c("Metabolite_short", "Batch_Id"), transform, intercept = coefficients(lm(Concentration ~ ChromIntensities))[1])
    
    #here we get the slope for each linear regression curve
    df = ddply(df, c("Metabolite_short", "Batch_Id"), transform, slope = coefficients(lm(Concentration ~ ChromIntensities))[2])
    
    #max and min value
    df = ddply(df, c('Metabolite_short'), transform, max = max(Concentration), min = min(Concentration))
    
    #let's write these data into a file
    write.table(df, file = paste0(path, "output/quant/top5_QMQcurveInfo.csv"), row.names = F)
    message('top5_QMQcurveInfo.csv generated!')
  } 
}


islinear_nacalc = function(met, test, cc){
  #'This functions checks if the determined peak area
  #'is within, below or above the linear range of
  #'the calibration curve
  #'Functions specifies "NaCal" if no calibration curve
  #'is available
  #'If there isn't any value for absconc it reports "na"
  #'
  
  if (cc == 'yes_cal') {
    if (!all(is.na(test))) {
      met.ch = as.character(met)
      
      curmin = min(qt$ChromIntensities[qt$Metabolite == met.ch], na.rm = T)
      curmax = max(qt$ChromIntensities[qt$Metabolite == met.ch], na.rm = T)
      answer = ifelse((test >= curmin & test <= curmax), 'linear', ifelse(test < curmin, 'below','above'))
    } else {
      answer = 'na'
    }
  } else {
    answer = 'NaCal'
  }
}


absconc = function(met, area){
  #'Absolute quantification based on calibration curves
  #'equation: y = absconc = intercept + (slope * area)
  #'
  
  
  if (!is.na(area)) {
    intercept = qt$intercept[qt$Lettercode == as.character(met)][1]
    slope = qt$slope[qt$Lettercode == as.character(met)][1]
    y = intercept + (slope * area)
  } else {y = NA}
  y
}



normalisation_calc = function(d_quant, ca = 1){
  #' This functions performs the calculation of normalised quantities 
  #' considering the sum of Area normalisation and the cinnamic acid factor
  #' Following units are implemented: ul (blood), mg (tissue) and count (cell extracts)
  #' Output states respectively: pmol / ml, pmol / mg or pmol/1e+6 cells
  #' 
  #' in the absence of cinnamic acid define ca = 0, instead of ca = 1
  
  
  #norm sum of all peak areas (1)
  d_quant$sumA_Conc = d_quant$corr_absconc / d_quant$area_fac
  
  #norm over cinnamic acid factor (2)
  if (ca == 1) {
    d_quant$ca_Conc = d_quant$corr_absconc / d_quant$CA_fac
  } else {
    d_quant$ca_Conc = rep(NA, length(d_quant$Lettercode))
  }
   
  #norm over cinnamic acid factor (2)
  if (ca == 1) {
    d_quant$ca_sumA_Conc = d_quant$ca_Conc / d_quant$area_fac 
  } else {
    d_quant$ca_sumA_Conc = rep(NA, length(d_quant$Lettercode))
  }
 
  #calculating  quantities: pmio =  (1) pmol / 1x106 cells, 
  #                                 (2) pmol/ug  OR pmol/ul, 
  #                                 (3) pmol/mg OR pmol/mg
  
  d_quant$sumA_Conc_pmio = ifelse(d_quant$Unit == "ul", d_quant$sumA_Conc * 1000 / d_quant$Extract_vol, 
                                  ifelse(d_quant$Unit == "mg",
                                         d_quant$sumA_Conc * 1 / d_quant$Extract_vol, 
                                         d_quant$sumA_Conc * 1e+6 / d_quant$Extract_vol))
  
  if (ca == 1) {
    d_quant$ca_Conc_pmio = ifelse(d_quant$Unit == "ul", d_quant$ca_Conc * 1000 / d_quant$Extract_vol, 
                                ifelse(d_quant$Unit == "mg",
                                       d_quant$ca_Conc * 1 / d_quant$Extract_vol,  
                                       d_quant$ca_Conc * 1e+6 / d_quant$Extract_vol))
  } else {
    d_quant$ca_Conc_pmio = rep(NA, length(d_quant$Lettercode))
  }
    
  
  if (ca == 1) {
    d_quant$ca_sumA_Conc_pmio = ifelse(d_quant$Unit == "ul", d_quant$ca_sumA_Conc * 1000 / d_quant$Extract_vol,
                                     ifelse(d_quant$Unit == "mg",
                                            d_quant$ca_sumA_Conc * 1 / d_quant$Extract_vol, 
                                            d_quant$ca_sumA_Conc * 1e+6 / d_quant$Extract_vol))
  } else {
    d_quant$ca_sumA_Conc_pmio = rep(NA, length(d_quant$Lettercode))
  }
  
  #without any kind of normalization
  d_quant$Conc_pmio = ifelse(d_quant$Unit == "ul", d_quant$corr_absconc * 1000 / d_quant$Extract_vol,
                             ifelse(d_quant$Unit == "mg", 
                                    d_quant$corr_absconc * 1 / d_quant$Extract_vol,
                                    d_quant$corr_absconc * 1e+6 / d_quant$Extract_vol)) 
  
  
  #calculation of Mole (M) for UNITS: UL
  if (ca == 1) {
    d_quant$Conc_microM = ifelse(d_quant$Unit == "ul", 
      d_quant$ca_Conc_pmio * d_quant$Extract_vol * 1000 / (1000 * 1000), "")
  } else {
    d_quant$Conc_microM = rep(NA, length(d_quant$Lettercode))
    #d_quant$Conc_microM = d_quant$ifelse(d_quant$Unit == "ul", 
     #                                    d_quant$Conc_pmio * d_quant$Extract_vol * 1000 / (1000 * 1000), "")
  }
  
  if (nrow(d_quant) == 0) {
    message("Empty data frame! Check column names for matching annotation (sample_extracts, data)!")
  }
  
  return(d_quant)
}


create_batchid = function(df){
  #' This function creates the batch id based on a dataframe
  #' containing the column "File" or "file"
  #' It is used for the automated merge of ManualQuantTables
  #' to run combined data analysis for multiple MAUI projects
  #' 
  #' input: data frame containing the column File with appropiate identifier
  #' , e.g., p17256xy for Phineas, 2017, 256th day measured by xy
  #' output: vector containing unique batch identifiers
  #' 
  
  #create Batch_ID
  temp = sapply(strsplit(as.character(df$File), "\\_"), "[[", 1)
  
  #extract unique elements
  temp_un = unique(temp)
  
  return(temp_un)
  
}



metabolic_profile = function(df, groups, value_mean) {
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
            scale_fill_manual(values = color_pathway) +
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
          scale_fill_manual(values = color_pathway) +
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


check_emptyfile <-  function(dataframe, essential = FALSE, tfile = NULL) {
  
  if (nrow(ann) != 0) {
    message("File imported! ", tfile)
  } else {
    message("Empty file detected!")
    
    if (essential == TRUE) {
      message("Essential file - please check. Processing canceled! ", tfile)
      knitr::knit_exit()
    }
  }
}


file_shaping = function(dataframe, shape = "long", file_annotation, type, mode = NA, complete_ann = FALSE) {
### Function: File-shaping
#' including:
#' melting of wide to long format
#' checkup File instead of file
#' adding annotation - File,Type
#' add Batch_Id based on File name
#' WHATS WITH METABOLITE_MANUAL
  

  #check File
  colnames(dataframe)[grepl("file", colnames(dataframe))] <- "File" 
  
  #metmax shaping
  test_row = "row.load" %in% names(dataframe)
  
  if (shape == "wide") {
  
    if (test_row == TRUE) {
      idx_ri = grep("ri", colnames(dataframe))
      idx_row = grep("row.load", colnames(dataframe))
      df = dataframe[,c(-idx_row, -idx_ri)]
      dataframe = reshape2::melt(df, id.vars = c("name", "mass"), value.name = "PeakArea",
        variable.name = "File", na.rm = TRUE)
      
      colnames(dataframe)[grepl("name", colnames(dataframe))] <- "Metabolite_manual"
      colnames(dataframe)[grepl("mass", colnames(dataframe))] <- "QuantMasses"
    }
    
    if (test_row == FALSE) {
      if (mode != "inc") {
        data_proc = dataframe
        dataframe = reshape2::melt(data_proc, id.vars = c("Metabolite", "QuantMasses"), 
                                   value.name = "PeakArea", variable.name = "File", na.rm = TRUE)
      } else
      {
        data_proc = dataframe
        dataframe = reshape2::melt(data_proc, id.vars = c("Metabolite"), 
                                   value.name = "PeakArea", variable.name = "File", na.rm = TRUE)
      }
    }
  }
  
  #add Batch_Id
  dataframe$Batch_Id  = sapply(strsplit(as.character(dataframe$File), "\\_"), "[[", 1)
  
  #merge with annotation to exclude project-unrelated data
  if (complete_ann == FALSE) {
    dataframe = merge(dataframe, file_annotation[,c("File", "Type")])
  }
  
  #
  if (complete_ann == TRUE) {
    #written for the generation of MQTable
    dataframe = merge(dataframe, file_annotation)
    colnames(dataframe)[grep('PeakArea', colnames(dataframe))] <- 'ChromIntensities'
  }
  
  if (type == "sample") {
    dataframe = subset(dataframe, dataframe$Type == as.character(type))
  }
  
  return(dataframe)

}

gc_metric_calc = function(dataframe, dataframe2 = NULL, title) {
#' Function calculating the qc metrices for GC-performance
#' title = alkanes, mz73, cinacid, sumArea  
  
  #ALKANE
  if (title == "alkanes") {
    df_stat  =  ddply(dataframe, c("File", "Batch_Id"), summarise,
                       N = length(intensity),
                       val = mean(intensity, na.rm = T),
                       val_sd = sd(intensity))
    
    #stats batch-wise
    batch_stats = ddply(df_stat, c("Batch_Id"), summarise, 
                        mean_batch = mean(val),
                        sd_batch = sd(val_sd))
  }
  
  #CINACID
  if (title == "cinacid") {
    df_stat <-  ddply(dataframe, c("Batch_Id"), transform, 
          CA_fac = PeakArea / mean(PeakArea))
    
    #Quality check - should be within a factors range of 0.65 and 1.45
    df_stat$CA_eval <-  ifelse(df_stat$CA_fac >= 0.65 & df_stat$CA_fac <= 1.45, 'within', 
                           ifelse(df_stat$CA_fac < 0.65, 'below','above'))
    
    #Export CA-normalization information
    write.csv(df_stat, paste0(path_setup,'output/gc/Cinacid_normfactors.csv'), row.names = F)
    
    #calculate cinnamic acid-factor
    batch_stats <-  ddply(df_stat, c("Batch_Id"), transform,
                        n_batch = length(PeakArea),
                        mean_batch = mean(PeakArea, na.rm = TRUE),
                        sd_batch = sd(PeakArea, na.rm = TRUE))
    
    write.csv(batch_stats, paste0(path_setup, 'output/gc/Cinacid_stats.csv'), row.names = F)
    }
  
  #SUMOFAREA
  if (title == 'sumofarea') {
    
    #in order to exclude multiple assignments of one file (T=0 or control files)
    df_unique = unique(dataframe)
    
    area_stats = ddply(df_unique, c("File", "Batch_Id"), summarise,
                    n_area = length(PeakArea),
                    sum_area = sum(PeakArea))
    
    batch_stats <- ddply(area_stats, c('Batch_Id'), transform,
                    n_total = max(n_area),
                    mean_batch = mean(sum_area, na.rm = TRUE),
                    sd_batch = sd(sum_area, na.rm = TRUE))
    
    batch_stats$area_fac = batch_stats$sum_area / batch_stats$mean_batch
    batch_stats = ddply(batch_stats, c('Batch_Id'), transform,
                        n_50 = n_total * 50 / 100)
    
    write.csv(unique(batch_stats[,c('Batch_Id', 'n_50')]), paste0(path_setup, 'output/gc/Min_Annotation.csv'), row.names = FALSE)
    message("Files with less than 50% of max(N) should be excluded from SumofArea normalisation.")
    
    #Export
    write.csv(batch_stats, paste0(path_setup, 'output/gc/SumArea_stats.csv'), row.names = F)
    
  }
  
  #DERIVATISATION
  if (title == 'mz73') {
    
    #Calculate N, mean and sd for target mass per file
    data73_stat = ddply(dataframe, c("File", "Batch_Id"), summarise, 
                      mean_73 = mean(targetmass_intensity, na.rm = TRUE),
                      sd_73 = sd(targetmass_intensity, na.rm = TRUE))         
    
    #Calculate total peak area per file
    tp_stat = ddply(dataframe2, c("File","Batch_Id"), summarise,
                      n_peaks = length(Area),
                      sum_area = sum(Area, na.rm = TRUE))
    
    #Combine 
    batch_prep = merge(data73_stat, tp_stat)
    batch_prep$ratio_total = batch_prep$mean_73 / batch_prep$sum_area
    
    write.csv(batch_prep, paste0(path_setup, 'output/gc/mz73_data.csv'), row.names = F)
    
    #Calculate mean and sd of the ratio per Batch_Id
    batch_stats = ddply(batch_prep, c("Batch_Id"), transform, 
                        mean_batch = mean(ratio_total),
                        sd_batch = sd(ratio_total))
  }
  
  #calc qc_metric
  batch_stats$qc_metric = 1 -  batch_stats$sd_batch / batch_stats$mean_batch
  qcm = unique(batch_stats[,c("Batch_Id", "qc_metric")])
  qcm$title = rep(title, length(qcm$Batch_Id))
  
  #export
  write.csv(qcm, paste0(path_setup,'output/gc/qcmetric_', title, '.csv'), row.names = F)
  message(paste0("QC-metric succesfully exported: ", title))
  return(qcm)
  
}


extract_addQ_annotation = function(annotation_file, phrase) {
  
  sel_q_ann = annotation_file[grepl(phrase, annotation_file$Type),]
  
  test_row = "Setup" %in% names(annotation_file)
  
  if (test_row == TRUE) {
    sel_q_ann = sel_q_ann[,c('File','Type','Setup')]
  } else {
    sel_q_ann = sel_q_ann[,c('File','Type')]
  }
  
  #Extract Dilution information
  sel_q_ann$Dilution_Step = sapply(strsplit(as.character(sel_q_ann$Type), "_"), "[", 2)
  sel_q_ann$Dilution = 1 / as.numeric(sel_q_ann$Dilution_Step)
  
  return(sel_q_ann)
}


create_manualquanttable = function(cal_dataframe, q1_values, met_translation, plot = FALSE) {
  
  #chromintensities
  cal_dataframe = cal_dataframe[,c("Batch_Id", "Metabolite", "Dilution", "ChromIntensities")]
  cal_dataframe = merge(cal_dataframe, met_translation[,c("Metabolite", "Lettercode")])
    
  #mqt_init  
  mqt_init = merge(q1_values, cal_dataframe)
  mqt_init = subset(mqt_init, mqt_init$ChromIntensities != 99)  
  
  #calculate Concentration Quant 1:1
  quant_idx = grep("Quant", colnames(mqt_init))
  mqt_init$Concentration = mqt_init[,quant_idx] * mqt_init$Dilution
  
  #export -> distinguish addQ and standardQ
  #var_quantversion = unlist(strsplit(colnames(mqt_init[q_idx]), split = "_"))[2]
  var_quantversion = unlist(strsplit(colnames(mqt_init[quant_idx]), split = "_"))[2]
  
 if (var_quantversion == 'ext') {
   #add flag-tag for quant-values origin
   mqt_init$Origin = rep("Qadd", length(mqt_init$Lettercode))
   
   write.csv(mqt_init, paste0(path_setup, set_input, "add_quant/ManualQuantTable_additionalQ.csv"), 
     row.names = F)
   
   message(paste0("ManualQuantTable for additional calibration curves has been generated. Quant1-values: ", 
     colnames(mqt_init[quant_idx])))
 } else {
   #add flag-tag for quant-values origin
   mqt_init$Origin = rep("Qstd", length(mqt_init$Lettercode))
   
   write.csv(mqt_init, paste0(path_setup, set_input, "quant/ManualQuantTable_calc", "_", 
                colnames(mqt_init[q_idx]),".csv"), row.names = F)
   
   message(paste0("ManualQuantTable for standard calibration curves has been generated. Quant1_", 
                var_quantversion))
 }
  
  if (plot == TRUE) {
    for (var in unique(mqt_init$Batch_Id)) {
      #ggplot(mqt_init, aes(Concentration, ChromIntensities)) +
      print(ggplot(subset(mqt_init, mqt_init$Batch_Id == var), aes(Concentration, ChromIntensities)) +
          geom_point() +
          geom_smooth(method = "lm", se = FALSE) +
          theme_bw() +
          ggtitle(paste0("Calibration curves: ", var)) +
          facet_wrap( ~ Lettercode, scales = "free") +
          theme(strip.text = element_text(size = 8)) #,
        #axis.text.y = element_blank(),
        #axis.text.x = element_blank())
      )
    }
  }

  return(mqt_init)
}

quant_metric_calc_new <- function(dataframe) {
  
  df_small <- dataframe[,c('Metabolite','Lettercode',"Batch_Id", "Origin",
    'adj_r_squared','intercept','slope')]
  
  batch_stats <- ddply(df_small, c("Lettercode", "Batch_Id", "Origin"), transform,
    Frac_calcurve = length(adj_r_squared) / 8)
  
  #remove duplicates
  batch_stats = unique(batch_stats, by = c("Metabolite", "Origin"))
  write.csv(batch_stats, paste0(path_setup,"output/quant/top5_CalibrationInfo_unique.csv"), row.names = F)
  
  #clean-up header
  colnames(batch_stats)[grepl("adj_r_squared", colnames(batch_stats))] <- "R2"
  
  #create table for plots
  qt_plot = melt(batch_stats, id.vars = c('Metabolite','Lettercode', "Batch_Id", "Origin"), 
    variable.name = 'Parameter', value.name =  'Par_value')
  
  qt_plot = subset(qt_plot, qt_plot$Parameter == 'R2' | 
      qt_plot$Parameter == 'Frac_calcurve')
  
  return(qt_plot)
}

quant_metric_calc <- function(dataframe) {
  
  df_small <- dataframe[,c('Metabolite','Lettercode',"Batch_Id",
    'adj_r_squared','intercept','slope')]
  
  batch_stats <- ddply(df_small, c("Lettercode", "Batch_Id"), transform,
    Frac_calcurve = length(adj_r_squared) / 8)
  
  #remove duplicates
  batch_stats = unique(batch_stats, by = 'Metabolite')
  write.csv(batch_stats, paste0(path_setup,"output/quant/top5_CalibrationInfo_unique.csv"), row.names = F)
  
  #clean-up header
  colnames(batch_stats)[grepl("adj_r_squared", colnames(batch_stats))] <- "R2"
  
  #create table for plots
  qt_plot = melt(batch_stats, id.vars = c('Metabolite','Lettercode', "Batch_Id"), 
    variable.name = 'Parameter', value.name =  'Par_value')
  
  qt_plot = subset(qt_plot, qt_plot$Parameter == 'R2' | 
      qt_plot$Parameter == 'Frac_calcurve')
  
  return(qt_plot)
}


eval_extractionfactor <- function(params) {
  
  #quantification factor determination
  #1/3 = 500 ul of 1500 ul quant mix polar phase dried
  quant_fullvol = as.numeric(1500)
  
  qsingle_idx = as.numeric(as.character((params[which(params$Parameter == "quant_vol"), "Value"])))
  quant_fac  = qsingle_idx / quant_fullvol  
  
  #sample factor determination
  #1: no technical backups 2: split into two technical backups  
  sample_fac = as.numeric(as.character(params[which(params$Parameter == "backups"), "Value"]))
  
  #combined = extraction factor
  extr_fac = quant_fac * sample_fac
  
  message('The quantification factor for that experimental setup: ', quant_fac)
  message('The sample factor for that experimental setup: ', sample_fac)
  message('The extraction factor for that experimental setup: ', extr_fac)
  
  return(extr_fac)
}


evaluate_qt_lin <- function(dataframe) {
  
  df_calcheck = ddply(dataframe, c('Metabolite','File', 'Batch_Id', 'Origin'), transform,
    islinear = islinear_nacalc(Metabolite, PeakArea, calc_curve))
  
  df_calcheck$islinear = factor(df_calcheck$islinear, 
    levels  =  c('below', 'linear', 'above', 'NA', 'NaCal'))
  
  write.csv(df_calcheck, paste0(path_setup, set_output, 'quant/calcheck_linearity.csv'), row.names = FALSE)
  
  #Calculate fraction of each level across Batches
  fraction_lincheck = ddply(df_calcheck, c("Lettercode", "islinear", "Batch_Id", "Origin"), 
    summarize, 
    count = length(PeakArea))
  
  fraction_lincheck = ddply(fraction_lincheck, c("Lettercode", "Batch_Id", "Origin"), transform,
    sum_lin = sum(count))
  
  fraction_lincheck$prop = fraction_lincheck$count / fraction_lincheck$sum_lin
  
  #Export list
  write.csv(fraction_lincheck, paste0(path_setup, set_output,
    "quant/pTop5_Calibration_Samples_lincheck.csv"), row.names = FALSE)
  
  message("Position of data points regarding calibration curves evaluated.")
  return(fraction_lincheck)
}



mid_metric_calc <- function(dataframe) {
  
  df_nascore  =  ddply(dataframe, c("File", "Lettercode", "Batch_Id"), summarise,
    na_SPA  =  sum(is.na(SamplePeakArea)),
    na_BPA  =  sum(is.na(BackupPeakArea)),
    #nb.mid = length(BackupPeakArea),
    na_frac  =  (na_SPA - na_BPA) / length(BackupPeakArea)) 
  
  #negative values in na.frac replaced with Zeros
  df_nascore[df_nascore <=  0]  =  0
  
  #round na.frac, define as factor
  df_nascore$na_frac  =  round(as.numeric(df_nascore$na_frac), 1)
  write.csv(df_nascore, paste0(path_setup, set_output, 'inc/SE_calculation_NAscore.csv'),
    row.names  =  F)
  
  
  df_nascore$na_frac_r  =  as.factor(df_nascore$na_frac)
  
  #statistics
  na_fraction_calc  =  ddply(df_nascore, c("Lettercode", "Batch_Id" ,"na_frac_r"), 
                      summarise,
                      N  =  length(File)) #,
                      #fracr_prop  =  round(as.numeric(N/length(unique(df_nascore$File))), 1))
  
  na_fraction_calc$fracr_prop = round(as.numeric(na_fraction_calc$N/length(unique(na_fraction_calc$File))), 1)
  
  na_frac_zero  =  subset(na_fraction_calc, na_fraction_calc$na_frac_r  ==  0)
  na_frac  =  na_frac_zero[,c('Lettercode', 'Batch_Id','fracr_prop')]
  write.csv(na_frac, paste0(path_setup, set_output, 'inc/SE_NAscore.csv'),
    row.names  =  F)
  
  return(df_nascore)

}

mid_metric_calc2 <- function(dataframe1, dataframe2) {
  
  SPA_low3A = ddply(dataframe1, c("Lettercode", "File"), subset, 
    SamplePeakArea %in% sort(SamplePeakArea, FALSE, na.last  =  FALSE)[1:3])
  
  spa_low3A_stat  =  ddply(SPA_low3A, c("Lettercode", "Batch_Id" ,"File"), summarise,
    sum_spa  =  sum(SamplePeakArea, na.rm = T))
  
  #BackupPeakArea
  BPA_sub  =  dataframe2[,c("Lettercode","File", "Batch_Id" ,"BackupPeakArea")]
  BPA_low3A  =  ddply(BPA_sub, c("Lettercode","File", "Batch_Id"), subset, 
    BackupPeakArea %in% sort(BackupPeakArea, FALSE, na.last = FALSE)[1:3])
  
  bpa_low3A_stat  =  ddply(BPA_low3A, c("Lettercode", "Batch_Id" ,"File"), summarise,
    sum_bpa  =  sum(BackupPeakArea, na.rm = T))
  
  #combine both
  low3A = merge(spa_low3A_stat, bpa_low3A_stat)
  #calculate ratio
  low3A$low3a_ratio = low3A$sum_spa / low3A$sum_bpa
  
  #sum of samplelw3 > backup3low
  low3A$rel_sb = ifelse(low3A$sum_spa > low3A$sum_bpa, 'higher','below')
  
  #low3A_ratio CLASSIFICATION
  low3A$val_score = ifelse(low3A$low3a_ratio < 2, ifelse(low3A$low3a_ratio == 0, 'bad', 'ok'), 'confident')
  
  #REQUIRED: if rel_sb  ==  Higher AND val_score ok / confident -> count as GOOD MID
  low3A$count_score = ifelse((low3A$rel_sb == 'below' & low3A$val_score == 'bad'),'lowQ','goodQ')
  
  #export
  write.csv(low3A, paste0(path_setup,'output/inc/SE_classification.csv'), row.names = F)
  
  #stats:
  low3A_stats  =  ddply(low3A, c("Lettercode", "Batch_Id" ,"count_score"), summarise,
    N_count = length(File))
  
  #calculate proportion
  low3A_stats  =  ddply(low3A_stats, c("Lettercode", "Batch_Id"), transform, 
    sum_files  =  sum(N_count))
  
  low3A_stats$prop  =  low3A_stats$N_count/low3A_stats$sum_files
  
  #export
  write.csv(low3A_stats, paste0(path_setup,'output/inc/SE_validation.csv'), row.names = F)
  
  return(low3A_stats)
}

extract_standards_export <- function(dataframe, met_names = con_se, ann_file = ann) {
 
   #dataframe = df_peakareas 
   #met_names = con_se
   #ann_file = ann
  
  
  #Extract internal standard
  is_subs = subset(met_names, met_names$Standards == "InternalStandard")
  internalstd = unique(is_subs$Lettercode)
  
  #Extract measures peak areas
  df_peak = merge(dataframe, met_names[,c("Metabolite_manual", "Metabolite","Lettercode")])
  df_standard = subset(df_peak, Lettercode %in% internalstd)
  
    if (is.null(dim(df_standard))) {
      
      message("Defined internal standard: ", internalstd)
      message("Peak areas detected for internal standard in peak area matrix.")
      
      #Selected columns for wide format
      df_sel = df_standard[,c("Metabolite", "QuantMasses" ,"File", "PeakArea")]
      df_wide = reshape2::dcast(df_sel, Metabolite + QuantMasses ~ File, value.var = "PeakArea")
      
      #Export
      write.csv(df_wide, paste0(path_setup, set_input, "gc/InternalStandard.csv"), row.names = FALSE)
      
      message("Data file for internal standard generated and exported to: input/gc/InternalStandard.csv")
      
    } else {
      
      message("Defined internal standard: ", internalstd)
      message("WARNING: No peak areas detected for internal standard in peak area matrix.")
      
    }

  
}

extraction_alkanes_export <- function(dataframe, phrase = "Alk", met_names = con_se, annotation_file = ann) {
  
  alk_set = subset(con_se, con_se$Standards == phrase)
  
  df_alk = merge(dataframe, alk_set)
  if (nrow(df_alk) == 0) {
    
    face_df = data.frame(file = as.character(), metabolite = as.character(), 
                          intensity = as.numeric(), quantMasses = as.character() )
    
    write.csv(face_df, paste0(path_setup, set_input, "gc/Alcane_Intensities.csv"), row.names = FALSE)
    message("No peak areas for alkanes in input table detectable!")
  } else {
    alk_table = df_alk[, c("File", "Lettercode", "PeakArea")]
    alk_table$quantMasses = rep("71;85;99", length(alk_table$File))
    
    colnames(alk_table) = c("file", "metabolite", "intensity", "quantMasses")
    write.csv(alk_table, paste0(path_setup, set_input, 
      "gc/Alcane_Intensities.csv"), row.names = FALSE)
    
    message("Alkane intensities have been exported to input/gc/Alcane_Intensities.csv")
  }
  
}


calculate_isotope_incorporation = function(dataframe, backups, mass_li, manval = FALSE) {
  
  if (manval == TRUE) {
    #clean column names
    colnames(dataframe)[grepl("ManVal_Intensity", colnames(dataframe))] <- "MID_Intensity"
  }

  #calculate MIDs
  dataframe = ddply(dataframe, c("File", "Metabolite"), transform, 
    SampleMID = MID_Intensity / sum(MID_Intensity))
  
  #calculate cor.natural = ref_mi / ref_m0 x sam_m0
  colnames(backups)[grepl("Mass.m.z.", colnames(backups))] <- "Mass_mz"
  
  #merge
  data_li = merge(dataframe, mass_li)
  data_li_export = merge(data_li, backups)
  
  #select specific columns
  data_li_sel = data_li_export[,c( 'File', 'Metabolite','Mass_mz', 'BackupMID', 
    'SampleMID', 'LI_MID')]
  
  #ratio_ref = ratio ref_mi / ref_m0
  ref_mid = data_li_sel[,c('File', 'Metabolite','Mass_mz', 'BackupMID','LI_MID')]
  ref_wide = dcast(ref_mid, File+Metabolite ~ LI_MID, value.var = 'BackupMID')
  
  ref_wide$ratio_ref = ref_wide$minc / ref_wide$m0
  
  ref_wide_sel = ref_wide[,c('File', 'Metabolite', 'ratio_ref')]
  
  #transform sample data
  sample_mid = data_li_sel[,c('File', 'Metabolite','Mass_mz', 'SampleMID','LI_MID')]
  sample_wide = dcast(sample_mid, File + Metabolite ~ LI_MID, value.var = 'SampleMID')
  sample_wide_s = merge(sample_wide, ref_wide_sel)
  
  #Calculate correction for natural MID: cor_natural = ratio * sample_m0
  sample_wide_s$cor_natural = sample_wide_s$ratio_ref * sample_wide_s$m0
  
  #Correct sample minc for natural MID = sample_minc - cor_natural
  sample_wide_s$cor_sample = sample_wide_s$minc - sample_wide_s$cor_natural
  
  #calculate LI
  ref_m0 = ref_wide[,c('File', 'Metabolite', 'm0')]
  colnames(ref_m0) = c('File', 'Metabolite', 'ref_m0')
  comb_mid = merge(sample_wide_s, ref_m0)
  comb_mid$LI_raw = comb_mid$cor_sample / (comb_mid$cor_sample + comb_mid$ref_m0)
  
  #replace negative with 0
  comb_mid$LI = ifelse(comb_mid$LI <= 0.0000, 0, comb_mid$LI)
  
  #LI only
  data_li = comb_mid[,c('File', 'Metabolite', 'LI_raw' ,'LI')]
  
  #wide format using cleaned LI
  data_inc_new = dcast(data_li, Metabolite ~ File, value.var = 'LI')	
  
  if (manval == TRUE) {
    #assumin maui-derived input after manual validation
    write.csv(data_li, paste0(path_setup, set_output, set_val,
                              'inc/LI_long_calculated.csv'), row.names = F)

    filen = paste0(set_val,'inc/LI_long_calculated.csv')
    message("Updated 13C-incorporation has been calculated and saved in: ", filen)
    
  } else {
    #assuming metmax-derived input
    write.csv(data_li, paste0(path_setup, set_input,
                              'metmax/Metmax-LI_long_calculated.csv'), row.names = F)
    
    write.csv(data_inc_new, paste0(path_setup, set_input, 'inc/DataMatrix.csv'), row.names = F)
    
    message("The Metmax-exported MIDS have been converted.")
    message("Determined 13C-incorporation has been saved: ", paste0(set_input, 'inc/DataMatrix.csv'))
  }
  #return(data_inc_new)
  return(data_li_export)
}

MID_export <- function(dataframe, backups = backup_mids) {
  
  #clean-up2 - restarting from the first part at line #
  table_conv = dataframe[,c( "File","Metabolite", "Mass_mz", "MID_Intensity", "SampleMID")]
  colnames(table_conv)[grepl("MID_Intensity", colnames(table_conv))] <- "SamplePeakArea"
  
  #Add backupMIDs
  colnames(backups)[grepl("Mass.m.z.", colnames(backups))] <- "Mass_mz"
  
  table_conv2 = merge(table_conv, backup_mids)
  table_conv2$UsedMID = rep("backup", length(table_conv2$File))
  
  #export
  write.csv(table_conv2, paste0(path_setup, set_input, "inc/pSIRM_SpectraData.csv"), row.names = F)
  
  message("Metmax-derived MIDs have been transformed into classical MTXQC input format.") 
  message("File saved in input/inc/pSIRM_SpectraData.csv.")

}


check_parvector <- function(params = params, annotation_file = ann) {
  
  pv = as.character()
  parvec = c(params$par1, params$par2, params$par3, params$par4)
  
  #Remove empty fields
  for (n in 1:4) {
    if (parvec[n] != "") {
      pv[n] = parvec[n]
    }
  }
  
  if (unique(pv %in% colnames(ann)) == TRUE) {
    message("All defined parameters detected!")
  } else {
    message("Please check your parameter selection!")
    knitr::knit_exit()
  }
  return(pv)
}



create_mauiparams <- function(paths) {
  
  def_maui = data.frame(Parameter = c("inststd", "alkanes", "peakchrom", "mqt", "inc"),
                        Value = c(TRUE, TRUE, TRUE, TRUE, FALSE))
  
  write.csv(def_maui, paste0(paths, "Maui_params.csv"), row.names = FALSE)
  return(def_maui)
  
}

MTXQCp1_checkinput <- function(data_extracts = data_extracts, ann = ann) {
  
  er <-  0
  
  ## Check sample extract: $Extract_vol, $Unit
  if (any(names(data_extracts) == 'Extract_vol' & any(names(data_extracts) == "Unit")) == TRUE) {
    message("Correct column names in file sample_extracts.csv")
  } else {
    er <- er + 1
    message("Sample extracts: missing columns - Extract_Vol, Unit")
  }
  
  ## Check sample annotation: $File, $Type
  
  if (any(names(ann) == 'File' & any(names(ann) == "Type")) == TRUE) {
    message("Correct column names in sample annotation")
  } else {
    er <-  er + 1
    message("Sample annotation: missing columns - File, Type")
  }
  
  ## Total error count
  if (er != 0) {
    message("Please check your input files regarding required input columns!")
    message("Input file check ups detected: ", er, " error(s)! ")
    knitr::knit_exit()
  } else {
    message("Input files checked!")
  }
  
  return(er)
}
