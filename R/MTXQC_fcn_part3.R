### MTXQC-functions: ManualValidation

transform_quant <- function(dataframe) {
  
  #### Quantities #########
  dataframe_trans = melt(dataframe, id = c("Metabolite", "QuantMasses"), 
                         value.name = "PeakAreas", variable.name = "File")
  
  #add a column to confirm check-up
  dataframe_trans$ManVal_PeakArea = rep("", length(dataframe_trans$Metabolite))
  
  #rearrangement of columns
  dataframe_file = dataframe_trans[,c("File", "Metabolite", "QuantMasses", 
                                      "PeakAreas", "ManVal_PeakArea")]
  
  #export
  write.csv(dataframe_file, paste0(path_setup, set_output, set_val, 
                                   "ManVal_PeakAreas.csv"), row.names = F)
  
  message("Transformation of quantMassAreaMatrix.csv is done! Check ManVal_PeakAreas.csv")
}



transform_inc <- function(df_se_val, df_mid, conversion_table = con_se) {
  
  #SpectraExport Quality
  df_se_val$ManVal_check = rep("", length(df_se_val$Lettercode))
  
  #### line modified:se_val_low = subset(df_se_val, df_se_val$count_score == "lowQ")
  #### Export all MIDs without preselection of lowQ only
  se_val_low = subset(df_se_val, df_se_val$count_score != "")
  
  se_low = se_val_low[,c("File","Lettercode","count_score", "ManVal_check")]
  
  #Export
  write.csv(se_low, paste0(path_setup, set_output, set_val, 
                           "MID_validation_lowQ_forCheck.csv"), row.names = FALSE)
  
  #SpectraExport
  colnames(df_mid)[3] = c('Mass_mz')
  data_inc_table = merge(df_mid, conversion_table)
  data_inc_table_merge = merge(data_inc_table, se_low)
  
  data_inc_table = data_inc_table_merge[, c("File","Lettercode" ,"Metabolite", 
                                            "Mass_mz", "SamplePeakArea", "BackupPeakArea", "BackupMID")]
  
  data_inc_table$ManVal_Intensity = rep("", length(data_inc_table$Metabolite))
  
  data_inc_export = arrange(data_inc_table, File, Lettercode, Mass_mz) 
  
  #Export
  write.csv(data_inc_export, paste0(path_setup, set_output, set_val, 
                                    "MID_validation_lowQ_values.csv"), row.names = F)
  
  message('Files for manual MID evaluation are generated!')
}


evaluate_peakareas <-  function(dataframe) {
  
  ## check if manual validated peak areas present or copy PeakAreas
  dataframe$ManVal_PeakArea = as.numeric(dataframe$ManVal_PeakArea)
  
  dataframe$Comb_PeakArea = ifelse(!is.na(dataframe$ManVal_PeakArea), 
                                   dataframe$ManVal_PeakArea,
                                   ifelse(!is.na(dataframe$PeakAreas), dataframe$PeakAreas, NA))
  
  dataframe$check_diff = ifelse(!is.na(dataframe$PeakArea), dataframe$ManVal_PeakArea - dataframe$PeakAreas, dataframe$ManVal_PeakArea)
  
  
  quant_nonzero = subset(dataframe, !is.na(dataframe$check_diff))
  
  #in case of empty dataframe: quant_nonzero
  if (nrow(quant_nonzero) == 0) {
    
    message("No modified peak areas detected. No data integration and export performed!")
  
  } else {
    
    nb_mod_peaks = length(quant_nonzero$check_diff)
    message('Number of modified peak areas: ', nb_mod_peaks)
    
    ##first: plot tendency of modified peak areas for each intermediate
    quant_nonzero = merge(quant_nonzero, con_se[,c("Lettercode", "Metabolite")])
    
    print(ggplot(quant_nonzero, aes(check_diff, fill = Lettercode)) +
            geom_histogram(size = .3, color = 'black') +
            theme_bw() +
            ggtitle('Check_diff = ManVal-PeakArea - Original-PeakArea') +
            theme(legend.position = "bottom"))
    
    ##second check: count modified peak areas per metabolite
    man_stat = ddply(quant_nonzero, c("Lettercode"), summarise,
                     n_mod = length(ManVal_PeakArea))
    
    print(ggplot(man_stat, aes(Lettercode, n_mod)) +
            geom_point() +
            coord_flip() +
            theme_bw() +
            ggtitle('Number of validated peak areas / intermediate') +
            theme(legend.position = "bottom") +
            ylab("Nb. of manually validated peaks") +
            xlab("Metabolite"))
    
    ## convert long into wide format and export into:
    ## output/quant/quantAreaMattrix_manVal.csv
    
    
    dataframe_wide = dcast(dataframe, Metabolite + QuantMasses ~ File, 
                           value.var = 'Comb_PeakArea')
    
      ##ERROR: aggregation function missing
        
        test_long = reshape2::melt(dataframe_wide, c("Metabolite", "QuantMasses"), na.rm = TRUE)
    
        test = ddply(test_long, c("Metabolite"), summarise, 
                     mean_val = mean(value))
        
        if (mean(test$mean_val) <= 250) {
          
          ## WARNING aggregation function missing
          message("WARNING: Be cautious if you get a warning like: aggregation function missing.")
          message("WARNING: If so - check your MassAreasMatrix_ManVal.csv - file.")
          message("WARNING: Most probably you have multiple values for one file and metabolite! dcast failes!")
          
          write.csv(test, paste0(path_setup, set_output, "ManVal-Aggregationfailure.csv"), row.names = FALSE)
          
          message("Have a look at ManVal-Aggregationfailure.csv in output-folder!")
          
        } else {
          message("Integration of validated peak areas seems to be fine!")
        }
    
        
    
    #Export
    write.csv(dataframe_wide, paste0(path_setup, set_input, 
                                     'quant/quantMassAreasMatrix_manVal.csv'), row.names = FALSE)
    
    file_name = paste0(set_input, "quant/MassAreasMatrix_ManVal.csv")
    
    message("Manual validated peak areas have been merged original data") 
    message("Updated peak areas saved in: ", file_name) 
    
    
   
  }
  
}


evaluate_modified_mids <-  function(df_check) {
  
  #how many MIDs were corrected
  df_stat = count(subset(df_check, df_check$ManVal_check == 'x'), 
                  vars = 'Lettercode')
  
  print(knitr::kable(df_stat, format = "markdown", 
                     caption = "Frequency of validated isotope incorporation values."))
  
  #extract corrected MIDs/Files
  cor = subset(df_check, df_check$ManVal_check == 'x')
  cor_lettercode_file = unique(cor[,c('File','Lettercode' ,'ManVal_check')])
  
  return(cor_lettercode_file)
  write.csv(cor_lettercode_file, paste0(path_setup, set_output, set_val, 
                                        'inc/Corrected_MIDs_summary.csv'), row.names = F)
}


integrate_manVal_MIDs <- function(dataframe, corrected_mids, original_mids) {
  
  temp = merge(dataframe, corrected_mids, all.x = TRUE)
  
  #corrected MIDS
  mid_corr = subset(temp, temp$ManVal_check == "x")
  mid_corr_sel = mid_corr[,c("File", "Metabolite", "Mass_mz", "ManVal_Intensity")]
  
  #merge with original MIDs
  colnames(original_mids)[grepl("Mass.m.z.", colnames(original_mids))] <- "Mass_mz"
  colnames(original_mids)[grepl("SamplePeakArea", colnames(original_mids))] <- "SamplePeakArea_orig"
  colnames(original_mids)[grepl("SampleMID", colnames(original_mids))] <- "SampleMID_orig"
  
  mid_merge = merge(original_mids, mid_corr_sel, all.x = TRUE)
  
  mid_merge$SamplePeakArea = ifelse(is.na(mid_merge$ManVal_Intensity), mid_merge$SamplePeakArea_orig,
                                    mid_merge$ManVal_Intensity)
  
  mid_merge = ddply(mid_merge, c("File", "Metabolite"), transform, 
                     SampleMID = SamplePeakArea / sum(SamplePeakArea))
  
  write.csv(mid_merge, paste0(path_setup, set_output, set_val, "inc/Fused_MIDs.csv"), row.names = FALSE)
  
  idx_ori1 = which(grepl("_orig", colnames(mid_merge)))
  fused_mids = mid_merge[,c(-idx_ori1)]
  
  write.csv(fused_mids, paste0(path_setup, set_input, "inc/pSIRM_SpectraData_manVal.csv"), row.names = FALSE)
  message("Manually validated MIDs have been incorporated and saved in: ", paste0(set_input, "inc/pSIRM_SpectraData_manVal.csv"))
  return(fused_mids)   
}