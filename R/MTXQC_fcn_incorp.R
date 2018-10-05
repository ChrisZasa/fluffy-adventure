####### MTXQC functions - Stable Isotope Incorporation ######


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
