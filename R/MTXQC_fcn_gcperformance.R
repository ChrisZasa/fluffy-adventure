####### MTXQC functions - GC-Performance ######


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
  if (title == "internalstandard") {
    df_stat <-  ddply(dataframe, c("Batch_Id"), transform, 
                      IntStd_fac = PeakArea / mean(PeakArea))
    
    #Quality check - should be within a factors range of 0.65 and 1.45
    df_stat$IntStd_eval <-  ifelse(df_stat$IntStd_fac >= 0.65 & df_stat$IntStd_fac <= 1.45, 'within', 
                               ifelse(df_stat$IntStd_fac < 0.65, 'below','above'))
    
    #Export CA-normalization information
    write.csv(df_stat, paste0(path_setup,'output/gc/IntStandard_normfactors.csv'), row.names = F)
    
    #calculate cinnamic acid-factor
    batch_stats <-  ddply(df_stat, c("Batch_Id"), transform,
                          n_batch = length(PeakArea),
                          mean_batch = mean(PeakArea, na.rm = TRUE),
                          sd_batch = sd(PeakArea, na.rm = TRUE))
    
    write.csv(batch_stats, paste0(path_setup, 'output/gc/IntStandard_stats.csv'), row.names = F)
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


extract_standards_export <- function(dataframe, met_names = con_se, ann_file = ann) {
  
  #Extract internal standard
  is_subs = subset(met_names, met_names$Standards == "InternalStandard")
  internalstd = unique(is_subs$Lettercode)
  
  #Extract measures peak areas
  df_peak = merge(dataframe, met_names[,c("Metabolite_manual", "Metabolite","Lettercode")])
  df_standard = subset(df_peak, Lettercode %in% internalstd)
  
  #if (!is.null(dim(df_standard))) {
  if (!is.data.frame(df_standard) && !nrow(df_standard)) {
    
    message("Defined internal standard: ", internalstd)
    message("Peak areas detected for internal standard in peak area matrix.")
    
    #Selected columns for wide format
    df_sel = df_standard[,c("Metabolite", "QuantMasses" ,"File", "PeakArea")]
    df_wide = reshape2::dcast(df_sel, Metabolite + QuantMasses ~ File, 
                              value.var = "PeakArea")
    
    #Export
    write.csv(df_wide, paste0(path_setup, set_input, "gc/InternalStandard.csv"), 
              row.names = FALSE)
    
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


