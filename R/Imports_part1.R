#### MTXQC_part1 IMPORTS to source()

#### import annotation ----- path_setup/input ---------------

ann_idx <-  as.character(setup_params[which(setup_params$Parameter == "ann"), "Value"])

if (file.exists(file.path(paste0(path_setup, set_input, ann_idx)))) { 
  ann  <-   read.csv(paste0(path_setup, set_input, ann_idx), T)
  
  if (ncol(ann) == 1) {
    ann = read.csv(paste0(path_setup, set_input, ann_idx), T, sep = ";")
    message("colon separator detected: ", ann_idx)
  }
} else {
  message("FATAL ERROR: Annotation file missing: ", ann_idx)
  knitr::knit_exit()
}

check_emptyfile(ann, TRUE, ann_idx)

#Cell count / extract
se_idx <-  as.character(setup_params[which(setup_params$Parameter == "sample_ext"), "Value"])

if (file.exists(file.path(paste0(path_setup, set_input, se_idx)))) { 
  
  data_extracts <-  read.csv(paste0(path_setup, set_input, se_idx), T)
  
  if (ncol(data_extracts) == 1) {
    data_extracts = read.csv(paste0(path_setup, set_input, se_idx), T, sep = ";")
    message("colon separator detected: ", se_idx)
  }
} else {
  message("FATAL ERROR: Sample_extracts file missing: ", se_idx)
  knitr::knit_exit()
}

check_emptyfile(data_extracts, TRUE, se_idx)


#check if addq-values exits if required
#Check if this whole addQ-file thing is needed!!!
nick_addq_file <- as.character(file_spec[which(file_spec$AssociatedFile == "addQ_file"), "Filename"])
idx_addq <- as.character(setup_params[which(setup_params$Parameter == "addQ"), "Value"])

if (idx_addq == "yes") {
  if (!file.exists(file.path(paste0(path_setup, set_input, "add_quant/", nick_addq_file)))) {
    message("Please specify the correct file containing additional Quant1-values in the ExperimentalSetup!")
    knitr::knit_exit()
  } else {
    message("Required table containing additional Quant1-values detected!", nick_addq_file)
  }
} else {
  message("Experimental setup does not include additional quantification standards!")
}

#### import: input/gc



### (1) GC-MS performance files ####
## import (i) cinnamic acid / internal standard peak areas, 
## (ii) alkane intensities,
## (iii) mz 73 intensities and 
## (iv) total peak densities from Chromatof

  # Cinnamic acid / Internal Standard ---------------------------------------------
  ca_idx <-  setup_params[which(setup_params$Parameter == "instd"), "Value"]
  mm_ca_idx <- id_settings[which(id_settings$Parameter == "intstd"), "Value"]
  
  if (ca_idx == TRUE) {
    if (mm_ca_idx == TRUE) {
      
      nick_ca <-  as.character(file_spec[which(file_spec$AssociatedFile == "cin_acid"), "Filename"])
      ca_path <- paste0(path_setup, set_input,'gc/', nick_ca)
      
      if (file.exists(ca_path) == TRUE) {
        cinacid  <-   read.csv(ca_path, T, sep = ',')
        
        if (ncol(cinacid) == 1) {
          cinacid  <-   read.csv(ca_path, T, sep = ';')
        }
        
        #Which internalstandard exported in InternalStandard.csv
        intstd = unique(cinacid$Metabolite)
        message("Defined Internal Standard: ", intstd)
        
        check_emptyfile(cinacid, FALSE, "InternalStandard.csv")
        is_exit = 0
        
      } else {
        message("WARNING: No file detected: InternalStandard.csv")
        is_exit = 1
      }
    } else {
      
      imp_error = imp_error + 1
      is_exit = 1
      message("WARNING: No internal standard for data input format defined!")
      
    }
  } else {
    is_exit = 1
    message("No internal standard in Experimental Setup defined!")
  }

  # Alkane intensities ---------------------------------------------------------
  alkane_idx <- as.character(id_settings[which(id_settings$Parameter == "alkanes"), "Value"])
  nick_alk <-  as.character(file_spec[which(file_spec$AssociatedFile == "alkane_int"), "Filename"])
  
  if (alkane_idx == TRUE) {
    data_alk  <-   read.csv(paste0(path_setup, set_input,'gc/', nick_alk), T, sep = ',')
    colnames(data_alk)[grepl("file", colnames(data_alk))] = "File"
    check_emptyfile(data_alk, FALSE, nick_alk)
    
  } else {
    imp_error = imp_error + 1
    message("WARNING: No alkane file defined for this input format!")
  }

  # m/z 73 Table ---------------------------------------------------------------
    nick_73 <-  as.character(file_spec[which(file_spec$AssociatedFile == "mz_73"), "Filename"])   
    mm_73_idx <-  id_settings[which(id_settings$Parameter == "mz"), "Value"]
    
    if (mm_73_idx != "") {
      
      if (!file.exists(file.path(paste0(path_setup, set_input,'gc/', nick_73)))) {
        
        message("WARNING: File missing: ", nick_73)
        message("Please review MTXQC_ExperimentalSetup!")
        soa_73_exit = 1
        
      } else {
        
        data_73  <-   read.csv(paste0(path_setup, set_input,'gc/', nick_73), T, sep  =  ",")
        check_emptyfile(data_73, FALSE, nick_73)
        soa_73_exit = 0
      }
    } else {
      
      imp_error = imp_error + 1
      soa_73_exit = 1
      message("WARNING: No file with m/z 73 values defined for this input format!")
      
    }

  #Peak-Chroma Table -----------------------------------------------------------
    nick_peaks <-  as.character(file_spec[which(file_spec$AssociatedFile == "peak_densities"), "Filename"])
    mm_peaks_idx <- id_settings[which(id_settings$Parameter == "peakchroma"), "Value"]
    
    if (mm_peaks_idx == TRUE) {
      if (!file.exists(file.path(paste0(path_setup, set_input,'gc/', nick_peaks)))) { 
        message("WARNING: File missing: ", nick_peaks)
        message("Please review MTXQC_ExperimentalSetup!")
        soa_exit = 1
        
      } else {
        total_peaks  <-   read.csv(paste0(path_setup, set_input,'gc/', nick_peaks), T)
        check_emptyfile(total_peaks, FALSE, nick_peaks)
        soa_exit = 0
      }
      
    } else {
      imp_error = imp_error + 1
      soa_exit = 1
      message("WARNING: No file defined for this input format! No sum of area normalisation! ", nick_peaks)
    }

#### (2) Absolute quantification ####
## (2.1) Exported peak areas according to top5 or pTop5 approach

if ((params$updated == "PeakArea") | (params$updated == "both")) {
  
  #incorporation of manualy validated data
  nick_samples <-  as.character(file_spec_manVal[which(file_spec_manVal$AssociatedFile == "sample_area"), "Filename"]) 
  if (!file.exists(file.path(paste0(path_setup,  set_input, 'quant/', nick_samples)))) {
    
    message("FATAL ERROR: Manual validated peak area file not detected: ", nick_samples)
    knitr::knit_exit()
    
  } else {
    data_area  <-   read.csv(paste0(path_setup,  set_input, 'quant/', nick_samples), T)
    
    if (ncol(data_area) == 1) {
      data_area  <-   read.csv(paste0(path_setup,  set_input, 'quant/', nick_samples), T, sep = ";")
    }
    check_emptyfile(data_area, TRUE, nick_samples)   
  }
  
} else {
  #only original data
  nick_samples <-  as.character(file_spec[which(file_spec$AssociatedFile == "sample_area"), "Filename"]) 
  
  if (!file.exists(file.path(paste0(path_setup,  set_input, 'quant/', nick_samples)))) {
    message("FATAL ERROR: Peak area table not detected: ", nick_samples)
    knitr::knit_exit()
    
  } else {
    data_area  <-   read.csv(paste0(path_setup,  set_input, 'quant/', nick_samples), T)
    
    if (ncol(data_area) == 1) {
      data_area  <-   read.csv(paste0(path_setup,  set_input, 'quant/', nick_samples), T, sep = ";")
    }
    check_emptyfile(data_area, TRUE, nick_samples)    
  }
  
}


#### (3) Incorporation data ####
## only needed for pSIRM experiment
## (i) MID exports and (ii) calculated isotope incorporation
inc_idx <-  setup_params[which(setup_params$Parameter == "data"), "Value"]
evalinc_idx <-  setup_params[which(setup_params$Parameter == "substr"), "Value"]

mm_inc_idx <- id_settings[which(id_settings$Parameter == "inc"), "Value"]

if ((params$updated == "none") | (params$updated == "PeakArea")) { #run if not Inc manually validated
  if (inc_idx != "qMTX") {
    if (mm_inc_idx != FALSE) {
      if (evalinc_idx != FALSE | evalinc_idx == "none") {
        #Mass isotopomer distributions
        nick_psirm <-  as.character(file_spec[which(file_spec$AssociatedFile == "pSIRM_se"), "Filename"])
        
        if (!file.exists(file.path(paste0(path_setup, set_input,'inc/', nick_psirm)))) {
          message("FATAL ERROR: Essential file not detected: ", nick_psirm)
          knitr::knit_exit()
          
        } else {
          data_mid  <-   read.csv(paste0(path_setup, set_input,'inc/', nick_psirm), T)
          
          if (ncol(data_mid) == 1) {
            data_mid  <-   read.csv(paste0(path_setup,  set_input, 'inc/', nick_psirm), T, sep = ";")
          }
          check_emptyfile(data_mid, TRUE, nick_psirm)
        }
        
        
        #13C-Incorporation
        nick_inc <-  as.character(file_spec[which(file_spec$AssociatedFile == "inc"), "Filename"]) 
        
        if (!file.exists(file.path(paste0(path_setup, set_input,'inc/', nick_inc)))) {
          message("FATAL ERROR: Essential file missing! ", nick_inc)
          message("Are you sure it's a pSIRM experiment???")
          knitr::knit_exit()
        } else {
          data_inc  <-   read.csv(paste0(path_setup, set_input,'inc/', nick_inc), T)
          
          if (ncol(data_inc) == 1) {
            data_inc  <-   read.csv(paste0(path_setup,  set_input, 'inc/', nick_inc), T, sep = ";")
          }
          check_emptyfile(data_inc, TRUE, nick_inc)
        }
        
      } else {
        message("No stable isotope incorporation evaluated.")
      }
    } else (
      message("No MID and 13Inc-calculation defined for this input format.")
    )
  } else {
    message("It's not a pSIRM experiment!")
  }
} else {
  #run if Inc == manually validated
  if (inc_idx != "qMTX") {
    if (mm_inc_idx != FALSE) {
      if (evalinc_idx != FALSE | evalinc_idx == "none") {
        
        #Mass isotopomer distributions
        nick_psirm <-  as.character(file_spec_manVal[which(file_spec_manVal$AssociatedFile == "pSIRM_se"), "Filename"])
        
        if (!file.exists(file.path(paste0(path_setup, set_input,'inc/', nick_psirm)))) {
          message("FATAL ERROR: Essential file missing: ", nick_psirm)
          message("Did you really validated the incorporation data?")
          knitr::knit_exit()
        } else {
          data_mid  <-   read.csv(paste0(path_setup, set_input,'inc/', nick_psirm), T)
          
          if (ncol(data_mid) == 1) {
            data_mid  <-   read.csv(paste0(path_setup,  set_input, 'inc/', nick_psirm), T, sep = ";")
          }
          check_emptyfile(data_mid, TRUE, nick_psirm)
        }
        
        
        #13C-Incorporation
        nick_inc <-  as.character(file_spec_manVal[which(file_spec_manVal$AssociatedFile == "inc"), "Filename"]) 
        
        if (!file.exists(file.path(paste0(path_setup, set_input,'inc/', nick_inc)))) {
          message("FATAL ERROR: Essential file missing: ", nick_inc)
          message("Did you really validated the incorporation data?")
          knitr::knit_exit()
        } else {
          data_inc  <-   read.csv(paste0(path_setup, set_input,'inc/', nick_inc), T)
          
          if (ncol(data_inc) == 1) {
            data_inc  <-   read.csv(paste0(path_setup,  set_input, 'inc/', nick_inc), T, sep = ";")
          }
          check_emptyfile(data_inc, TRUE, nick_inc)
        }
        
      } else {
        message("No stable isotope incorporation evaluated.")
      }
    } else {
      message("No MID and 13Inc-calculation defined for this input format.")
    }
  } else {
    message("It's not a pSIRM experiment!")
  }
}


#### Check for import errors #####
imp_error <- MTXQCp1_checkinput(data_extracts, ann)

if (imp_error != 0) {
  message("Check the imported files messages! Number of files without import: ", imp_error)
} else {
  message("Annotation and Sample_extract.csv correctly imported!")
}