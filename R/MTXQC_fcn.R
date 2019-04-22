####### MTXQC functions - general #################

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


file_shaping = function(dataframe, shape = "long", file_annotation, type, inc = FALSE, complete_ann = FALSE) {
### Function: File-shaping
#' including:
#' melting of wide to long format
#' checkup File instead of file
#' adding annotation - File,Type
#' add Batch_Id based on File name
#' WHATS WITH METABOLITE_MANUAL
  
  #check File
  colnames(dataframe)[grepl("file", colnames(dataframe))] <- "File" 
  
  #metmax shaping: test if row.load and ri are still present in the input
  test_row = "row.load" %in% names(dataframe)
  test_metabolite_metmax = "Metabolite" %in% names(dataframe)
  
  if (shape == "wide") {
  
    #row.load and ri still present
    if (test_row == TRUE) {
      
      idx_ri = grep("ri", colnames(dataframe))
      idx_row = grep("row.load", colnames(dataframe))
      df = dataframe[,c(-idx_row, -idx_ri)]
     
      dataframe = reshape2::melt(df, id.vars = c("name", "mass"), 
                                  value.name = "PeakArea",
                                  variable.name = "File", 
                                  na.rm = TRUE)
  
      if (inc == TRUE) {
        colnames(dataframe)[grepl("name", colnames(dataframe))] <- "Metabolite_manual"
        colnames(dataframe)[grepl("PeakArea", colnames(dataframe))] <- "MID_Intensity"
        colnames(dataframe)[grepl("mass", colnames(dataframe))] <- "Mass_mz"
      } else {
        colnames(dataframe)[grepl("name", colnames(dataframe))] <- "Metabolite_manual"
        colnames(dataframe)[grepl("mass", colnames(dataframe))] <- "QuantMasses"
      }
    }
    
    #row.load and ri already deleted, but name and mass still present
    if (test_row == FALSE) {
      
      colnames(dataframe)[grepl("name", colnames(dataframe))] <- "Metabolite_manual"
      colnames(dataframe)[grepl("mass", colnames(dataframe))] <- "Mass_mz"
      
      if (inc == FALSE) {
        
        data_proc = dataframe
        dataframe = reshape2::melt(data_proc, id.vars = c("Metabolite", "QuantMasses"), 
                                   value.name = "PeakArea", variable.name = "File", na.rm = TRUE)
      } else {
        
        if (test_metabolite_metmax == TRUE) {
          data_proc = dataframe
          dataframe = reshape2::melt(data_proc, id.vars = c("Metabolite"), 
                                   value.name = "MID_Intensity", variable.name = "File", na.rm = TRUE)
          
          colnames(dataframe)[grepl('MID_Intensity', colnames(dataframe))] <- "Inc"
          
        } else {
          data_proc = dataframe
          dataframe = reshape2::melt(data_proc, id.vars = c("Metabolite_manual", "Mass_mz"), 
                                     value.name = "MID_Intensity", variable.name = "File", na.rm = TRUE)
        }
        
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
  
  #include only data from files == sample
  if (type == "sample") {
    dataframe = subset(dataframe, dataframe$Type == as.character(type))
  }
  
  return(dataframe)

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
