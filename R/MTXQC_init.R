initiate_MTXQC = function(subfolder = "") {
#' This script creates the required folder structure for 
#' the setup of an MTXQC project, e.g., defining input and
#' output folders including subfolders.
#' In addition it checks installed packages that are required
#' for all MTXQC modules.
#' 
#' input: Subfolder name if required, default

  
#set directory of the experimental setup
  if (subfolder == "") {
    path_setup = ""
  } else {
    path_setup = paste0(subfolder, "/")
  }
  set_input  =  'input/'
  set_output  =  'output/'

  
#### Create input folder ####  
  ifelse(!dir.exists(file.path(paste0(path_setup))), 
    dir.create(paste0(path_setup)), 
    FALSE)

#### Create input folder ####  
  ifelse(!dir.exists(file.path(paste0(path_setup, set_input))), 
         dir.create(paste0(path_setup, set_input)), 
         FALSE)
  
  
  ifelse(!dir.exists(file.path(paste0(path_setup, set_input, 'quant'))), 
         dir.create(paste0(path_setup, set_input, 'quant')), 
         FALSE)
  
  ifelse(!dir.exists(file.path(paste0(path_setup, set_input, 'gc'))), 
         dir.create(paste0(path_setup, set_input, 'gc')), 
         FALSE)
  
  ifelse(!dir.exists(file.path(paste0(path_setup, set_input, 'inc'))), 
         dir.create(paste0(path_setup, set_input, 'inc')), 
         FALSE)

  ifelse(!dir.exists(file.path(paste0(path_setup, set_input ,'add_quant'))),
         dir.create(paste0(path_setup, set_input, 'add_quant')), 
         FALSE)
  
  ifelse(!dir.exists(file.path(paste0(path_setup, 'metmax'))),
    dir.create(paste0(path_setup, set_input, 'metmax')), 
    FALSE)
  
#### Create output folder ####

  ifelse(!dir.exists(file.path(paste0(path_setup, set_output))), 
         dir.create(paste0(path_setup, set_output)), 
         FALSE)
  
  
  ifelse(!dir.exists(file.path(paste0(path_setup, set_output, 'quant'))), 
         dir.create(paste0(path_setup, set_output, 'quant')), 
         FALSE)
  
  ifelse(!dir.exists(file.path(paste0(path_setup, set_output, 'gc'))), 
         dir.create(paste0(path_setup, set_output, 'gc')), 
         FALSE)
  
  ifelse(!dir.exists(file.path(paste0(path_setup, set_output, 'inc'))), 
         dir.create(paste0(path_setup, set_output, 'inc')), 
         FALSE)

#### Create figure folder ####
  
  ifelse(!dir.exists(file.path(paste0(path_setup, "figure"))), 
         dir.create(paste0(path_setup, "figure")), 
         FALSE)


#### Libraries check ####
  list_packages <- c("ggplot2", "reshape2", "RColorBrewer", "plyr", "scales",
                     "kableExtra", "magrittr", "grid", "gridExtra", "shiny", "knitr",
                     "rmarkdown", "formatR", "readr", "tinytex", "gplots", "devtools",
                      "dplyr", "ggrepel", "tufte")
  
  new_packages <- list_packages[!(list_packages %in% 
                                     installed.packages()[,"Package"])]
  
  if (length(new_packages)) {
    install.packages(new_packages, repos = 'http://cran.us.r-project.org')
  }

}
