#### MTXQC_fcn - Import data
# 
# 
# path <- "test/"
# files <- list.files(path = path, pattern = "*.csv")
# 
# for (file in files) {
#   perpos <- which(strsplit(file, "")[[1]] == ".")
#   assign(
#     gsub(" ","",substr(file, 1, perpos - 1)), 
#     read.csv(paste(path,file,sep = "")))
# }
# 
# 
# ### import params
# files <- list.files(path = path_setup, pattern = "*.csv")
# 
# for (file in files) {
#   extr <- which(strsplit(file, "")[[1]] == ".")
#   assign(
#     gsub(" ","",substr(file, 1, extr - 1)), 
#     read.csv(paste(path_setup, file, sep = "")))
#   message("Imported: ", file)
# }
# 
# 
# 
# ### List version
# library(foreign)
# ## create and view an object with file names and full paths
# (f <- file.path("https://stats.idre.ucla.edu/stat/data", c("auto.dta","cancer.dta", "efa_cfa.dta")))
# 
# d <- lapply(f, read.dta)
# str(d, give.attr = FALSE)
# lapply(d, names)
# 
# names(d) <- gsub(".*/(.*)\\..*", "\\1", f)
# 
# summary(d$cancer)
# 
# 
#   ##version for MTXQC
# 
#   path <- "test/"
#   files <- list.files(path = path, pattern = "*.csv")
#   
#   #import
#   f <- file.path("test/", files)
#   d <- lapply(f, read.csv)
#   #add names to each dataframe in the list
#   names(d) <- gsub(".*/(.*)\\..*", "\\1", f)
#   
#   
#   ## create a function
#   import_content <- function(dir, mode){
#     
#     files <- list.files(path = dir, pattern = "*.csv")
#     
#     ## check number of files imported (params = 2, input = 2)
#     n_files = length(files)
#     
#     if (mode == "params") {
#       if (n_files >= 1) {
#         message("Parameter files imported!", dir)
#       } else {
#         message("Missing parameter file detected - check project folder")
#       }
#     }
#     
#     if (mode == "input") {
#       if (n_files == 2) {
#         message("All files imported - annotation and sample extracts file!")
#       } else {
#         message("Missing parameter file detected - annotation and sample extracts file")
#       }
#     }  
# 
#     #import
#     f <- file.path(dir, files)
#     d <- lapply(f, read.csv)
#     
#     #add names
#     names(d) <- gsub(".*/(.*)\\..*", "\\1", f)
#   
#     return(d)
#   }
#   
#   
#   
