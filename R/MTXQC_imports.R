#### MTXQC_fcn - Import data


path <- "test/"
files <- list.files(path = path, pattern = "*.csv")

for (file in files) {
  perpos <- which(strsplit(file, "")[[1]] == ".")
  assign(
    gsub(" ","",substr(file, 1, perpos - 1)), 
    read.csv(paste(path,file,sep = "")))
}
