###### Combine dataset derived from Metmax

subbatch = "exp173_metmax/"
rawextr = "raw_data"

### Peak Areas

area1 = read.csv(paste0(subbatch, rawextr, "/area_4.csv"), TRUE, sep = ';')
area2 = read.csv(paste0(subbatch, rawextr, "/area_9.csv"), TRUE, sep = ';')
area3 = read.csv(paste0(subbatch, rawextr, "/area_21.csv"), TRUE, sep = ';')

area_all = merge(area1, area2, all = TRUE)
area_all = merge(area_all, area3, all = TRUE)

write.csv(area_all, paste0(subbatch, 'input/metmax/areamatrix.csv'), row.names = FALSE)


### mids

mid1 = read.csv(paste0(subbatch, rawextr, "/mids_4.csv"), TRUE, sep = ';')
mid2 = read.csv(paste0(subbatch, rawextr, "/mids_9.csv"), TRUE, sep = ';')
mid3 = read.csv(paste0(subbatch, rawextr, "/mids_21.csv"), TRUE, sep = ';')

mid_all = merge(mid1, mid2, all = TRUE)
mid_all = merge(mid_all, mid3, all = TRUE)

write.csv(mid_all, paste0(subbatch, 'input/metmax/mid_all.csv'), row.names = FALSE)


### mz
mz1 = read.csv(paste0(subbatch, rawextr, "/mz73_4.csv"), TRUE, sep = ';')
mz2 = read.csv(paste0(subbatch, rawextr, "/mz73_9.csv"), TRUE, sep = ';')
mz3 = read.csv(paste0(subbatch, rawextr, "/mz73_21.csv"), TRUE, sep = ';')

mz_all = merge(mz1, mz2, all = TRUE)
mz_all = merge(mz_all, mz3, all = TRUE)

write.csv(mz_all, paste0(subbatch, 'input/metmax/mz73_all.csv'), row.names = FALSE)
