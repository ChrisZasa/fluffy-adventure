library(reshape2)

#clean environment
rm(list = ls())

######## Combine files of two batches

subbatch = "brain/"
rawextr = "raw_exp/"

subfolder1 = "sf1/"
subfolder2 = "sf2/"

#location of additional calibration curves (setup: heart and blood)
addq = "raw_addQ/"

### Alkanes
alk1 = read.csv(paste0(subbatch, rawextr, subfolder1, 'Alcane_Intensities.csv'), T)
alk2 = read.csv(paste0(subbatch, rawextr, subfolder2, 'Alcane_Intensities.csv'), T)

alkc = rbind(alk1, alk2)
write.csv(alkc, paste0(subbatch,'input/gc/Alcane_Intensities.csv'), row.names = F)


### CinAcid
cin1 = read.csv(paste0(subbatch, rawextr, subfolder1, 'InternalStandard.csv'), T)
cin2 = read.csv(paste0(subbatch, rawextr, subfolder2, 'InternalStandard.csv'), T)

cinc = merge(cin1, cin2, all = TRUE)

write.csv(cinc, paste0(subbatch,'input/gc/InternalStandard.csv'), row.names = F)


### PeakDensities

pd1 = read.csv(paste0(subbatch, rawextr, subfolder1, 'PeakDensities-Chroma.csv'), T)
pd2 = read.csv(paste0(subbatch, rawextr, subfolder2, 'PeakDensities-Chroma.csv'), T)


pdc = rbind(pd1, pd2)
write.csv(pdc, paste0(subbatch,'input/gc/PeakDensities-Chroma.csv'), row.names = F)


### mz73


mz1 = read.csv(paste0(subbatch, rawextr, subfolder1, 'MassSum-73.csv'), T)
mz2 = read.csv(paste0(subbatch, rawextr, subfolder2, 'MassSum-73.csv'), T)

mzc = rbind(mz1, mz2)
write.csv(mzc, paste0(subbatch,'input/gc/MassSum-73.csv'), row.names = F)


#### quantities

q1 = read.csv(paste0(subbatch, rawextr, subfolder1, 'quantMassAreasMatrix.csv'), T)
q2 = read.csv(paste0(subbatch, rawextr, subfolder2, 'quantMassAreasMatrix.csv'), T)

addq = read.csv("raw_addQ/addq_rawdata.csv", T)

qc = merge(q1, q2)
qc2 = merge(qc, addq, all = TRUE)
write.csv(qc2, paste0(subbatch,'input/quant/quantMassAreasMatrix.csv'), row.names = F)



#### incorporation

inc1 = read.csv(paste0(subbatch, rawextr, subfolder1, 'DataMatrix.csv'), T)
inc2 = read.csv(paste0(subbatch, rawextr, subfolder2, 'DataMatrix.csv'), T)

incc = merge(inc1, inc2)
write.csv(incc, paste0(subbatch,'input/inc/DataMatrix.csv'), row.names = F)



#### MIDS

mid1 = read.csv(paste0(subbatch, rawextr, subfolder1, 'pSIRM_SpectraData.csv'), T)
mid2 = read.csv(paste0(subbatch, rawextr, subfolder2, 'pSIRM_SpectraData.csv'), T)

midc = rbind(mid1, mid2)
write.csv(midc, paste0(subbatch,'input/inc/pSIRM_SpectraData.csv'), row.names = F)


