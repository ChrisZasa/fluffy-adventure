######## Combine files of two batches

### Alkanes

alk1 =read.csv('clone2646/input_MTXQC_2646/gc/Alcane_Intensities.csv', T)
alk2 = read.csv('clone2646/input_MTXQC_2646SL/gc/Alcane_Intensities.csv', T)

alkc = rbind(alk1, alk2)
write.csv(alkc, 'clone2646/_forMTXQCvX/Alcane_Intensities.csv', row.names = F)


### CinAcid

cin1 = read.csv('clone2646/input_MTXQC_2646/gc/cinAcid.csv', T)
cin2 = read.csv('clone2646/input_MTXQC_2646SL/gc/cinAcid.csv', T)

cinc = merge(cin1, cin2, all = TRUE)
write.csv(cinc, 'clone2646/_forMTXQCvX/cinacid.csv', row.names = F)


### PeakDensities

pd1 = read.csv('clone2646/input_MTXQC_2646/gc/PeakDensities-Chroma.csv', T)
pd2 = read.csv('clone2646/input_MTXQC_2646SL/gc/PeakDensities-Chroma.csv', T)

pdc = rbind(pd1, pd2)
write.csv(pdc, 'clone2646/_forMTXQCvX/PeakDensities-Chroma.csv', F)


### mz73

mz73 = read.csv('clone2646/input_MTXQC_2646/gc/MassSum-73.csv', T)
mz73_2 = read.csv('clone2646/input_MTXQC_2646SL/gc/MassSum-73.csv', T)

mz_c = rbind(mz73, mz73_2)
write.csv(mz_c, 'clone2646/_forMTXQCvX/MassSum-73.csv', F)

#### quantities

q1 = read.csv('clone2646/input_MTXQC_2646/quant/quantMassAreasMatrix_2646split.csv', T)
q2 = read.csv('clone2646/input_MTXQC_2646SL/quant/quantMassAreasMatrix_2646sl.csv', T)

q_c = merge(q1, q2, all = T)
write.csv(q_c, 'clone2646/_forMTXQCvX/quantMassAreasMatrix.csv', F)

#### incorporation

i1 =read.csv('clone2646/input_MTXQC_2646/inc/DataMatrix.csv', T)
i2 =read.csv('clone2646/input_MTXQC_2646SL/inc/DataMatrix.csv', T)

i_c = merge(i1, i2)
write.csv(i_c, 'clone2646/_forMTXQCvX/DataMatrix.csv', F)

#### MIDS

s1 = read.csv('clone2646/input_MTXQC_2646/inc/pSIRM_SpectraData.csv', T)
s2 = read.csv('clone2646/input_MTXQC_2646SL/inc/pSIRM_SpectraData.csv', T)

s_c = rbind(s1, s2)
write.csv(s_c, 'clone2646/_forMTXQCvX/pSIRM_SpectraData.csv', F)


