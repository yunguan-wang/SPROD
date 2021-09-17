if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rhdf5",update=F)
library("rhdf5")
h5closeAll()

setwd("/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_190921_19/")
# if failure, see if other sessions are reading this h5 file at the same time
h5f = H5Fopen("denoised_cts.h5df") 
h5f

denoised=h5f$denoised
names(denoised)

h5closeAll()
