library(scImpute)
# library(SAVER)
#both need raw counts as input
#both functions have other parameters
#run saver in RStudio

input_dir = "~/work_xiao/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/patches/"
# cts_patches = Filter(function(x) grepl("Counts",x), list.files(input_dir))
# for (patch in cts_patches) {
#     cts = read.table(
#         paste(input_dir,patch,sep=''),row.names = 1, header = T
#     )
#     saver_cts=saver(t(as.matrix(cts)), ncores = 16, size.factor = 1)$estimate
#     write.table(saver_cts, paste(input_dir,gsub(".txt", ".saver.txt", patch), sep=''), sep = '\t')
# }


#The input of scimpute is a csv file
cts_patches = Filter(function(x) grepl("raw_cts",x), list.files(input_dir))
for (patch in cts_patches) {
    cts_fn = paste(input_dir,patch,sep='')
    cts = read.table(cts_fn,row.names = 1, header = T)
    if (dim(cts)[2]>dim(cts)[1]) {
        write.table(t(cts), cts_fn, sep = '\t')
        }
    scimpute(cts_fn, infile='txt',outfile='txt', 
             out_dir="~/work_xiao/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/scimpute",
             labeled=FALSE,drop_thre=0.5,Kcluster=5, ncores = 16)
}