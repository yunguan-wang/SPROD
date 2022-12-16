# Rversion>=4, need distances, dplyr, optparse,
# if diag=T, need ggplot2, gridExtra
######  load environment  ###########
suppressPackageStartupMessages(library(distances))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
#######  read arguments  ##############
option_list = list(
  make_option(c("-e", "--Exp"), action="store", default=NA,
              type='character',help="Expression matrix"),
  make_option(c("-c", "--Cspot"), action="store",
              default=NA, type='character',
              help="Spot metadata, contains X, Y coordinates."),
  make_option(c("-f", "--ImageFeature"), action="store",
              default=NA,type = 'character',
              help="Extracted image features"),
  make_option(c("-n", "--numberPC"), action="store",
              default = -1,type = 'double',
              help="# Number of PCs to use. positive integers. by default set to -1 to use the number of PCs the same as the numebr of features"),
  make_option(c("-x","--umap"),action="store_true",dest = "umap",
              default = FALSE,type = "logical",
              help ="# Use UMAP of Image features or NOT. By default set to FALSE, which the PCA is used to transform IF."),
  make_option(c("-g","--diagnose"),action="store_true",dest = "diagnose",
              default = FALSE,type = "logical",
              help ="# Generate diagnostic figures or NOT. By default set to FALSE."),
  make_option(c("-r", "--Rratio"), action="store",
              default=0.08,type = 'double',
              help="Spot neighborhood radius ratio, 0-1, radius=R*min(xmax-xmin,ymax-ymin)"),
  make_option(c("-u", "--U"), action="store",
              default = 250,type = 'double',
              help="# perplexity, Tunable"),
  make_option(c("-k", "--K"), action="store",
              default=10,type = 'double',
              help="latent space dimension, Tunable"),
  make_option(c("-l", "--lambda"), action="store",
              default = 1,type = 'double',
              help="# regularizer, tunable"),
  make_option(c("-t", "--L_E"), action="store",
              default=0.4,type = 'double',
              help="regularizer for denoising"),
  make_option(c("-m", "--margin"), action="store",
              default = 0.001,type = 'double',
              help="Margin for bisection search, smaller = slower => accuracy u"),
  make_option(c("-d","--diagNeighbor"),action="store_true",
              default = FALSE, type = "logical",
              help = "To prevent over smooth."),
  make_option(c("-s", "--scriptPath"), action="store",
              default=NA, type ='character',
              help="path to the software folder"),
  make_option(c("-o", "--outputPath"), action="store",
              default=NA, type ='character',
              help="Output path"),
  make_option(c("-p", "--projectID"), action="store",
              default=NA, type = 'character',
              help="project name, First part of names in the output"),
  make_option(c("-q", "--TopN"), action="store",
              default=1000, type = 'double',
              help="Number of edges selected, only used when diagnose mode is on.")
  )
opt = parse_args(OptionParser(option_list=option_list))

###########  define input  ############
counts_fn <- opt$Exp
spot_meta_fn <- opt$Cspot
image_features_fn <- opt$ImageFeature
N_PC=opt$numberPC
diagnose = opt$diagnose
um =opt$umap
R_ratio=opt$Rratio # Spot neighborhood radius ratio, 0-1, radius=R*min(xmax-xmin,ymax-ymin)
U=opt$U # perplexity, Tunable
K=opt$K # latent space dimension, Tunable
LAMBDA=opt$lambda # regularizer, tunable
L_E=opt$L_E # regularizer for denoising
margin=opt$margin # Margin for bisection search, smaller = slower => accuracy up
d = opt$diagNeighbor
script_path=opt$scriptPath # path to the software folder
output_path=opt$outputPath # output path
diagnose=opt$diagnose
project_name=opt$projectID
TopN=opt$TopN

log = file(
  file.path(
    output_path,
    paste(project_name,'.log', sep=''))
  ,open='wt')
sink(log, type='out')
sink(log, type='message')
cat(paste("N_PC:",N_PC,"umap",um,"R_ratio:",R_ratio,"U:",U,"K",K,"LAMBDA:",LAMBDA,"L_E:",L_E,"margin:",margin,"W_diag:",d,"project ID:",project_name,"diagnose mode on:",diagnose,"TopN:",TopN))
cat("\n\n")

source(file.path(script_path,"bisection/bisection.R"))
source(file.path(script_path,"denoise_functions.R"))

# (1) in practice, we specify an even stronger sufficient condition
# regarding the relationship between d(y) and 1-(p_n_n + t(p_n_n))/2
# this part is tricky. need to think about how to justify
# for this theorectically, in the best way
# (2) test this in simulation carefully but we want to be vague
# about this in the paper
Power_tsne_factor=7

##########read input######
cat("Inputing...\n\n")
E=as.matrix(read.table(counts_fn,row.names = 1,header = T))
C=read.csv(spot_meta_fn,row.names = 1,header = T,stringsAsFactors = F)
if ("Z" %in% colnames(C)){
  cat("Z column is detected!\n")
  R=min(max(C$X)-min(C$X),max(C$Y)-min(C$Y),max(C$Z)-min(C$Z))*R_ratio
}else{
  if ("X" %in% colnames(C) & "Y" %in% colnames(C)){
    cat("XY columns are detected but not Z!\n")
    R=min(max(C$X)-min(C$X),max(C$Y)-min(C$Y))*R_ratio
  }else{
    cat ("No correct colnumns are detected! Please check the colnames of spot_metadata!\n")
  }
}
# Make it flexible to accept csv and txt.
if (substr(image_features_fn, nchar(image_features_fn)-3, nchar(image_features_fn)) == '.csv') {
  IF=as.matrix(read.csv(image_features_fn,row.names = 1,header = T))
} else {
  IF=as.matrix(read.csv(image_features_fn,row.names = 1,header = T, sep='\t'))
}

if (any(rownames(E) != rownames(C)) || any(rownames(E)!=rownames(IF))) {
  stop("Spot IDs mismatch!")
  }

#######  initialization  #############
cat("Constructing latent graph based on position and image features...\n\n")
# proximity matrix
if ("Z" %in% colnames(C)){
  P = as.matrix(dist(C[,c("X","Y","Z")],diag = FALSE)) <= R
}else{
  P = as.matrix(dist(C[,c("X","Y")],diag = FALSE)) <= R
}
table(P)
diag(P)=FALSE # set self proximity to 0
P[lower.tri(P)]=F # for ensuring alpha_n1_n2=alpha_n2_n1
# Initialize graph weight matrix, symmetrical.
# enforce a_n1_n2 == a_n2_n1
ALPHA=diag(dim(IF)[1])
#S =1
ALPHA[] = 1/2*P

# Preprocessing the features with PCA or UMAP
IF=scale(IF,center=TRUE,scale=TRUE)
IF0=IF #save the raw image features for plots#
if (um) {
  cat('Image features preprocessing: Umap...\n\n')
  IF=umap::umap(IF)$layout
  cat("IF: UMAP done!\n\n")
} else {
    cat('Image features preprocessing: PCA...\n\n')
    if (N_PC > 0) {
      IF = prcomp(IF)$x[,1:min(N_PC,dim(IF)[2])]
    } else {
      IF=prcomp(IF)$x[,1:dim(IF)[2]]
    }
}

cat("Top 5 rows of input image features:\n\n")
head(IF[1:5,1:2])
cat("\n\n")
#######  spatial denoise  ##############
ptm <- proc.time()
# similarity between spots in the original image feature space.
sigma_n=find_sigma_n_all(t(IF),U,margin)
euc_dist2=as.matrix(dist(IF,diag = FALSE))^2
# probably will need some work to improve the computational efficiency
# of p_n_n. future work
p_n_n=sapply(1:dim(IF)[1],
             function(n) calculate_spot_dist(n,sigma_n,euc_dist2))
if (!all(complete.cases(p_n_n))){
  cat("Error: NA in distance matrix! Please check your image features!\n")
}
p_nn_tsne=1-(p_n_n + t(p_n_n))/2
p_nn_tsne=p_nn_tsne^(dim(IF)[1]/Power_tsne_factor)
# p_nn_tsne could have null values using simulated data.
{
  if(!all(complete.cases(p_nn_tsne))) 
  stop("Error: Latent space proximity matrix contains Nan!\n")

# Build a graph based on image features and spot closeness
while (1==1) {
  ALPHA[]=optim(as.vector(ALPHA),fn=fr,gr=gr,method='L-BFGS-B',
                lower=0,upper=P*1,data=list(LAMBDA, p_nn_tsne, P , K),
                control=list(trace=T))$par

  if (sum(ALPHA>0)>0) {break}
  cat("Lambda is set too large (trivial solution)!\n\n")
  LAMBDA=LAMBDA*0.8
  cat(paste("Resetting to",LAMBDA,"\n\n"))
}

cat("Graph construction succeeded.\n\n")
cat("Number of non-zero edges: \n")
table(ALPHA>0)

if (sum(ALPHA>0)/sum(P==T)>0.9) {warning("Lambda may be set too small!")}
if (sum(ALPHA>0)/sum(P==T)<0.1) {warning("Lambda may be set too large!")}
ALPHA=ALPHA+t(ALPHA)
# Denoise expression matrix based on the graph
G = ALPHA/1
rownames(G)=colnames(G)=rownames(E)
#######WetRun Going on below#######
W = solve(diag(1, dim(G)[1]) + L_E * (diag(rowSums(G)) - G))
if (d) {
  cat("Weight matrix is forced to be diagonal\n\n")
  diag_W=diag(W)
  W[G==0]=0
  diag(W)=diag_W
  W=W/rowSums(W)
}

E_denoised =  W %*% E

# retrieve the embedding space
cat("Retrieving the embedding space...\n\n")
Q = diag(1, dim(ALPHA)[1])+4*diag(rowSums(ALPHA))-4*ALPHA
Q_inv=solve(Q)
Q_inv=(Q_inv+t(Q_inv))/2

tmp=Q
tmp[]=1/dim(Q_inv)[1]
tmp=diag(dim(Q_inv)[1])-tmp

decomposition=eigen(tmp %*% Q_inv %*% tmp)
Y=decomposition$vectors %*% diag(decomposition$values)^0.5
Y=Y[,1:K]
proc.time() - ptm
}
########## The following is only executed with diagnose mode on ##########
if (diagnose) {
  cat("Diagnose mode is on, calculating diagnostic measures...\n\n")
  # only needed for diagnose mode
  #######Output sum of dist#######
  Stack <- data.frame(S4Vectors::stack(G))
  Stack <- Stack[which(Stack$value>0),]
  Stack$x1 <- C[Stack$row,"X"]
  Stack$y1 <- C[Stack$row,"Y"]
  Stack$x2 <- C[Stack$col,"X"]
  Stack$y2 <- C[Stack$col,"Y"]
  if ("Z" %in% colnames(C)){
    Stack$z1 <- C[Stack$row,"Z"]
    Stack$z2 <- C[Stack$col,"Z"]
  }
  ## Generate the t_SNE plot
  set.seed(0)
  tsne <- Rtsne::Rtsne(IF0,check_duplicates = FALSE)
  rownames(tsne$Y)= rownames(IF0)
  colnames(tsne$Y) = c("tsne1","tsne2")
  Stack$tsne11 <- tsne$Y[Stack$row,"tsne1"]
  Stack$tsne12 <- tsne$Y[Stack$row,"tsne2"]
  Stack$tsne21 <- tsne$Y[Stack$col,"tsne1"]
  Stack$tsne22 <- tsne$Y[Stack$col,"tsne2"]
  
  if (um) {
    colnames(IF) = c("umap1","umap2")
    Stack$umap11 <- IF[Stack$row,"umap1"]
    Stack$umap12 <- IF[Stack$row,"umap2"]
    Stack$umap21 <- IF[Stack$col,"umap1"]
    Stack$umap22 <- IF[Stack$col,"umap2"]
    
  } else {
    
    IFump=umap::umap(IF0)$layout
    colnames(IFump) = c("umap1","umap2")
    Stack$umap11 <- IFump[Stack$row,"umap1"]
    Stack$umap12 <- IFump[Stack$row,"umap2"]
    Stack$umap21 <- IFump[Stack$col,"umap1"]
    Stack$umap22 <- IFump[Stack$col,"umap2"]
  }
  
  Stack$rk = rank(Stack$value,ties.method = "random")
  if (nrow(Stack)>=TopN){
    TopN = TopN
  }else{
      TopN = nrow(Stack)
      cat("Number of edges is smaller than TopN, and all edges are included!\n")
    }
  Stack %>%
    arrange(desc(rk)) %>%
    head(TopN) -> Stack_top
  head(Stack_top)
  if ("Z" %in% colnames(C)){
    spa = t(Stack_top[,c("x1","y1","z1","x2","y2","z2")])
  }else{
    spa = t(Stack_top[,c("x1","y1","x2","y2")])
  }
  tsn = t(Stack_top[,c("tsne11","tsne12","tsne21","tsne22")])
  ump = t(Stack_top[,c("umap11","umap12","umap21","umap22")])
  if ("Z" %in% colnames(C)){
    dist_spa=sapply(data.frame(spa), function(x){dist(rbind(x[1:3],x[4:6]))})
  }else{
    dist_spa=sapply(data.frame(spa), function(x){dist(rbind(x[1:2],x[3:4]))})
  }
  dist_tsn=sapply(data.frame(tsn), function(x){dist(rbind(x[1:2],x[3:4]))})
  dist_ump=sapply(data.frame(ump), function(x){dist(rbind(x[1:2],x[3:4]))})
  
  dist_out=sapply(list(dist_spa,dist_tsn,dist_ump), mean,na.rm=T)
  cat(paste("Average distance: \n-spatial:",dist_out[1],
            "\n-image tsne:",dist_out[2],
            "\n-image umap:",dist_out[3],"\n\n"))
  
  # Prepare for diagnostic figures
  library(ggplot2)
  cellNames= data.frame(t(Stack[,1:2]))
  Stack$pair = sapply(cellNames,function(x) {paste0(sort(x)[1],sort(x)[2])})
  
  
  # Plotting diagnostic measures
  # Only plot 5000 edges.
  if (length(table(Stack$pair)) > 5000) {
    Stack_plot = filter(Stack, pair %in% sample(names(table(Stack$pair)),5000))
  } else {
    Stack_plot = Stack
  }
  
  pdf(file.path(output_path, paste(project_name, '_Diagnose_Spatial.pdf', sep='')))
  p = ggplot(Stack_plot,aes(x=x1,y=y1))+
    geom_line(aes(group=pair,color=value,alpha=value))+
    xlab("X")+
    ylab("Y")+
    ggtitle("Spatial")
  print(p)
  dev.off()
  
  pdf(file.path(output_path, paste(project_name, '_Diagnose_TSNE.pdf', sep='')))
  p = ggplot(Stack_plot,aes(x=tsne11,y=tsne12)) +
    geom_line(aes(group=pair
                  ,color=value
                  ,alpha=value))+
    xlab("TSNE1")+
    ylab("TSNE2")+
    ggtitle("TSNE")
  print(p)
  dev.off()
  
  pdf(file.path(output_path, paste(project_name, '_Diagnose_UMAP.pdf', sep='')))
  p = ggplot(Stack_plot,aes(x=umap11,y=umap12)) +
    geom_line(aes(group=pair
                  ,color=value
                  ,alpha=value))+
    xlab("UMAP1")+
    ylab("UMAP2")+
    ggtitle("UMAP")
  print(p)
  dev.off()
}

########## End of diagnose part ########

##########Export correlations #######
#cor1=cor(as.matrix(dist(C)),G)
#cor2=cor(as.matrix(dist(IF)),G)
#Scor1=cor(as.matrix(dist(C)),G,method = "spearman")
#Scor2=cor(as.matrix(dist(IF)),G,method = "spearman")
#Cors=sapply(list(cor1,cor2,Scor1,Scor2),mean,na.rm=T)
#cat(paste("Correlations: \n-physical:",Cors[1],
#          "\n-image:",Cors[2],"\n-physical spearman:",Cors[3],
#          "\n-image spearman:",Cors[4],
#          "\n\n"))

#########  output  ##############
cat("Saving outputs...\n\n")
rownames(E_denoised)=rownames(E)
rownames(Y)=rownames(E)
colnames(Y)=paste0("Component",1:dim(Y)[2])

# setting file names
E_denoised_fn = paste(project_name,'_Denoised_matrix.txt', sep='')
G_fn = paste(project_name,'_Detected_graph.txt', sep='')
Y_fn = paste(project_name,'_Latent_space.txt', sep='')

write.table(round(E_denoised,d=5),
            file.path(output_path,E_denoised_fn),
            sep='\t',row.names = T, col.names = T, quote=F)
write.table(round(G,d=5),
            file.path(output_path,G_fn),
            sep='\t',row.names = T, col.names = T, quote=F)
write.table(round(Y,d=5),
            file.path(output_path,Y_fn),
            sep='\t',row.names = T, col.names = T, quote=F)
sink()


