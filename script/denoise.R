# R>=4, need distances, Rcpp

# Rscript /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/script/denoise.R \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/Counts.txt \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/Spot_metadata.csv \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/Spot_level_haralick_features.csv \
# 5 0.2 1000 5 5 10 1 0.001 \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/

######  read arguments  ################

args=commandArgs(trailingOnly=TRUE)

if (Sys.getenv("RSTUDIO") == "1") # for dev purpose
{
  setwd("/project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data")
  args=c("Counts.txt","Spot_metadata.csv",
    "Spot_level_haralick_features.csv",5,0.2,1000,5,5,10,1,0.001,
    "/project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise",
    ".")  
}

counts_fn=args[1] # Expression matrix.
spot_meta_fn=args[2] # Spot metadata, contains X, Y coordinates.
image_features_fn=args[3] # Extracted image features.
N_PC=as.numeric(args[4]) # Number of PCs to use. positive integers. set to -1 to disable, and use the original IF matrix
R_ratio=as.numeric(args[5]) # Spot neighborhood radius ratio, 0-1, radius=R*min(xmax-xmin,ymax-ymin)
U=as.numeric(args[6]) # perplexity, Tunable
K=as.numeric(args[7]) # latent space dimension, Tunable
S=as.numeric(args[8]) # graph edge weight upper bound, tunable
LAMBDA=as.numeric(args[9]) # regularizer, tunable
L_E=as.numeric(args[10]) # regularizer for denoising
margin=as.numeric(args[11]) # Margin for bisection search, smaller = slower => accuracy up 
script_path=args[12] # path to the software folder
output_path=args[13] # output path

######  load environment  #############

library(distances)
library(Rcpp)

if (Sys.getenv("RSTUDIO") == "1")
{
  source(paste(script_path,"/script/bisection/bisection.R",sep=""))
}else
{
  sourceCpp(paste(script_path,"/script/bisection/bisection_par.cpp",sep=""))
}

source(paste(script_path,"/script/denoise_functions.R",sep=""))

# (1) in practice, we specify an even stronger sufficient condition
# regarding the relationship between d(y) and 1-(p_n_n + t(p_n_n))/2
# this part is tricky. need to think about how to justify 
# for this theorectically, in the best way
# (2) test this in simulation carefully but we want to be vague 
# about this in the paper
Power_tsne_factor=7

###########  read input  ##############

E=as.matrix(read.table(counts_fn,row.names = 1,header = T))
C=read.csv(spot_meta_fn,row.names = 1,header = T,stringsAsFactors=F)
R=min(max(C$X)-min(C$X),max(C$Y)-min(C$Y))*R_ratio
IF=as.matrix(read.csv(image_features_fn,row.names=1,header=T))
if (any(rownames(E)!=rownames(C)) || any(rownames(E)!=rownames(IF)))
  {stop("Spot IDs mismatch!")}
if (any(is.na(E))
  {stop("Input counts matrix contains null values!")}

dir.create(output_path,recursive = TRUE)
#######  initialization  ###############

# proximity matrix
P=as.matrix(dist(C[,c("X","Y")],diag = FALSE))<=R 
table(P)
diag(P)=FALSE # set self proximity to 0
P[lower.tri(P)]=F # for ensuring alpha_n1_n2=alpha_n2_n1

# Initialize graph weight matrix, symmetrical.
# enforce a_n1_n2 == a_n2_n1
ALPHA=diag(dim(IF)[1])
ALPHA[]=S/2*P

# Preprocessing F matrix
IF=scale(IF,center=TRUE,scale=TRUE)
if (N_PC>0) {IF=prcomp(IF)$x[,1:min(N_PC,dim(IF)[2])]} # Preprocessing the features with PCA, optional

#######  spatial denoise  ################

# similarity between spots in the original image feature space.
sigma_n=find_sigma_n_all(t(IF),U,margin)
euc_dist2=as.matrix(dist(IF,diag = FALSE))^2
# probably will need some work to improve the computational efficiency
# of p_n_n. future work
p_n_n=sapply(1:dim(IF)[1], 
  function(n) calculate_spot_dist(n,sigma_n,euc_dist2))
p_nn_tsne=1-(p_n_n + t(p_n_n))/2
p_nn_tsne=p_nn_tsne^(dim(IF)[1]/Power_tsne_factor)
# hist(p_nn_tsne)

# Build a graph based on image features and spot closeness
while (1==1)
{
  ALPHA[]=optim(as.vector(ALPHA),fn=fr,gr=gr,method='L-BFGS-B',
    lower=0,upper=P*S,data=list(LAMBDA, p_nn_tsne, P , K),
    control=list(trace=T))$par 
  
  if (sum(ALPHA>0)>0) {break}
  cat("Lambda is set too large (trivial solution)!\n")
  LAMBDA=LAMBDA*0.8
  cat(paste("Resetting to",LAMBDA,"\n"))
}

table(ALPHA>0)
table(P==T)
if (sum(ALPHA>0)/sum(P==T)>0.9) {warning("Lambda may be set too small!")}
if (sum(ALPHA>0)/sum(P==T)<0.1) {warning("Lambda may be set too large!")}
ALPHA=ALPHA+t(ALPHA)

# Denoise expression matrix based on the graph
G = ALPHA/S
W = solve(diag(1, dim(G)[1]) + L_E * (diag(rowSums(G)) - G))
E_denoised =  W %*% E

# retrieve the embedding space
Q = diag(1, dim(ALPHA)[1])+4*diag(rowSums(ALPHA))-4*ALPHA
Q_inv=solve(Q)
Q_inv=(Q_inv+t(Q_inv))/2

tmp=Q
tmp[]=1/dim(Q_inv)[1]
tmp=diag(dim(Q_inv)[1])-tmp

decomposition=eigen(tmp %*% Q_inv %*% tmp)
Y=decomposition$vectors %*% diag(decomposition$values)^0.5
Y=Y[,1:K]

#########  output  ################

rownames(E_denoised)=rownames(E)
write.table(round(E_denoised,d=5),
  paste(output_path,'/Denoised_matrix.txt',sep=""),
  sep='\t',row.names = T,col.names = T,quote=F)

rownames(G)=colnames(G)=rownames(E)
write.table(round(G,d=5),
  paste(output_path,'/Detected_graph.txt',sep=""),
  sep='\t',row.names = T,col.names = T,quote=F)

rownames(Y)=rownames(E)
colnames(Y)=paste("Component",1:dim(Y)[2],sep="")
write.table(round(Y,d=5),
  paste(output_path,'/Latent_space.txt',sep=""),
  sep='\t',row.names = T,col.names = T,quote=F)
