#########  required R packages  ############

packages <- c("dirichletprocess")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(
        setdiff(packages, rownames(installed.packages())), 
        repos = "http://cran.us.r-project.org")
}
library("dirichletprocess")

##########  input and parameter  #############

# genes (actually components) on cols
# must be reduced to components by some dimension reduction method.
# use the first 2-3 components to make it fast - and this is enough
# if not, the speed is too slow
# spots on rows

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Error: Not the correct number of arguments!")
}
f_matrix_fn <- args[1]
output_path <- args[2]


f_matrix <- read.table(f_matrix_fn, row.names = 1, header = T, sep = ",")

# Tunable parameters
max_iter <- 200
up_alpha_every <- 10
max_spots <- 500
set.seed(0)
########  sub-sample f_matrix  ############
# DP is slow. So we sub-sample a portion of f_matrix to detect clusters
# but calculate probabilities for all spots

f_matrix <- t(t(f_matrix) - colMeans(f_matrix))
f_matrix_orig <- f_matrix
f_matrix <- f_matrix[
    sample(1:dim(f_matrix)[1], min(dim(f_matrix)[1], max_spots)), ]

########  DP  ##################

dp_cluster <- DirichletProcessMvnormal(f_matrix, numInitialClusters = 2)
cat("Starting dirichlet process clustering...\n")
for (s in seq_len(max_iter)) {
    dp_cluster <- ClusterComponentUpdate(dp_cluster)
    dp_cluster <- ClusterParameterUpdate(dp_cluster)

    if (s %% up_alpha_every == 0) {
        cat(paste("Iteration", s, "\n"))
        dp_cluster <- UpdateAlpha(dp_cluster)
    }
}

########  cluster and cluster prob.  #################

f_matrix <- f_matrix_orig
parameters <- dp_cluster$clusterParameters

# this is the output, probability of each spot belonging to each cluster
cluster_prob <- list()

if (dim(parameters$mu)[3] < 2) {
    cat("Dirichilet process failed to detect >1 clusters!\n")

    km <- kmeans(f_matrix, centers = 2)
    cluster_prob[[1]] <- 1 * (km$cluster == 1)
    cluster_prob[[2]] <- 1 * (km$cluster == 2)
} else {
    plot(dp_cluster)

    for (cluster in 1:dim(parameters$mu)[3])
    {
        mu <- t(f_matrix) - parameters$mu[, , cluster]

        cluster_prob[[cluster]] <-
            exp(-0.5 * apply(mu, 2, function(x) t(x) %*% solve(parameters$sig[, , cluster]) %*% x)) /
                (2 * 3.1415926)^(dim(mu)[1] / 2) /
                determinant(parameters$sig[, , cluster], log = F)$modulus[1]^0.5
    }
}

#########  output  ###############

cluster_prob_matrix <- matrix(NA,
    ncol = length(cluster_prob), nrow =
        length(cluster_prob[[1]])
)
for (i in 1:length(cluster_prob))
{
    cluster_prob_matrix[, i] <- cluster_prob[[i]]
}

rownames(cluster_prob_matrix) <- rownames(f_matrix)
colnames(cluster_prob_matrix) <- paste(
    "cluster", 1:dim(cluster_prob_matrix)[2], sep = "")
cluster_prob_matrix=cluster_prob_matrix/rowSums(cluster_prob_matrix)
write.table(
    cluster_prob_matrix,
    file = paste(output_path, "/dp_cluster_prob.csv", sep = ""),
    quote = F, row.names = T, col.names = T, sep = ",")