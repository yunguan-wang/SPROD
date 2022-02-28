calculate_spot_dist <- function(n2, sigma_n, euc_dist2) 
{
    # calculate spot distance based on gaussian kernal regulated euclidean distance
    # n1, n2 is indicies, sigma_n is a vector of N and euc_distance is a matrix of N x N.
    # will return a vector of N, which is P_n_n2
    tmp = -euc_dist2[,n2]/2/sigma_n[n2]^2
    tmp = tmp - max(tmp[-n2])
    d_n_n2 = exp(tmp)
    p_n_n2 =d_n_n2/sum(d_n_n2[-n2])
    p_n_n2[n2] = 0
    return(p_n_n2)
}

# The following functions are just wrappers to be used in optim.
fr <- function(par, data) 
{
    # take parameters
    alpha=matrix(par,ncol=ncol(data[[3]]))
    alpha=alpha+t(alpha)
    
    lambda=data[[1]]
    p_nn_tsne=data[[2]]
    proximity=data[[3]]
    k=data[[4]]
    
    # Q
    Q = diag(1, dim(alpha)[1])+4*diag(rowSums(alpha))-4*alpha
    
    # cost
    cost=sum(alpha*lambda*p_nn_tsne)- 
        k/2*determinant(Q, logarithm=TRUE)$modulus[1]
    #cat(paste("cost=",round(cost,d=10),"\n",sep=""))
    cost
}

gr <- function(par, data) 
{
    # take parameters
    alpha=matrix(par,ncol=ncol(data[[3]]))
    alpha=alpha+t(alpha)
    
    lambda=data[[1]]
    p_nn_tsne=data[[2]]
    proximity=data[[3]]
    k=data[[4]]
    
    # Q^-1
    Q = diag(1, dim(alpha)[1])+4*diag(rowSums(alpha))-4*alpha
    Q_inv=solve(Q)
    Q_inv=(Q_inv+t(Q_inv))/2
    
    tmp=Q_inv
    tmp[]=diag(Q_inv)
    
    # gradient
    #print(mean(abs(lambda*p_nn_tsne))/mean(abs(2*k*(tmp+t(tmp)-2*Q_inv))))
    lambda*p_nn_tsne-2*k*(tmp+t(tmp)-2*Q_inv)   
}


