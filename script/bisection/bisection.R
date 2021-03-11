# calculate difference between u and 2^sum(-p_n*log2(p_n))
# https://en.wikipedia.org/wiki/Perplexity
compare_u<-function(u,euc_dist,sigma_n)
{
  tmp=exp(-euc_dist/2/sigma_n^2)
  log_tmp=-euc_dist/2/sigma_n^2
  if (sum(tmp)==0) {return("P0")} # positive or division by 0
  p_n=tmp/sum(tmp)
  log_p_n=log_tmp-log(sum(tmp))
  test=sum(p_n*log_p_n)+log(u)
  if (test>0) {return("P0")} else {return("N")}
}

# find the best sigma for each spot
find_sigma_n<-function(I,n,u,margin)
{
  euc_dist=colSums((I[n,]-t(I[-n,]))^2)
  sigma_n_small=0
  sigma_n_large=sqrt(max(euc_dist)*10)
  
  while (sigma_n_large-sigma_n_small>margin*sqrt(mean(euc_dist)))
  {
    sigma_n_new=(sigma_n_large+sigma_n_small)/2
    results=compare_u(u,euc_dist,sigma_n_new)
    if (results=="N") {sigma_n_large=sigma_n_new} else {sigma_n_small=sigma_n_new}
  }
  
  (sigma_n_large+sigma_n_small)/2
}

# find the sigma for all spots
find_sigma_n_all<-function(I,u,margin)
{
  I=t(I)
  sapply(1:dim(I)[1],function(n) find_sigma_n(I,n,u,margin))
}
