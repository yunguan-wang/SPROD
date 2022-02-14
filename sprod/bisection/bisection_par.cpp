// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>  
#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

// calculate difference between u and 2^sum(-p_n*log2(p_n))
// https://en.wikipedia.org/wiki/Perplexity
int compare_u(double u,arma::vec euc_dist,double sigma_n) 
{
  int i;
  double sum_tmp,test;
  arma::vec p_n,log_p_n,tmp,log_tmp;
  
  log_tmp=-euc_dist/(2*sigma_n*sigma_n);
  tmp=exp(log_tmp);
  
  sum_tmp=0;
  for (i=0;i<tmp.n_rows;++i) {sum_tmp+=tmp[i];}
  if (sum_tmp==0) {return 1;} // positive or division by 0
  
  p_n=tmp/sum_tmp;
  log_p_n=log_tmp-log(sum_tmp);
  
  test=log(u);
  for (i=0;i<tmp.n_rows;++i) {test+=(p_n[i]*log_p_n[i]);}
  
  return test>0 ? 1 : 0;
}

// find the best sigma for each spot
// note that n starts from 0 (C++ convention)
double find_sigma_n(arma::mat I,int n,double u,double margin)
{
  arma::vec euc_dist(I.n_cols-1),I_n(I.n_cols-1);
  int results,i,j,k,max_euc_dist,ave_euc_dist;
  double tmp,sigma_n_small,sigma_n_large,sigma_n_new;
  
  max_euc_dist=-1;
  ave_euc_dist=0;
   
  for (j=0;j<I.n_rows;++j) {I_n[j]=I(j,n);}
   
  for (i=0;i<I.n_cols;++i)
  {
    if (i==n) {continue;}
    k=i<n ? i : i-1;
    euc_dist[k]=0;
    
    for (j=0;j<I.n_rows;++j)
    {
      tmp=I(j,i)-I_n[j];
      euc_dist[k]+=(tmp*tmp);
    }
  
    if (euc_dist[k]>max_euc_dist) {max_euc_dist=euc_dist[k];}
    ave_euc_dist+=euc_dist[k];
  }
   
  sigma_n_small=0;
  sigma_n_large=sqrt(max_euc_dist*10);
  ave_euc_dist=ave_euc_dist/euc_dist.n_rows;
  
  while (sigma_n_large-sigma_n_small>margin*sqrt(ave_euc_dist))
  {
    sigma_n_new=(sigma_n_large+sigma_n_small)/2;
    results=compare_u(u,euc_dist,sigma_n_new);
    if (results==0) {sigma_n_large=sigma_n_new;} else {sigma_n_small=sigma_n_new;}
  }
  
  return (sigma_n_large+sigma_n_small)/2;
}

struct Worker_find_sigma_n : public Worker
{
  // source 
  arma::mat I;
  double u;
  double margin;
  int size;
  IntegerVector ns;

  // destination 
  NumericVector sigma_n_all;
  
  // initialize with source and destination
  Worker_find_sigma_n(const arma::mat I,const double u,
    const double margin,NumericVector sigma_n_all,const int size) 
    : I(I),u(u),margin(margin),sigma_n_all(sigma_n_all),size(size) 
  {
    ns=seq_len(size)-1;
  }
  
  // operator
  void operator()(std::size_t begin, std::size_t end) 
  {
    std::transform(ns.begin() + begin, ns.begin() + end, 
      sigma_n_all.begin() + begin, 
      [&](int n) -> double {return find_sigma_n(I,n,u,margin);});
  }
};

// find the sigma for all spots
// [[Rcpp::export]]
NumericVector find_sigma_n_all(arma::mat I,double u,double margin)
{
  int size=I.n_cols;
  NumericVector sigma_n_all(size); // output
  Worker_find_sigma_n worker_find_sigma_n(I,u,margin,sigma_n_all,size);

  parallelFor(0,size,worker_find_sigma_n);

  return sigma_n_all;
}
