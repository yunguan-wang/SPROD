#include <Rcpp.h>
#include <math.h>  

using namespace Rcpp;

// calculate difference between u and 2^sum(-p_n*log2(p_n))
// https://en.wikipedia.org/wiki/Perplexity
int compare_u(double u,NumericVector euc_dist,double sigma_n) 
{
  int i;
  double sum_tmp,test;
  NumericVector p_n,tmp;
  
  tmp=exp(-euc_dist/(2*sigma_n*sigma_n));

  sum_tmp=0;
  for (i=0;i<tmp.size();++i) {sum_tmp+=tmp[i];}
  if (sum_tmp==0) {return 1;} // positive or division by 0
  
  p_n=tmp/sum_tmp;
  test=log(u);
  for (i=0;i<tmp.size();++i) {test+=(p_n[i]*log(p_n[i]));}
  
  return test>0 ? 1 : 0;
}

// find the best sigma for each spot
// note that n starts from 0 (C++ convention)
double find_sigma_n(NumericMatrix I,int n,double u,double margin)
{
  NumericVector euc_dist(I.ncol()-1),I_n(I.ncol()-1);
  int results,i,j,k,max_euc_dist,ave_euc_dist;
  double tmp,sigma_n_small,sigma_n_large,sigma_n_new;
  
  max_euc_dist=-1;
  ave_euc_dist=0;
   
  for (j=0;j<I.nrow();++j) {I_n[j]=I(j,n);}
   
  for (i=0;i<I.ncol();++i)
  {
    if (i==n) {continue;}
    k=i<n ? i : i-1;
    euc_dist[k]=0;
    
    for (j=0;j<I.nrow();++j)
    {
      tmp=I(j,i)-I_n[j];
      euc_dist[k]+=(tmp*tmp);
    }
  
    if (euc_dist[k]>max_euc_dist) {max_euc_dist=euc_dist[k];}
    ave_euc_dist+=euc_dist[k];
  }
   
  sigma_n_small=0;
  sigma_n_large=sqrt(max_euc_dist*10);
  ave_euc_dist=ave_euc_dist/euc_dist.size();
  
  while (sigma_n_large-sigma_n_small>margin*sqrt(ave_euc_dist))
  {
    sigma_n_new=(sigma_n_large+sigma_n_small)/2;
    results=compare_u(u,euc_dist,sigma_n_new);
    if (results==0) {sigma_n_large=sigma_n_new;} else {sigma_n_small=sigma_n_new;}
  }
  
  return (sigma_n_large+sigma_n_small)/2;
}

// find the sigma for all spots
// [[Rcpp::export]]
NumericVector find_sigma_n_all_Cpp(NumericMatrix I,double u,double margin)
{
  int n;
  NumericVector sigma_n_all(I.ncol());
  
  for (n=0;n<sigma_n_all.size();++n)
  {
    sigma_n_all[n]=find_sigma_n(I,n,u,margin);
  }
  
  return sigma_n_all;
}

