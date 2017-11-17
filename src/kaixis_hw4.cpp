#include <Rcpp.h>
using namespace std;
using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector ctmcViterbi(NumericVector ts, double theta, NumericMatrix obs)
{
  int m = (int)ts.size();
  int n = (int)obs.nrow();

  NumericMatrix obs1(n,m);
  
  double max1=-(1/0.0);

  double index_max1=0;
  double max2=-(1/0.0);
  double index_max2=0;
  
  double diag_t=0;
  double off_t=0;
  
  for(int i=0;i<n;++i){
    for(int j=0;j<m;++j){
      obs1(i,j)=log(obs(i,j));
    }

  }

  NumericMatrix delta(n,m);
  NumericVector delta_1(n);
  NumericVector delta_2(n);
  NumericMatrix phi(n,m);
  
  IntegerVector viterbiPath(m);
  
  
  if ( obs1.ncol() != m ){
    stop("The input matrix does not conform to the other parameters");
  }
  
  for (int i=0; i<n;++i){
    delta(i,0)=log(pow(n,-1))+obs1(i,0);
    
    delta_1(i)=delta(i,0);
    delta_2(i)=delta(i,0);
    
  }
  
  for(int i=0;i<n;++i){
    
    if(max1<delta_1(i)){
      
      max1=delta_1(i);
      index_max1=i;
    }
  }
  
  for (int i=0;i<n;++i){
    if(i!=index_max1){
      
      if(max2<delta_2(i)){
        max2=delta_2(i);
        index_max2=i;
      }
      
    }
  }
  
  
  
  for(int t=1;t<m;t++){
    diag_t=log(1-((n-1)*pow(n,-1))*(1-exp(-1*theta*(ts[t]-ts[t-1]))));
    off_t=log((pow(n,-1))*(1-exp(-1*theta*(ts[t]-ts[t-1]))));
    
    for(int i=0;i<n;i++){
      if(i==index_max1){
        
        if((max1+diag_t)>(max2+off_t)){
          delta(i,t)=(max1+diag_t);
          phi(i,t)=index_max1;
        }
        else{
          delta(i,t)=(max2+off_t);
          phi(i,t)=index_max2;
        }
        
      }
      else{
        if((max1+off_t)>(delta(i,t-1)+diag_t)){
          delta(i,t)=(max1+off_t);
          phi(i,t)=index_max1;
        }
        else{
          delta(i,t)=(delta(i,t-1)+diag_t);
          phi(i,t)=i;
        }
        
      }
      
      delta(i,t)+=obs1(i,t);
      delta_1(i)=delta(i,t);
      delta_2(i)=delta(i,t);
    }
    
    max1=-(1/0.0);
    max2=-(1/0.0);
    
    
    for(int i=0;i<n;++i){
      if(max1<delta_1(i)){
        
        max1=delta_1(i);
        index_max1=i;
      }
    }
    
    for (int i=0;i<n;++i){
      if(i!=index_max1){
        
        if(max2<delta_2(i)){
          max2=delta_2(i);
          index_max2=i;
        }
        
      }
    }
    
    
  }
  
  double ml = -(1/0.0);
  for(int i=0; i < n; ++i) {
    if ( ml < delta(i,m-1) ) {
      ml = delta(i,m-1);
      viterbiPath[m-1] = i;
    }
  }  
  for(int i=m-1; i > 0; --i) {
    viterbiPath[i-1] = phi(viterbiPath[i],i);
  }
  
  
  
  // TODO : Implement your function here
  return viterbiPath;
}



// [[Rcpp::export]]
void backwardLoop(NumericMatrix& obs, NumericMatrix& beta, NumericVector& ts, double theta){
  int T=(int)ts.size();
  int n = (int)obs.nrow();
  
  double diag=0;
  double off=0;
  
  NumericVector heng_diag(T);
  NumericVector heng_off(T);
  
  
  for(int i=0;i<n;++i){
    beta(i,T-1)=1;
  }
  
  
  
  
  
  for(int t=T-2;t>=0;--t){
    diag=1-((n-1)*pow(n,-1))*(1-exp((-1)*theta*(ts[t+1]-ts[t])));
    off=(pow(n,-1))*(1-exp(-1*theta*(ts[t+1]-ts[t])));
    for(int i=0;i<n;++i){
      heng_off(t)+=beta(i,t+1)*obs(i,t+1)*off;
    }
    
    for(int i=0;i<n;++i){
      beta(i,t)=heng_off(t)+beta(i,t+1)*obs(i,t+1)*(diag-off);
    }
    
    
  }
 
  
}



// [[Rcpp::export]]
void forwardLoop( NumericMatrix& obs,NumericMatrix& alpha, NumericVector& ts, double theta){
  
  int T=(int)ts.size();
  int n = (int)obs.nrow();
  
  NumericVector heng(T);
  
  
  double diag=0;
  double off=0;
  
  
  
  for(int i=0;i<n;++i){
    alpha(i,0)=pow(n,-1)*obs(i,0);
    
  }
  // for(int i=0;i<n;++i){
  //   heng(0)+=alpha(i,0);
  // }
  
  for(int t=1;t<T;++t){
    diag=1-((n-1)*pow(n,-1))*(1-exp(-1*theta*(ts[t]-ts[t-1])));
    off=(pow(n,-1))*(1-exp(-1*theta*(ts[t]-ts[t-1])));
    
    
    
    for(int i=0;i<n;++i){
      
      heng(t)+=(alpha(i,t-1)*off);
      
    }
    
    for (int i=0;i<n;++i){
      alpha(i,t)=heng(t)+(alpha(i,t-1)*(diag-off));
      alpha(i,t)*=obs(i,t);
    }
    
  }
  

}


// Do not forget to add documentation for your package using roxygen2
// [[Rcpp::export]]
NumericVector ctmcForwardBackward(NumericVector ts, double theta, NumericMatrix obs) {
  int m = (int)ts.size();
  int n = (int)obs.nrow();
  if ( obs.ncol() != m ) {
    stop("The input matrix does not conform to the other parameters");
  }
  
  NumericMatrix alpha(n,m);
  NumericMatrix beta(n,m);
  
  NumericVector pi(m);
  NumericMatrix condProb(n,m);
  
  
  forwardLoop(obs, alpha, ts, theta);
  backwardLoop(obs, beta, ts, theta);
  
  
  for(int t=0;t<m;++t){
    double sum=0;
    for(int i=0;i<n;++i){
      sum+=(alpha(i,t)*beta(i,t));
    }
    for(int i=0;i<n;++i){
      condProb(i,t)=alpha(i,t)*beta(i,t)/sum;
    }
    
  }
  
  
  
  
  // TODO : Implement your function here
  return condProb;
}