

setwd("C:/Users/Bru/Desktop/metaptuxiako/2o/dedramis/2i ergasia dedramis")

source("manual_functions.R")


likel_ar1_garch = function(theta, y){
  T = length(y);
  ss=matrix(0,T,1);
  ss[2]=var( y );
  lik=0;
  #cat(theta)
  for (l in 2:T)
  {
    er[l]= y[l]- theta[1]-theta[2]*y[l-1];
    #lik=lik+log( norm_pdf1( y[l],theta[1]+theta[2]*y[l-1], abs(ss[l]) ) )
    lik=lik+log_norm_pdf1( y[l],theta[1]+theta[2]*y[l-1], abs(ss[l]) )
    ss[l+1]=theta[3]+theta[4]*(er[l]^2)+theta[5]*ss[l];
  }
  lik1=-lik;
  return(lik1)
}

ar1arch1 <- function(yy) {
  T = length(yy);
  y=matrix( yy[2:T],T-1)
  xo=matrix(yy[1:T-1],T-1)
  c=matrix( rep(1,T-1),T-1)
  x1=cbind(c,xo)
  q1=ols_1(y,x1);
  z=matrix(unlist(q1[1]),2,1)
  er = y-x1%*%z;
  er2=er^2;
  T1 = length(er2);
  w=matrix( er2[2:T1],T1-1);
  wo=matrix(er2[1:T1-1],T1-1);
  wc=matrix( rep(1,T1-1),T1-1);
  w1=cbind(wc,wo);
  q2=ols_1(w,w1);
  z2=matrix(unlist(q2[1]),2,1)
  w=rbind(z,z2);
  return(w)
}


for_ar1_garch = function(theta, y,p){
  T = length(y);
  ss=matrix(0,T+1,1);
  m = matrix(0,T+1,1);
  ss[2]=var( y );
  lik=0;
  #cat(theta)
  for (l in 2:T)
  {
    er[l]= y[l]- theta[1]-theta[2]*y[l-1];
    m[l+1]= theta[1]+theta[2]*y[l];
    ss[l+1]=theta[3]+theta[4]*(er[l]^2)+theta[5]*ss[l];
  }
  var=m[T+1]+sqrt(ss[T+1])*qnorm(p)
  return(var)
}


likel_ar1_egarch=function(theta,y)
{
  T=length(y)
  
  ss=matrix(0,T,1);
  ss[2]=var(y)
  lik=0;
  
  for(l in 2:T)
  {
    er[l]=y[l]-theta[1]-theta[2]*y[l-1];
    
    lik=lik+log_norm_pdf1(y[l],theta[1]+theta[2]*y[l-1],abs(ss[l])) 
    
    ss[l+1]=exp(theta[3]+theta[4]*(log(ss[l]))
                 +theta[5]*(er[l]/sqrt(ss[l]))
                 +theta[6]*((abs(er[l])/sqrt(ss[l]))-sqrt(2/pi)));
    
  }
  lik1=-lik;
  
  return(lik1)
  }


for_ar1_egarch = function(z,y,p)
  
{
  
  T=length(y);
  ss=matrix(0,T+1,1);
  m= matrix(0,T+1,1);
  ss[2]=var(y);
  lik=0;
  
  for (l in 2:T)
  {
    er[l]= y[l]-z[1]-z[2]*y[l-1];
    
    m[l+1]= z[1]+z[2]*y[l];
    
    ss[l+1]=exp(z[3]+z[4]*(log(ss[l]))
                 +z[5]*(er[l]/sqrt(ss[l]))
                 +z[6]*(((er[l])/sqrt(ss[l]))));
    
    
  }
  
  var=m[T+1]+sqrt(ss[T+1])*qnorm(p)
  
  return(var)
  
  
}

