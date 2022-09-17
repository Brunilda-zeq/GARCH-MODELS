

setwd("C:/Users/Bru/Desktop/metaptuxiako/2o/dedramis/2i ergasia dedramis")

autocorr_1<- function(y,l) {
  T<-length(y)
  mm<-mean(y)
  yy<- y-matrix(rep(mm,T),T)
  v<-matrix(yy[c((l+1):T)],T-l)
  v_l<- matrix(yy[1:(T-l)],T-l)
  g_l<-crossprod(v,v_l)/(T-l)
  g_0<-var(y)
  rho_l<- g_l/g_0
  w<- list(rho_l,g_l)
  return(w)
}

ols_1<-function(y,x) {
  #Inverse matrix of x'x
  z1<-solve(t(x)%*%x)
  b<- z1%*%(crossprod(x,y))
  e<- y-x%*%b
  T<-length(e)
  k<- ncol(x)
  ss<- c(crossprod(e,e)/(T-k))
  s1<-z1*ss
  std<-matrix(sqrt(diag(s1)),k)
  w<-list(b,std)
  return(w)}

likel_ar1<- function(theta,y){
  T<-length(y);
  lik<-0;
  for (l in 2:T) {
    lik<-lik +log(norm_pdf1(y[l],theta[1]+theta[2]*y[l-1],theta[3]))
    
  }
  lik1=-lik;
  return(lik1)
}

norm_pdf1<- function(x,m,ss) {
  #normal pdf
  f1<- ((2*pi*ss)^(-1/2))
  f2<- exp(-((x-m)^2)/(2*ss))
  f<-f1*f2;
  return(f)
}

log_norm_pdf1<- function(x,m,ss) {
  #normal pdf
  f1<-(-1/2)*log(ss);
  f2= -((x-m)^2)/(2*ss);
  f3= -(1/2)*log(2*pi);
  f<- f1+f2+f3;
  return(f)
}



likel_ar2<- function(theta,y){
  T<-length(y);
  lik<-0;
  for (l in 3:T) {
    lik<-lik +log(norm_pdf1(y[l],theta[1]+theta[2]*y[l-1]+theta[3]*y[l-2],theta[4]))
    
  }
  lik2=-lik;
  return(lik2)
}


likel_ar3<- function(theta,y){
  T<-length(y);
  lik<-0;
  for (l in 4:T) {
    lik<-lik +log(norm_pdf1(y[l],theta[1]+theta[2]*y[l-1]+theta[3]*y[l-2]+theta[4]*y[l-3],theta[5]))
    
  }
  lik3=-lik;
  return(lik3)
}

