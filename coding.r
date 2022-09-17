
setwd('C:\\Users\\Ioanna\\Documents\\MSc in BUSINESS ECONOMICS WITH ANALYTICS\\B SEMESTER\\Econometrics - Ioannis_Dendramis\\2nd')

source('xfunctions.R')


T=5000
yy=matrix(0,T)
er=matrix(0,T)
ss=matrix(0,T)
a0=0.0078
phi=0.0290
c0=0.0013
c1=0.0365
c2=0.9598
for (l in 2:T)
{
  ss[l]=c0+c1*( er[l-1]^2 )+c2*(ss[l-1])
  yy[l]= a0+phi*yy[l-1] +rnorm(1,mean=0,sd=1)*sqrt(ss[l])
  er[l] = yy[l]-a0-phi*yy[l-1]
}
ss1=matrix( ss[101:T],T-100)
yy1=matrix(yy[101:T],T-100)
plot(yy1,type='l', main="Simulated AR(1) data",ylab = 'time series' )
dev.new()
plot(ss1,type='l',main="Simulated GARCH(1,1) data",ylab = 'time series')

thet0=ar1arch1(yy1)
theta0=rbind(thet0,0);
theta0=c(thet0[1:2],var(yy),0,0);
test=likel_ar1_garch(theta0, yy1)
low=c(0.0,0.0,0.001,0.001,0)
h=c(1,1,1,1,1)
theta_opt<-optim(theta0,likel_ar1_garch, gr = NULL,y=yy1,
                 method = c("L-BFGS-B"),
                 lower = low, upper = h,
                 hessian = TRUE)
theta_opt$par
H =theta_opt$hessian
invH<- solve(H)
w=H%*%invH
std=sqrt(diag(invH))
# normal inverse
w=qnorm(1/100)
# For VaR use the theta_opt to compute the mu(T+1), ss(T+1) and the normal inverse
VaR=for_ar1_garch(theta_opt$par, yy1,1/100)
VaR
#theta0=c(a0,phi,0.01,0.01,0.01,0.8)
#theta0=c(0,0.1,var(yy1),0.1,0.1)

thet0=ar1logarch1(yy1)
theta0=rbind(thet0[c(1)],thet0[c(2)],thet0[3],0,0,thet0[4])
test=likel_ar1_egarch(theta0, yy1)
low=c(0,0, -10, -1, -1, 0.0)
h = c(1,1, 10, 1, 1, 0.999)
theta_opt<-optim(theta0,likel_ar1_egarch, gr = NULL,y=yy1,
                 method = c("L-BFGS-B"),
                 lower = low, upper = h,
                 hessian = TRUE)
theta_opt$par
H =theta_opt$hessian
invH<- solve(H)
w=H%*%invH
std=sqrt(diag(invH))
VaR2=for_ar1_egarch (theta_opt$par, yy1,1/100)


#egarch

theta2<-c(0,0,log(var(yy)),0,0,0)
test2=likel_ar1_egarch(theta2,yy1)
low2=c(-5,0, -10, -1, -1, 0.0)#(0,0...)
h2 = c(1,1, 10, 0.1, 1, 0.999)

theta_opt2<-optim(theta2,likel_ar1_egarch,gr = NULL,y=yy1,
                  method = c("L-BFGS-B"),
                  lower = low2, upper = h2,
                  hessian = TRUE)
theta_opt2$par


VaR2=for_ar1_egarch(theta_opt2$par,yy1,1/100)
VaR2
#GJR-GARCH

library(rugarch)
egarchsnp.spec <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1)),mean.model=list(armaOrder=c(0,0))) 
egarchsnp.fit = ugarchfit(egarchsnp.spec,yy1) 
coef(egarchsnp.fit)
library(cvar)
VaR(yy1, p=.99)
