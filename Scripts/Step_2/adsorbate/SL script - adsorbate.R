library("R.matlab")
library("glmnet")
library("LaplacesDemon")
library("GIGrvg")
library("MASS")
library("scales")
library("mvtnorm")
library("expm")

set.seed(1)

DL<-function(x,y){
  p=ncol(x)
  n=nrow(x)
  #calculate hyperparameter
  xtx=t(x)%*%x
  d=eigen(xtx/n)$values
  P=sum(d)
  Q=4*sum(d^2)-sum(d)^2
  R=-sum(d)^3
  C=P^2/9-Q/3
  A=P*Q/6-P^3/27-R/2
  B=A^2-C^3
  hyper=sqrt(2/((A+sqrt(B))^(1/3)+sign(A-sqrt(B))*abs(A-sqrt(B))^(1/3)-P/3))
  hyper[hyper<1/p]=1/p
  #initial parameters
  a=rep(hyper,p)
  psi=rexp(p,rate=1/2)
  psi1=rep(0,p)
  phi=rdirichlet(n=1,alpha=a)
  phi[phi <= (1e-40)]<-(1e-40)
  tau=rgamma(n=1,shape=p*a,rate=1/2)
  Ti=rep(1,p)
  beta=rep(0,p)
  hi=rep(1,p)
  betamatrix<-matrix(rep(NA,5000*p),nrow=5000)
  
  #Niter iterations
  for(i in 1:10000){
    #step1:sample sigma^2
    s=c(psi*phi^2*tau^2)
    E_1=max(t(y-x%*%beta)%*%(y-x%*%beta),1e-8)
    E_2=max(sum(beta^2*s),1e-8)
    sigma2=1/stats::rgamma(1,(n+p)/2,rate=(E_1+E_2)/2)
    sigma1=sqrt(sigma2)
    if(sigma1>1e20) print("Please choose a better hyperparameter, it is too big")
    #step2:sample beta
    u=rnorm(p)*sqrt(s)
    delta=rnorm(n)
    v=x%*%u+delta
    stx=as.numeric(s)*t(x)
    w=ginv(x%*%stx+diag(n))%*%(y/sigma1-v)
    beta=(u+(stx%*%w))*sigma1
    
    mix1=abs(beta)/sigma1
    mix2=mix1/c(phi)
    #step3:sample psi
    mu=tau/mix2
    pv=(rnorm(p))^2
    pu=runif(p)
    temp2=mu*pv
    temp3=sqrt(4*temp2+temp2^2)
    temp4=mu+0.5*(pv*(mu^2))-0.5*(mu*temp3)
    locs=(pu<=mu/(mu+temp4))
    psi1=locs*temp4+(1-locs)*(mu^2/temp4)
    psi=1/psi1
    
    #step4:sample tau
    tau=GIGrvg::rgig(n=1,lambda=p*a-p,psi=1,chi=2*sum(mix2))
    
    #step5:sample phi
    hu=runif(p,0,exp(-1/(2*hi)))
    hl=1/(2*log(1/hu))
    hf=pgamma(hl,shape=1-a,rate=mix1)
    hr=pmin(runif(p,hf,1),(1-(1e-20)))
    hi=qgamma(hr,shape=1-a,rate=mix1)
    Ti=1/hi
    phi=Ti/sum(Ti)
    phi[phi<=(1e-20)]=(1e-20)
    
    #beta output
    if(i>5000) betamatrix[i-5000,]=beta
    if(i%%1000==0) print(i)
  }
  return(betamatrix)
}


hs<-function(x,y)
{
  n=nrow(x)
  p=ncol(x)
  ## parameters ##
  Beta=rep(0,p)
  lambda=rep(1,p)
  sigma_sq=1
  sigma1=1
  tau=1
  U=matrix(rep(0,p*n),nrow=p)
  ## output ##
  betaout=matrix(rep(0,p*5000),nrow=5000)
  
  ## matrices ##
  l0=rep(0,p)
  l1=rep(1,n)
  l2=rep(1,p)
  ## start Gibb's sampling ##
  for(i in 1:10000)
  {
    ## update beta ##
    lambda_star=tau*lambda
    U=as.numeric(lambda_star^2)*t(x)
    ## step 1 ##
    u=stats::rnorm(l2,l0,lambda_star)
    v=x%*%u+stats::rnorm(n)
    ## step 2 ##
    sigma1=sqrt(sigma_sq)
    v_star=ginv(x%*%U+diag(n))%*%((y/sigma1)-v)
    Beta=sigma1*(u+U%*%v_star)
    
    
    ## update lambda_j's in a block using slice sampling ##
    eta = 1/(lambda^2)
    upsi = stats::runif(p,0,1/(1+eta))
    tempps = Beta^2/(2*sigma_sq*tau^2)
    ub = (1-upsi)/upsi
    # now sample eta from exp(tempv) truncated between 0 & upsi/(1-upsi)
    Fub = 1 - exp(-tempps*ub) # exp cdf at ub
    Fub[Fub < (1e-4)] = 1e-4;  # for numerical stability
    up = stats::runif(p,0,Fub)
    eta = -log(1-up)/tempps
    lambda = 1/sqrt(eta)
    
    ## update tau ##
    ## Only if prior on tau is used
    tempt = sum((Beta/lambda)^2)/(2*sigma_sq)
    et = 1/tau^2
    utau = stats::runif(1,0,1/(1+et))
    ubt = (1-utau)/utau
    Fubt = stats::pgamma(ubt,(p+1)/2,scale=1/tempt)
    Fubt = max(Fubt,1e-8) # for numerical stability
    ut = stats::runif(1,0,Fubt)
    et = stats::qgamma(ut,(p+1)/2,scale=1/tempt)
    tau = 1/sqrt(et)
    
    
    ## update sigma_sq ##
    E_1=max(t(y-x%*%Beta)%*%(y-x%*%Beta),(1e-8))
    E_2=max(sum(Beta^2/((tau*lambda))^2),(1e-8))
    sigma_sq=1/stats::rgamma(1,(n+p)/2,scale=2/(E_1+E_2))
    
    if(i%%1000==0)  {print(i)}
    if(i>5000) betaout[i-5000,]=Beta
  }
  return(betaout)
}

LS<-function(x,y){
  p=ncol(x)
  xc<-scale(x, scale=F)
  yc<-y-mean(y)
  model<-glmnet(xc,yc,standardize=FALSE,alpha=1,family="gaussian")
  lam=cv.glmnet(xc,yc,standardize=FALSE,alpha=1,family="gaussian")$lambda.1se
  beta=coef(model,s=lam)[2:(p+1)]
  return(beta)
}


countnormalize<-function(betamatrix){
  p=ncol(betamatrix)
  betamean<-apply(betamatrix,2,mean)
  betacov<-cov(betamatrix)
  D<-betamean^2
  cov1=expm::sqrtm(betacov)
  rm(betacov)
  covinv<-ginv(cov1)
  rm(cov1)
  xstar=t(as.numeric(D)*covinv)
  ystar=covinv%*%betamean
  #scale 
  xc<-scale(xstar,scale=T,center=T)
  yc<-ystar-mean(ystar)
  sc=attributes(xc)$"scaled:scale"
  model<-glmnet(xc,yc,standardize=F,alpha=1,family="gaussian")
  lam=cv.glmnet(xc,yc,standardize=F,alpha=1,family="gaussian")$lambda.1se
  betastar=c(coef(model,s=lam)[2:(p+1)])/sc
  betatil=D*betastar
  return(betatil)
}

countdirect<-function(result){
  #use penalized credible region to do variable selection
  #matrix computation to make the solutions be accomplished by LASSO
  p=ncol(result)
  betamean<-apply(result,2,mean)
  betacov<-cov(result)
  D<-betamean^2
  cov1=expm::sqrtm(betacov)
  covinv<-MASS::ginv(cov1)
  xstar=t(as.numeric(D)*covinv)
  ystar=covinv%*%betamean
  #solve LASSO problem by glmnet package
  model<-glmnet::glmnet(xstar,ystar,standardize=FALSE,alpha=1,family="gaussian")
  lam=glmnet::cv.glmnet(xstar,ystar,standardize=FALSE,alpha=1,family="gaussian")$lambda.1se
  betastar=coef(model,s=lam)[2:(p+1)]
  betatil=D*betastar
  return(betatil)
}



datapreprocessing<-function(x1,y1){
len1=ncol(x1)
#remove x1 cols with duplicated data
namex1=as.character(1:len1)
datax1=data.frame(t(x1))
x11=x1[,!duplicated(datax1)]
namex11=namex1[!duplicated(datax1)]
len11=ncol(x11)


#remove x without variability
s=rep(0,len11)
for(i in 1:len11){
  s[i]=sd(x11[,i])
}
x111=x11[,which(s!=0)]
namex111=namex11[which(s!=0)]
len111=ncol(x111)



#check the categorical property
num=rep(0,len111)
#categorical x and factor
for(i in 1:len111){
  a=unique(x111[,i])
  num[i]=length(a)
}



#xmat divide into two
xmat11=x111[,which(num>4)]
xmat2=x111[,which(num<=4)]
namexmat11=namex111[which(num>4)]
namexmat2=namex111[which(num<=4)]



lenxmat11=ncol(xmat11)
#multiplication
xx1=matrix(rep(0,nrow(xmat11)*lenxmat11*(lenxmat11+1)/2),nrow=nrow(xmat11))
xx1[,(1:lenxmat11)]=xmat11
namexx1=as.character(rep(0,lenxmat11*(lenxmat11+1)/2))
namexx1[(1:lenxmat11)]=namexmat11
for(i in 1:(lenxmat11-1))
  for(j in (i+1):lenxmat11){
    xx1[,((2*lenxmat11-i-1)*i/2+j)]=xx1[,i]*xx1[,j]
    namexx1[((2*lenxmat11-i-1)*i/2+j)]=paste0(namexx1[i],"*",namexx1[j])
  }
lenxx1=ncol(xx1)
s=rep(0,lenxx1)
for(i in 1:lenxx1){
  s[i]=sd(xx1[,i])
}
xmat1=xx1[,which(s!=0)]
lenxmat1=ncol(xmat1)
namexmat1=namexx1[which(s!=0)]


# to get rid of redundant categorical variable 
xmat21=matrix(rep(NA,nrow(xmat2)*ncol(xmat2)),nrow=nrow(xmat2))
for(i in 1:ncol(xmat2))
{
  n=xmat2[,i]
  a=unique(n)
  for(j in 1:length(a)){
    xmat21[which(n==a[j]),i]<-letters[j]
  }
}
datamat21=data.frame(t(xmat21))
mat2=xmat21[,!duplicated(datamat21)]
namemat2=namexmat2[!duplicated(datamat21)]
colnames(mat2)=namemat2


#generate design matrix with (k-1) dummy variables
mat2f=model.matrix(~.,data.frame(mat2))
mat2fs=mat2f[,-1]
namemat2fs=attributes(mat2fs)$"dimnames"[[2]]



#screen down by correlation
lenxmat1=ncol(xmat1)
co=rep(0,lenxmat1)
co=cor(xmat1,y1)
or=order(abs(co))
xmat1d=xmat1[,which(or<=1000)]
namexmat1ds=namexmat1[which(or<=1000)]



#scale xmat1d
xmat1ds<-scale(xmat1d,scale=T,center=T)
sc=attributes(xmat1ds)$"scaled:scale"
ce=attributes(xmat1ds)$"scaled:center"

#combine the matrix
xmat=cbind(xmat1ds,mat2fs)
namexmat=c(namexmat1ds,namemat2fs)


xmatn=a<-data.frame(name=namexmat,mat=t(xmat))
return(xmatn)
}



## Adsorbate-modified MgO ##

#transfer
data<-readMat("../../Step_1/adsorbate/data_set_adsorbate.mat")

#save
save(data,file="dataapplication2.RData")

load("dataapplication2.RData")

x1=data$Data.A.p.train
y1=data$Data.A.train

#data$Data.A.p.train [1:52, 1:4823]
#data$Data.A.train [1:52, 1]

#data preprocessing
xmatn=datapreprocessing(x1,y1)
namexmat=as.character(xmatn$name)
xmat=t(xmatn[,2:ncol(xmatn)])
yc<-y1-mean(y1)


#laaso
set.seed(20191)
model<-glmnet(xmat,yc,standardize=FALSE,alpha=1,family="gaussian")
lam=cv.glmnet(xmat,yc,standardize=FALSE,alpha=1,family="gaussian")$lambda.1se
lsbeta=coef(model,s=lam)
lsbeta[which(lsbeta!=0)]
lsname=namexmat[which(lsbeta!=0)]
write.table(lsname, "./MgO_adsorbate_LS.txt", sep="\t", row.names = F, col.names = F)
lsname

#horseshoe sampling
set.seed(20192)
hsresult=hs(xmat,yc)

#horsesshoe without normailize
hsbeta1=countdirect(hsresult)
write.table(namexmat[which(hsbeta1!=0)], "./MgO_adsorbate_HS.txt", sep="\t", row.names = F, col.names = F)
namexmat[which(hsbeta1!=0)]



#dirichlet laplace
set.seed(20193)
dlresult=DL(xmat,yc)

#DL without normalization
dlbeta1=countdirect(dlresult)
write.table(namexmat[which(dlbeta1!=0)], "./MgO_adsorbate_DL.txt", sep="\t", row.names = F, col.names = F)
namexmat[which(dlbeta1!=0)]

