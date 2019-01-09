#-------------------------------------------------------------------------------------#
###################################    Funciones  #####################################
#-------------------------------------------------------------------------------------#

# 1er y 2do Momentos de una normal truncada -------------------------------------------
# u  = vector media
# S  = Matris de varianza
# qc = datos

MomemNT <- function(u=c(0,0),S=diag(2),qc=c(1,2)) {
  
  nic=length(u)
  
  if (nic==1) {
    
    qq <- (1/sqrt(S))*(-qc+u)
    R<-1
    alpha <- pnorm(-qq)
    
    dd <- dnorm(-qq)
    H <- qq*dd
    EX <- (1/alpha)*dd   # a vector with a length of nic
    EXX <- 1+1/alpha*H
    varX <- EXX-EX^2
    Eycens <- -sqrt(S)*EX+u
    varyic<- varX*S
    E2yy<-varyic+Eycens^2
    
  }
  
  else {
    
    qq <- diag(1/sqrt(diag(S)))%*%(-qc+u)
    R <-  diag(1/sqrt(diag(S)))%*%S%*%diag(1/sqrt(diag(S)))
    alpha <- pmvnorm(upper=as.vector(-qq), corr=R)
    #print(qq)
    dd <- rep(0, nic)   #derivative vector
    
    for (j in 1:nic){
      V <- R[-j, -j, drop=F]-R[-j,j, drop=F]%*%R[j,-j, drop=F]
      nu <- -qq[-j]+R[-j,j, drop=F]%*%qq[j]
      dd[j] <- dnorm(-qq[j])*pmvnorm(upper=as.vector(nu), sigma=V)
    }
    
    H <- matrix(rep(0, nic*nic), nrow=nic)
    RH <- matrix(rep(0, nic*nic), nrow=nic)
    
    if(nic==2)     {
      H[1,2] <- H[2,1] <- dmvnorm(-qq[c(1, 2)],sigma=matrix(c(1, R[1,2], R[2,1], 1), nrow=2))
      #sigma==R since qq is standardized
      RH[1,2] <- RH[2,1] <- R[1,2]*H[1,2]
    }
    
    else {
      for( s in 1:(nic-1)){
        for (t in (s+1):nic){
          invR <- solve(R[c(s,t), c(s,t), drop=F])
          nu <- -qq[-c(s,t)]+R[-c(s,t), c(s,t), drop=F]%*%invR%*%qq[c(s,t),,drop=F]
          V <-  R[-c(s,t), -c(s,t), drop=F]- R[-c(s,t), c(s,t), drop=F]%*%invR%*%R[c(s,t), -c(s,t), drop=F]
          H[s,t] <- H[t,s] <- pmvnorm(upper=as.vector(nu), sigma=V)*dmvnorm(-qq[c(s, t)],sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2))
          RH[s,t] <- RH[t,s] <- R[s,t]*H[s,t]
        }
      }
    }
    
    h <- qq*dd-apply(RH, 1, sum)
    diag(H) <- h
    EX <- (1/alpha)*R%*%dd   # a vector with a length of nic
    EXX <- R+1/alpha*R%*%H%*%R
    varX <- EXX-EX%*%t(EX)
    Eycens <- -diag(sqrt(diag(S)))%*%EX+u
    varyic <- diag(sqrt(diag(S)))%*%varX%*%diag(sqrt(diag(S)))
    E2yy <- varyic+Eycens%*%t(Eycens)
    
  }
  
  return(list(Ey=Eycens,Eyy=E2yy,Vary=varyic))
  
}

# Construcion de los phi --------------------------------------------------------------
# pit = pi's

estphit = function(pit){
  p = length(pit)
  Phi = matrix(0,ncol=p,nrow=p)
  if (p>1) {
    diag(Phi) = pit
    for (j in 2:p) {
      for (k in 1:(j-1)) {
        Phi[j,k] = Phi[j-1,k] - pit[j]*Phi[j-1,j-k]
      }
    }
    return(Phi[p,])
  }
  else return(pit)
}



# Construcion de la Matriz Mn(phi) ----------------------------------------------------
# pi     = pi's
# n      = tamnho
# sigma2 = varianza del error

MatArp<-function(pii,n,sigma2){
  p = length(pii)
  phi=estphit(pii) 
  if (n==1) Rn = 1
  else Rn = toeplitz(ARMAacf(ar=phi, ma=0, lag.max = n-1))
  rhos = ARMAacf(ar=phi, ma=0, lag.max = p)[(1:p)+1]
  Rnx<-sigma2*Rn/(1-sum(rhos*phi))
  Rnx<-(Rnx+t(Rnx))/2
  return(Rnx)
}




# Funcion para maximizar el EM : Estima phis de Arp -----------------------------------
# pi     = vector de pi's
# beta1  = vector de betas
# sigmae = varianza del error
# ubi    = Esperana de b (efecto misto)
# ubbi   = Esperana de bb'
# uybi   = Experanza de yb
# uyyi   = Esperanza de yy'
# uyi    = esperanza de y'
# x      = matrix corespondiente de los efectos fijos beta
# z      = matrix corespondiente de los efectos aletorios b
# nj     = observaciones por individuo

FCiArp<-function(pi,beta1,sigmae, ubi,ubbi,uybi,uyyi,uyi,x,z,nj){
  m<-length(nj)[1]
  N<-sum(nj)
  p<-dim(x)[2]
  q1<-dim(z)[2]
  m1<-m*p
  m2<-m*q1
  
  soma=0
  for (j in 1:m ){
    
    ub<-ubi[(((j-1)*q1)+1) : (j*q1), j]
    ubb<-ubbi[(((j-1)*q1)+1) : (j*q1), (((j-1)*q1)+1) : (j*q1)]
    uyb<-uybi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(((j-1)*q1)+1) : (j*q1)]
    uyy<-uyyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
    uy<-uyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),j]
    x1=matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
    z1=matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
    #tt1=tt[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ]
    gammai=x1%*%beta1                                                
    Cii<- MatArp(pi,nj[j],sigmae) #MatAr1(tt1,rho,gamma,sigmae)
    
    soma<- soma - 0.5*log(det(as.matrix(Cii)))-0.5*(sum(diag(uyy%*%solve(Cii)))-t(uy)%*%solve(Cii)%*%gammai-t(gammai)%*%solve(Cii)%*%uy-sum(diag(solve(Cii)%*%((uyb)%*%t(z1))))-sum(diag(solve(Cii)%*%((uyb)%*%t(z1))))
                                                    +t(gammai)%*%solve(Cii)%*%z1%*%ub+t(ub)%*%t(z1)%*%solve(Cii)%*%gammai+t(gammai)%*%solve(Cii)%*%gammai+sum(diag(ubb%*%t(z1)%*%solve(Cii)%*%z1)))
  }
  
  return(-soma)
}

##########################################################################
## Derivadas da matriz do modelo Yi=XiB+Zbi+ei, ei~N(0,sigma2*Mn(phi)), ##
##  bi~N(0, D(alpha)), vector de parametros teta=(betas,phi,sigma,alfa) ##
##########################################################################


# Funcion para el calculo de D(y,beta) -------------------------------
# Matriz usada para la reparametrizacion del modelo 
# en que :
# D[i,j] = (yi-xi*beta)*(yj-xj*beta)+...+(yn+1-j-xn+1-j*beta)*(yn+1-i-xn+1-i*beta)
# beta   = betas do modelo
# y      = vector y (variable independediente)
# Z      = matriz Z (variables correspondiente a los afectos aleatorios)
# b      = efectos aleatorios
# x      = Matriz X (variables dependedientes)
# p      = orden del autoregresivo

Dbeta = function(beta,y,x,z,b,p) {
  n = length(y)
  D = matrix(0,p+1,p+1)
  for (ii in 1:(p+1)) {
    for (jj in 1:(p+1)) {
      D[ii,jj] = sum((y-x%*%beta-z%*%b)[ii:(n+1-jj)]*(y-x%*%beta-z%*%b)[jj:(n+1-ii)])
      D[is.na(D)]<-0
    }
  }
  return(D)
}

# FUNCION DE LA PRIMERA DERIVADA DE LA MATRIZ Mn(phi) (Errores auto correlacionados)------------------------------------
# (vector de derivaras de phi)
# se para la deriabada de phi se usa la matriz particionada D[D11,Dphi1;D1phi,Dphiphi]
#lamda*D(y,beta)*landa= D11-2*phi*Dphi1+2*phi*Dphiphi*phi 
# en que landa=(1,phi)
# M = Matriz de variaza de los efectos aleaoreos 

#Matriz gp=log(det(Mp(phi)))
#gp = function(phi,sigma=sig2){if(length(phi)==1){log(MatArp(phi,length(phi),sigma))
#}else{log(det(MatArp(phi,length(phi),sigma)))}}


D11   = function(beta,y,x,z,b,p) matrix(Dbeta(beta,y,x,z,b,p)[1,1])
Dphi = function(beta,y,x,z,b,p) matrix(Dbeta(beta,y,x,z,b,p)[2:(p+1),1])
Dphiphi = function(beta,y,x,z,b,p) Dbeta(beta,y,x,z,b,p)[2:(p+1),2:(p+1)]



# FUNCION DE LA PRIMERA DERIVADA DE LA MATRIZ D (efectos aleatoreos)------------------------------------
# (vector de derivaras de alpha)
# M = Matriz de variaza de los efectos aleaoreos 

dD<-function(M){
  
  m1<-dim(M)[1]
  m2<-dim(M)[2]  
  d<-list()
  for(i in 1:m1){
    d[[i]]<-list()
    for(j in 1:(m2+1-i)){
      d[[i]][[j]]<-matrix(0,m1,m2)
      if(j==1){d[[i]][[j]][i,i]<-1}   
      else{
        d[[i]][[j]][i,i+(j-1)]<-d[[i]][[j]][i+(j-1),i]<-1}
    }
  }
  
  return(d=d)
  
}

#MATRIS DE INFORACION PARA MI MODELO----------------------------------------

#MAtris jacobiana (primeras derivadas)
# theta = parametros del modelo ordenado de la siguitente forma (betas,sigma,phis)
# yest     = vector y (variable independediente) uyi
# x     = matriz X (variables dependedientes)
# Z    = matriz Z (variables correspondiente a los afectos aleatorios)
# best     = Ubi
# bb    =ubb
# D1= matriz D

#theta=tetaMI;y=uy;b=ub;bb=ubb; x=x1;z=z1

Jt = function(theta,y,x,z,b,bb,p,Arp,D1) {
  l      = p
  n      = length(y)
  beta   = matrix(theta[1:l])
  sig2   = theta[l+1]
  
  if(Arp==0){
    Mn     = diag(1,n)*sig2
    invMn  = solve(Mn/sig2)
    spi= (t(y-x%*%beta-z%*%b)%*%invMn%*%(y-x%*%beta-z%*%b))
          }
  if(Arp!=0){  phi    = theta[(l+2):length(theta)]
  p      = length(phi)
  Mn     = MatArp(phi,n,sig2)
  invMn  = solve(Mn/sig2)
  lambda = matrix(c(-1,phi))
  spi    = t(lambda)%*%Dbeta(beta,y,x,z,b,p)%*%lambda}
  # if(length(invMn)==1){
  #   dbeta  = 1/sig2*(x%*%invMn%*%(y-z%*%b)- x%*%invMn%*%t(x)%*%beta)
  #   }
  #if(length(invMn)!=1){
    dbeta  = 1/sig2*(t(x)%*%invMn%*%(y-z%*%b)- t(x)%*%invMn%*%x%*%beta)
   # }
  
  dsig2  = -n/(2*sig2) +spi/(2*sig2^2)
  if(length(D1)==1){
    dD_alp = 1
    dalpha<-rep(0,1)
    md2<-1
  }
  if(length(D1)>1){
    dD_alp = dD(D1) 
    dalpha<-rep(0,length(D1[upper.tri(D1, diag = T)]))
    md2<-dim(D1)[1]  }
  kont <- 0
  for(i1 in 1:md2){
    for(i2 in 1:(md2+1-i1)){
      kont <- kont+1
      di <- dD_alp[[i1]][[i2]]
      dalpha[kont] <- (-0.5)*sum(diag(solve(D1)%*%di-solve(D1)%*%di%*%solve(D1)*bb))     
    }
  }
  derivadas=cbind(t(dbeta),t(dsig2),t(dalpha))
  if(Arp!=0){
    gp = function(phi,sigma=sig2)ifelse(length(phi)==1,log(MatArp(phi,length(phi),sigma))
    ,log(det(MatArp(phi,length(phi),sigma))))
    dgp    = matrix(jacobian(gp,phi))
    dphi   = -1/sig2*(-Dphi(beta,y,x,z,b,p) + Dphiphi(beta,y,x,z,b,p)%*%phi)-1/2*dgp
    derivadas=cbind(t(dbeta),t(dsig2),t(dphi),t(dalpha))
    }
  return(derivadas)
}



# #x matriz com todas las covariables por columna
 nlf<-function(x,u,betas){
   resp<-log(exp(betas[1]+u[1])*exp(-(betas[2]+u[2])*x[,1])+
               exp(betas[3]+u[3])*exp(-(betas[4]+u[4]+betas[5]*x[,2])*x[,1]))/log(10)
   return(resp)
 }
# 
# exp<- expression(log(exp(bta1+b1)*exp(-(bta2+b2)*xx1)+
#                              exp(bta3+b3)*exp(-(bta4+b4+bta5*xx2)*xx1)/log(10)))
# nomes=c("xx1","xx2","xx3","xx4","xx5","xx6","xx7","xx8")
# nEA=2
# nEF=8
# names(x)=nomes
#  attach(x,warn.conflicts = FALSE)
#   
# betass=as.data.frame(t(beta1))
# names(betass)=paste("bta",1:nEF,sep ="")
# attach(betass,warn.conflicts = FALSE)
# 
# 
# b=(ub1[(((j-1)*q1)+1) : (j*q1), j])
# bs=as.data.frame(t(b))
# names(bs)=paste("b",1:nEA,sep ="")
# attach(bs,warn.conflicts = FALSE)
# 

# nlf<-function(exp,){
#   resp<-log(exp(betas[1]+u[1])*exp(-(betas[2]+u[2])*x[,1])+
#               exp(betas[3]+u[3])*exp(-(betas[4]+u[4]+betas[5]*x[,2])*x[,1]))/log(10)
#   return(resp)
# }


# for ( k in 1:dim(x)[2]) {
# x=as.data.frame(x)
#   }

#h=derivada de la funcion no lineal respecto de cada b
funcH<-function(x,cd4,u,betas){
  ni<-length(x)
  resp<-matrix(0,ni,4)
  resp[,1]<- 1/(exp(betas[1]+u[1]-x*betas[2]-x*u[2])+exp(betas[3]+u[3]-x*betas[4]-x*u[4]-x*betas[5]*cd4))/(log(2)+log(5))*exp(betas[1]+u[1]-x*betas[2]-x*u[2])
  resp[,2]<- -x/(exp(betas[1]+u[1]-x*betas[2]-x*u[2])+exp(betas[3]+u[3]-x*betas[4]-x*u[4]-x*betas[5]*cd4))/(log(2)+log(5))*exp(betas[1]+u[1]-x*betas[2]-x*u[2])
  resp[,3]<- 1/(exp(betas[1]+u[1]-x*betas[2]-x*u[2])+exp(betas[3]+u[3]-x*betas[4]-x*u[4]-x*betas[5]*cd4))/(log(2)+log(5))*exp(betas[3]+u[3]-x*betas[4]-x*u[4]-x*betas[5]*cd4)
  resp[,4]<- -x/(exp(betas[1]+u[1]-x*betas[2]-x*u[2])+exp(betas[3]+u[3]-x*betas[4]-x*u[4]-x*betas[5]*cd4))/(log(2)+log(5))*exp(betas[3]+u[3]-x*betas[4]-x*u[4]-x*betas[5]*cd4)
  return(resp)
}

#w=derivada de la funcion no lineal respecto de cada beta
funcW<-function(x,cd4,u,betas){
  ni<-length(x)
  resp<-matrix(0,ni,length(betas))
  
  resp[,1]=1/(exp(betas[1]+u[1]-x*betas[2]-x*u[2])+exp(betas[3]+u[3]-x*betas[4]-x*u[4]-x*betas[5]*cd4))/(log(2)+log(5))*exp(betas[1]+u[1]-x*betas[2]-x*u[2])
  resp[,2]=-x/(exp(betas[1]+u[1]-x*betas[2]-x*u[2])+exp(betas[3]+u[3]-x*betas[4]-x*u[4]-x*betas[5]*cd4))/(log(2)+log(5))*exp(betas[1]+u[1]-x*betas[2]-x*u[2])
  resp[,3]=1/(exp(betas[1]+u[1]-x*betas[2]-x*u[2])+exp(betas[3]+u[3]-x*betas[4]-x*u[4]-x*betas[5]*cd4))/(log(2)+log(5))*exp(betas[3]+u[3]-x*betas[4]-x*u[4]-x*betas[5]*cd4)
  resp[,4]=-x/(exp(betas[1]+u[1]-x*betas[2]-x*u[2])+exp(betas[3]+u[3]-x*betas[4]-x*u[4]-x*betas[5]*cd4))/(log(2)+log(5))*exp(betas[3]+u[3]-x*betas[4]-x*u[4]-x*betas[5]*cd4)
  resp[,5]=-cd4*x/(exp(betas[1]+u[1]-x*betas[2]-x*u[2])+exp(betas[3]+u[3]-x*betas[4]-x*u[4]-x*betas[5]*cd4))/(log(2)+log(5))*exp(betas[3]+u[3]-x*betas[4]-x*u[4]-x*betas[5]*cd4)
  return(resp)
}


# W=funcW(x1,(ub1[(((j-1)*q1)+1) : (j*q1), j]),beta1) #X
# H=funcH(x1,(ub1[(((j-1)*q1)+1) : (j*q1), j]),beta1) #Z
# yt=y1-nlf(x1,beta1,(ub1[(((j-1)*q1)+1) : (j*q1), j]))+W%*%beta1+H%*%(ub1[(((j-1)*q1)+1) : (j*q1), j])


