#setwd("C:/Users/rommy/OneDrive/Escritorio/Rommy/Pacote/Funciones")
#source('FunGeral.R')

# install.packages("mvtnorm")
# install.packages("mnormt")
# install.packages("lmec")
# install.packages("numDeriv")
# install.packages("xtable")
#install.packages("tmvtnorm")

#library(mvtnorm)
#library(mnormt)
#library(lmec)
#library(numDeriv)
#library(xtable)
#library(tmvtnorm)
#library(tcltk)
#library(MASS)

EMCensArpL=function(cc,y,x,z,nj,Arp, betai,sigmaei,D1i,pisi,cens.type="left", LI,LS,MaxIter,ee, Prev,step,isubj,xpre,zpre){
  start.time <- Sys.time()
  pb = tkProgressBar(title = "AR(p)MMEC by EM", min = 0,max = MaxIter, width = 300)
  setTkProgressBar(pb, 0, label=paste("Iter ",0,"/",MaxIter,"     -     ",0,"% done",sep = ""))
 
   ## the program admit left, righ and intervalar censored
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  if(cens.type=="left"){
    LI=rep(-Inf,length(cc))
    LS=rep(Inf,length(cc))
    LS[cc==1]=y[cc==1]
    LI=as.vector(LI)
    LS=as.vector(LS)
  }
  
  if(cens.type=="right"){
    LI=rep(-Inf,length(cc))
    LI[cc==1]=y[cc==1]
    LS=rep(Inf,length(cc))
    LI=as.vector(LI)
    LS=as.vector(LS)
  }
  
  if(cens.type=="interval"){
    LI=LI
    LS=LS
    LI=as.vector(LI)
    LS=as.vector(LS)
  }
  
  
  m<-length(nj)[1]
  N<-sum(nj)
  p<-dim(x)[2]
  q1<-dim(z)[2]
  m1<-m*p
  m2<-m*q1
  
  
  
  # #valores iniciais
  if(is.null(betai)){ betai=solve(t(x)%*%x)%*%t(x)%*%y  }
  if(is.null(sigmaei)){sigmaei= 0.25 }
  if(is.null(D1i)){ D1i=0.1*diag(dim(z)[2]) }
  beta1<-betai
  sigmae<- sigmaei
  D1<-D1i
  if(is.null(pisi)){ pisi = as.numeric(pacf((y - x%*%beta1),lag.max=Arp,plot=F)$acf) }
  
  pis<-pisi
  iD1<- solve(D1)
  #phi=phii
  
  if(Arp!="UNC"){
    phi = estphit(pis)}
  if(Arp=="UNC"){
    Arp=0
    pis = 0
    phi = 0}
  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],phi)
  teta1<- c(beta1,sigmae,phi)

 
  criterio<-1
  count<-0
  
  while(criterio > ee){
    
    count <- count + 1
   # print(count)
    soma1<- matrix(0,q1,q1)
    soma2<-0
    soma3<- matrix(0,p,p)
    soma4<- matrix(0,p,1)
    soma5<- matrix(0,p,p)
    MI <- matrix(0,p+1+length(D1[upper.tri(D1, diag = T)])+Arp,
                 p+1+length(D1[upper.tri(D1, diag = T)])+Arp) 
    
    ubi=matrix(0,m2,m)
    ubbi=matrix(0,m2,m2)
    uybi=matrix(0,N,m2)
    uyyi=matrix(0,N,N)
    uyi=matrix(0,N,m)
    xi=matrix(0,N,m1)
    zi=matrix(0,N,m2) 
    
    ver<-matrix(0,m,1)
    
    for (j in 1:m ){
      
      cc1=cc[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      
      y1=y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      
      #if(typeModel=="L"){
        
        x1=matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
        z1=matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
        
      #}
      
      # if(typeModel=="NL"){
      #   
      #   # W=funW[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      #   # H=funH[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      #   # yt=y1-funlf[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      #   # 
      #   # 
      #   # 
      #   # 
      #   # x1<-W
      #   # z1<-H
      #   # y1<-yt
      # }
      # 
      
      LI1<- LI[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      LS1<- LS[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      gammai=x1%*%beta1
      
      
      if(Arp==0){eGama=diag(1,nj[j])
      Gama=eGama*sigmae}
      if(Arp!=0){
        Gama<- MatArp(pis,nj[j],sigmae) 
        eGama<-Gama/sigmae
      }
      
      Psi<-(Gama+(z1)%*%D1%*%t(z1))
      Psi<-(Psi+t(Psi))/2
      delta<- solve(iD1+(t(z1)%*%solve(eGama)%*%(z1*(1/sigmae))))
      
      if(sum(cc1)==0){
        
        uy<- matrix(y1,nj[j],1)
        uyy<- y1%*%t(y1)
        ub<- delta%*%(t(z1)*(1/sigmae))%*%solve(eGama)%*%(uy-gammai)
        ubb<- delta+(delta%*%(t(z1)*((1/sigmae)^2))%*%solve(eGama)%*%(uyy-uy%*%t(gammai)-gammai%*%t(uy)+gammai%*%t(gammai))%*%solve(eGama)%*%z1%*%delta)
        uyb<- (uyy-uy%*%t(gammai))%*%solve(eGama)%*%(z1*(1/sigmae))%*%delta
        ver[j,]<- dmvnorm(as.vector(y1),gammai,Psi)
        
      }
      
      if(sum(cc1)>=1){
        
        
        if(sum(cc1)==nj[j]){
          muc=x1%*%beta1
          Sc<-Psi
          #aux<- MomemNT(muc,Sc,y1)
          #uy<-aux$Ey
          #uyy<- aux$Eyy
          aux<- mtmvnorm(mean=c(muc), sigma=Sc, lower=LI1, upper=LS1)
          uy<- aux$tmean
          uyy<- aux$tvar+uy%*%t(uy)
          
          ub<- delta%*%(t(z1)*(1/sigmae))%*%solve(eGama)%*%(uy-gammai)
          ubb<- delta+(delta%*%(t(z1)*((1/sigmae)^2))%*%solve(eGama)%*%(uyy-uy%*%t(gammai)-gammai%*%t(uy)+gammai%*%t(gammai))%*%solve(eGama)%*%z1%*%delta)
          uyb<- (uyy-uy%*%t(gammai))%*%solve(eGama)%*%(z1*(1/sigmae))%*%delta
          ver[j,]<- pmvnorm(LI1, LS1, mean=c(muc),sigma=Sc,algorithm = GB)#pmnorm(y1,c(muc),Sc)
          
        }
        
        else {
          
          muc=x1[cc1==1,]%*%beta1+Psi[cc1==1,cc1==0]%*%solve(Psi[cc1==0,cc1==0])%*%(y1[cc1==0]-x1[cc1==0,]%*%beta1)
          Sc <-Psi[cc1==1,cc1==1]-Psi[cc1==1,cc1==0]%*%solve(Psi[cc1==0,cc1==0])%*%Psi[cc1==0,cc1==1]
          #aux <-MomemNT(muc,Sc,y1[cc1==1])
          aux<- mtmvnorm(mean=c(muc), sigma=Sc, lower=LI1[cc1==1], upper=LS1[cc1==1])
          uy <-matrix(y1,nj[j],1)
          #uy[cc1==1]<-aux$Ey
          uy[cc1==1]<- aux$tmean
          uyy<-matrix(0,nj[j],nj[j])
          #uyy[cc1==1,cc1==1]<-aux$Vary
          uyy[cc1==1,cc1==1]<- aux$tvar
          uyy<- uyy+uy%*%t(uy)
          ub<- delta%*%(t(z1)*(1/sigmae))%*%solve(eGama)%*%(uy-gammai)
          ubb<- delta+(delta%*%(t(z1)*((1/sigmae)^2))%*%solve(eGama)%*%(uyy-uy%*%t(gammai)-gammai%*%t(uy)+gammai%*%t(gammai))%*%solve(eGama)%*%z1%*%delta)
          uyb<- (uyy-uy%*%t(gammai))%*%solve(eGama)%*%(z1*(1/sigmae))%*%delta
          ver[j,]<- dmvnorm(c(y1[cc1==0]),mean=c(gammai[cc1==0]),sigma=as.matrix(Psi[cc1==0,cc1==0]))*(pmvnorm(LI1[cc1==1],LS1[cc1==1],mean=as.vector(muc),sigma=Sc,algorithm = GB))[1]
          
          
          #dmnorm(y1[cc1==0],gammai[cc1==0],Psi[cc1==0,cc1==0])*pmnorm(y1[cc1==1],c(muc),Sc)
          
        }
        
      }
      
      soma1<- soma1 + ubb
      soma2<- soma2 + (sum(diag(uyy%*%solve(eGama)))-t(uy)%*%solve(eGama)%*%gammai-t(gammai)%*%solve(eGama)%*%uy-sum(diag(solve(eGama)%*%((uyb)%*%t(z1))))-sum(diag(solve(eGama)%*%((uyb)%*%t(z1))))
                       +t(gammai)%*%solve(eGama)%*%z1%*%ub+t(ub)%*%t(z1)%*%solve(eGama)%*%gammai+t(gammai)%*%solve(eGama)%*%gammai+sum(diag(ubb%*%t(z1)%*%solve(eGama)%*%z1)))
      soma3<- soma3 + (t(x1)%*%solve(eGama)%*%x1)
      soma4<- soma4 + (t(x1)%*%solve(eGama)%*%(uy-z1%*%ub))
      soma5<- soma5 + (t(x1)%*%solve(Psi)%*%x1-t(x1)%*%solve(Psi)%*%(uyy-uy%*%t(uy))%*%solve(Psi)%*%x1)
      
      
      ubi[(((j-1)*q1)+1) : (j*q1), j]<-ub
      ubbi[(((j-1)*q1)+1) : (j*q1), (((j-1)*q1)+1) : (j*q1)]<-ubb
      uybi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(((j-1)*q1)+1) : (j*q1)]<-uyb
      uyyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]<-uyy
      uyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),j]<-uy
      zi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])), (((j-1)*q1)+1) : (j*q1)]<-z1
      xi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])), (((j-1)*p)+1) : (j*p)]<-x1
      
      tetaMI=c(beta1,sigmae,phi)
      si<-Jt(tetaMI,uy,x1,z1,ub,ubb,p,Arp,D1)
      MI <- MI + t(si)%*%si
      
    }
    
    beta1<- solve(soma3)%*%soma4
    sigmae<- (1/N)*(soma2)
    sigmae<-as.numeric(sigmae)
    D1<- (1/m)*(soma1)
    iD1<- solve(D1) 
    teta1 <- c(beta1,sigmae)
    if(Arp!=0){
      pis <- optim(pis, method = "L-BFGS-B", FCiArp, lower =rep(-.999,Arp), upper =rep(.999,Arp), beta1=beta1,sigmae=sigmae, ubi=ubi,ubbi=ubbi,uybi=uybi,uyyi=uyyi,uyi=uyi,x=x,z=z,nj=nj,hessian=TRUE)$par
      phi=estphit(pis) 
      teta1 <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)], phi)  }
    
    #print(teta1)
    logver <- sum(log(ver))
    #print(logver)
    varbeta<-solve(soma5)
    
    ## if (count>1){
    ##   criterio <- (abs(logver1-logver))
    ## }
    
    if (count>=1){
      criterio <- sqrt((teta1/teta-1)%*%(teta1/teta-1))
      setTkProgressBar(pb, count, label=paste("Iter ",count,"/",MaxIter,"     -     ",floor((count)/(MaxIter)*100),"% done",sep = ""))
    
      }
    
    if (count==MaxIter){
      criterio <- 0.0000000000001
    }
    
    teta<-teta1
    logver1<-logver
    
  }
  
  dd<-D1[upper.tri(D1, diag = T)]
  
  npar<-length(c(teta1))
  
  ni<-sum(nj)
  
  loglik<-logver1
  
  AICc<- -2*loglik +2*npar
  AICcorr<- AICc + ((2*npar*(npar+1))/(ni-npar-1))
  BICc <- -2*loglik +log(ni)*npar
  
  if(Prev)
  { contt=0
    Predicao<- matrix(0,length(isubj),1+step)
    for (j in isubj ){
      #
      contt=contt+1
      IndPred=c(rep(0,nj[j]),rep(1,step))
      xobs=x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ]
      xprei=xpre[(step*contt-(step-1)) : (step*contt),  ]
      xobspre=rbind(xobs,xprei)
      
      zobs=z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ]
      zprei=zpre[(step*contt-(step-1)) : (step*contt),  ]
      zobspre=rbind(zobs,zprei)
      
      gammai = xobs%*%beta1
      yobs=uyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),j]
      
      if(Arp==0){ Gama=diag(1,nj[j]+step)*sigmae}
      if(Arp!=0){Gama<- MatArp(pis,nj[j]+step,sigmae)}
      PsiPred<-(Gama+(zobspre)%*%t(D1)%*%t(zobspre))
      Aux1Pred <- xprei%*%beta1
      Aux2Pred <- PsiPred[IndPred==1,IndPred==0]%*%solve(PsiPred[IndPred==0,IndPred==0])
      Aux3Pred <- (yobs-gammai)
      
      # Predicao[((step*contt-(step-1)) : (step*j)),1] <- j               
      #  Predicao[((step*contt-(step-1)) : (step*j)),2] <- Aux1Pred + Aux2Pred%*%Aux3Pred               
       Predicao[contt,1] <- j               
        Predicao[contt,2:(step+1)] <- Aux1Pred + Aux2Pred%*%Aux3Pred               
      
    }
    Predicao=as.data.frame(Predicao)
    colnames(Predicao) = c("subj",paste("step",1:step))
  }
  
  if(!Prev){Predicao=NULL}

  SE=round(sqrt(diag(ginv(MI))),3)
  intPar=round(1.96*SE,3)

  tableB  = data.frame(round(beta1,3),SE[1:p],paste("<",round(beta1,3)-round(intPar[1:p],3),",",round(beta1,3)+round(intPar[1:p],3),">"))
  rownames(tableB) = paste("beta",1:p)
  colnames(tableB) = c("Est","SE","IConf(95%)")
  
  
  if((round(sigmae,3)-round(intPar[p+1],3))<0) tableS  = data.frame(round(sigmae,3),SE[p+1],paste("<",0,",",round(sigmae,3)+round(intPar[p+1],3),">"))
  if((round(sigmae,3)-round(intPar[p+1],3))>=0) tableS  = data.frame(round(sigmae,3),SE[p+1],paste("<",round(sigmae,3)-round(intPar[p+1],3),",",round(sigmae,3)+round(intPar[p+1],3),">"))
  rownames(tableS) = "Sigma^2"
  colnames(tableS) = c("Est","SE","IConf(95%)")
  
  
  tableP  = data.frame(round(phi,3),SE[(p+2):(p+1+Arp)],paste("<",round(phi,3)-round(intPar[(p+2):(p+1+Arp)],3),",",round(phi,3)+round(intPar[(p+2):(p+1+Arp)],3),">"))
  rownames(tableP) = paste("Phi",1:Arp)
  colnames(tableP) = c("Est","SE","IConf(95%)")
  
  nnp=0
  for(al in 1:dim(D1)[1]) 
  {noa=paste(1:al,al,sep = "")
  nnp=c(nnp,noa)
  }
  nnp=nnp[-1]
  ici=round(dd,3)-round(intPar[(p+2+Arp):(p+1+length(D1[upper.tri(D1, diag = T)])+Arp)],3)
  ics=round(dd,3)+round(intPar[(p+2+Arp):(p+1+length(D1[upper.tri(D1, diag = T)])+Arp)],3)
  ici[as.numeric(nnp)%%11==0&ici<0]=0
  tableA  = data.frame(round(dd,3),SE[(p+2+Arp):(p+1+length(D1[upper.tri(D1, diag = T)])+Arp)],paste("<",ici,",",ics,">"))
  rownames(tableA) = paste("Alpha",nnp)
  colnames(tableA) = c("Est","SE","IConf(95%)")
  
  
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  obj.out <- list(beta1 = beta1, sigmae= sigmae, phi=phi, dd = dd, loglik=loglik,
                  AIC=AICc, BIC=BICc, AICcorr=AICcorr, iter = count, varbeta=varbeta,
                  ubi = ubi, ubbi = ubbi, uybi = uybi, uyi = uyi, uyyi = uyyi , MI=MI,
                  Prev= Predicao, time=time.taken, SE=SE,tableB=tableB,tableS=tableS,tableP=tableP,
                  tableA=tableA)
  
  
  if  (count == MaxIter)
  {
    setTkProgressBar(pb, MaxIter, label=paste("MaxIter reached ",count,"/",MaxIter,"    -    100 % done",sep = ""))
    close(pb)
  }
  else
  {
    setTkProgressBar(pb, MaxIter, label=paste("Convergence at Iter ",count,"/",MaxIter,"    -    100 % done",sep = ""))
    close(pb)
  }
  
  
  class(obj.out) <- "EM_NCens"
  
  return(obj.out)
  
}



  