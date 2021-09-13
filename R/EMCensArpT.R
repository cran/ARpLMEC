
EMCensDECT<- function(cc,y,x,z,ttc,nj,struc,initial,cens.type,LL,LU,nu.fixed,iter.max,precision)
{
  start.time <- Sys.time()
  pb = tkProgressBar(title = "DEC-tLMEC by EM", min = 0,max = iter.max, width = 300)
  setTkProgressBar(pb, 0, label=paste("Iter ",0,"/",iter.max,"     -     ",0,"% done",sep = ""))
  
  
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  if(cens.type=="left"){
    LL=rep(-Inf,length(cc))
    LU=rep(Inf,length(cc))
    LU[cc==1]=y[cc==1]
    LL=as.vector(LL)
    LU=as.vector(LU)
  }
  
  if(cens.type=="right"){
    LL=rep(-Inf,length(cc))
    LL[cc==1]=y[cc==1]
    LU=rep(Inf,length(cc))
    LL=as.vector(LL)
    LU=as.vector(LU)
  }
  
  if(cens.type=="interval"){
    LL=LL
    LU=LU
    LL=as.vector(LL)
    LU=as.vector(LU)
  }
  
  m <- length(nj)[1]
  N <- sum(nj)
  q1 <- dim(z)[2]
  m2 <- m*q1
  p <- dim(x)[2]
  
  if(!is.null(initial)){
    beta1 <- matrix(initial$betas,p,1)
    sigmae <- initial$sigma2
    D1 <- initial$alphas
    iD1 <- solve(D1)
    iD1 <- (iD1 + t(iD1))/2
    nu <- initial$nu
    phis<-initial$phi
    
    if(struc=="DEC")
    {
      phi1 <- initial$phi1
      phi2 <- initial$phi2
      
      teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],phi1,phi2,nu)
      
    } else if(struc=="DEC(AR)"){
      phi1 <- initial$phi1
      phi2 <- 1
      teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],phi1,nu)
      
    } else if(struc=="SYM"){
      phi1 <- initial$phi1
      phi2 <- 0
      teta <- c(beta1, sigmae,D1[upper.tri(D1, diag = T)],phi1,nu)
      
    } else {
      phi1 <- NULL
      phi2 <- NULL
      teta <- c(beta1, sigmae,D1[upper.tri(D1, diag = T)],nu)
    }
    
  }
  
  if(is.null(initial)){ 
    beta1=solve(t(x)%*%x)%*%t(x)%*%y  
    sigmae= 0.25 
    D1=0.1*diag(dim(z)[2])
    nu=3
     iD1<- solve(D1)
    if(struc=="DEC")
    {
      phi1 <- 0.1
      phi2 <- 0.1
      
      teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],phi1,phi2,nu)
      
    } else if(struc=="DEC(AR)"){
      phi1 <- 0.1
      phi2 <- 1
      teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],phi1,nu)
      
    } else if(struc=="SYM"){
      phi1 <- 0.1
      phi2 <- 0
      teta <- c(beta1, sigmae,D1[upper.tri(D1, diag = T)],phi1,nu)
      
    } else {
      phi1 <- NULL
      phi2 <- NULL
      teta <- c(beta1, sigmae,D1[upper.tri(D1, diag = T)],nu)
    }
  }
  
  qr <- length(D1[lower.tri(D1, diag = T)])
  W <- x
  gamma1 <- as.vector(c(beta1))
  
  criterio <- 1
  count <- 0

  
  loglik <- logliktslmec(nu=nu,y=y,x=x,z=z,cc=cc,ttc=ttc,nj=nj,LL=LL,LU=LU,betas=beta1,sigmae=sigmae,D1=D1,phi1=phi1,phi2=phi2,struc=struc)
  
  loglikp <- loglik 
  
  
  while(criterio > precision){
    
    count <- count + 1
  
    
    soma1 <- matrix(0,q1,q1)
    soma2 <- 0
    soma3 <- matrix(0,p,p)
    soma4 <- matrix(0,p,1)
    
    soma7 <- matrix(0,p,p)
    
    Infbetas <- matrix(0,p,p)
    
    ui <- rep(0,m) 
    uyi <- matrix(0,N,m) 
    uyyi <- matrix(0,N,N) 
    ubi <- matrix(0,m2,m)  
    ubbi <- matrix(0,m2,m2) 
    uybi <- matrix(0,N,m2)  
    yest <- matrix(0,N,1)     
    biest <- matrix(0,m2,m) 
    yhi <- matrix(0,N,1)   
    
   
    for (j in 1:m){
   
      cc1 <- cc[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      y1 <- y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      x1 <- matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
      z1 <- matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
      tt1 <- ttc[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      W1 <- x1
      
      LL1 <- LL[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      LU1 <- LU[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      
      muii <- W1%*%gamma1
      Gama <- MatDec(tt1,phi1,phi2,struc)
      invGama <- solve(Gama)
      SIGMA <- (sigmae*Gama + (z1)%*%D1%*%t(z1)) 
      SIGMA <-(SIGMA+t(SIGMA))/2
      SIGMAinv <- solve(SIGMA)
      Lambda1 <- solve(iD1 + (t(z1)%*%invGama%*%z1)*(1/sigmae))
      Lambda1 <- (Lambda1 + t(Lambda1))/2
      
      dm <- as.numeric(t(y1 - muii)%*%SIGMAinv%*%(y1-muii))
      cdm <- as.numeric((nu+nj[j])/(nu+dm))
      
      if(sum(cc1)==0)
      {
        u <- cdm
        uy <- matrix(y1,nj[j],1)*cdm
        uyy <- (y1%*%t(y1))*cdm
        
        ub <- (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uy - u*muii)
        ubb <- Lambda1 + (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uyy - uy%*%t(muii) - muii%*%t(uy) + u*muii%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)
        uyb <- (uyy - uy%*%t(muii))%*%(invGama%*%(z1*(1/sigmae))%*%Lambda1)
        
        yh <- matrix(y1,nj[j],1)
        best <- (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(yh - muii)
    
        Eu2yy <- (cdm^2)*(y1%*%t(y1))
        Eu2y <- (cdm^2)*(matrix(y1,nj[j],1))
        Eu2 <- cdm^2 
        
        E2 <- Eu2yy - Eu2y%*%t(muii) - muii%*%t(Eu2y) + Eu2*(muii)%*%t(muii)
        E1 <- (uy - u*(muii))%*%t(uy - u*(muii))
        
       
      }
      if(sum(cc1)>=1)
      {
        if(sum(cc1)==nj[j])
        {
          
             aux1U <- MomTrunc::pmvnormt(lower = as.vector(LL1),upper = as.vector(LU1),mean = as.vector(muii),sigma = as.matrix((nu/(nu + 2))*SIGMA),nu = (nu+2))
          aux2U <- MomTrunc::pmvnormt(lower = as.vector(LL1),upper=as.vector(LU1), mean = as.vector(muii),sigma = as.matrix(SIGMA),nu=nu)
          
          u <- as.numeric(aux1U/aux2U)
          
          auxy <- MomTrunc::meanvarTMD(lower = as.vector(LL1),upper=as.vector(LU1),mu = as.vector(muii), Sigma = as.matrix((nu/(nu + 2))*SIGMA),nu=(nu+2),dist = "t")
          
          uy <- u*auxy$mean
          uyy <- u*auxy$EYY
          
          ub <- (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uy - u*muii)
          ubb <- Lambda1 + (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uyy - uy%*%t(muii) - muii%*%t(uy) + u*muii%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)
          uyb <- (uyy - uy%*%t(muii))%*%(invGama%*%(z1*(1/sigmae))%*%Lambda1)
          
          
          auxb <- MomTrunc::onlymeanTMD(lower = as.vector(LL1),upper=as.vector(LU1),mu = as.vector(muii), Sigma = as.matrix(SIGMA),nu=(nu),dist = "t")
          yh <- auxb
          best <- (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(yh - muii)
           
          cp <- (((nu+nj[j])/nu)^2)*((gamma((nu+nj[j])/2)*gamma((nu+4)/2))/(gamma(nu/2)*gamma((nu+nj[j]+4)/2)))
          
          auxw <- MomTrunc::meanvarTMD(lower = as.vector(LL1),upper=as.vector(LU1),mu = as.vector(muii), Sigma = as.matrix((nu/(nu + 4))*SIGMA),nu=(nu+4),dist = "t")
          
           auxEU <- MomTrunc::pmvnormt(lower = as.vector(LL1),upper=as.vector(LU1), mean = as.vector(muii),sigma = as.matrix((nu/(nu + 4))*SIGMA),nu=(nu+4))
          
          
          Eu2yy <- cp*(auxEU/aux2U)*auxw$EYY
          Eu2y <- cp*(auxEU/aux2U)*auxw$mean
          Eu2 <- cp*(auxEU/aux2U)
          
          E2 <- Eu2yy - Eu2y%*%t(muii) - muii%*%t(Eu2y) + Eu2*(muii)%*%t(muii)
          E1 <- (uy - u*(muii))%*%t(uy - u*(muii))
            }
        else{
          
          muiic <-  W1[cc1==1,]%*%gamma1 + SIGMA[cc1==1,cc1==0]%*%solve(SIGMA[cc1==0,cc1==0])%*%(y1[cc1==0]-W1[cc1==0,]%*%gamma1)
          Si <- SIGMA[cc1==1,cc1==1]-SIGMA[cc1==1,cc1==0]%*%solve(SIGMA[cc1==0,cc1==0])%*%SIGMA[cc1==0,cc1==1]
          Si <- (Si+t(Si))/2
          
          Qy0 <- as.numeric(t(y1[cc1==0]-W1[cc1==0,]%*%gamma1)%*%solve(SIGMA[cc1==0,cc1==0])%*%(y1[cc1==0]-W1[cc1==0,]%*%gamma1))
          
          auxQy0 <- as.numeric((nu + Qy0)/(nu + length(cc1[cc1==0])))
          auxQy02 <- as.numeric((nu + Qy0)/(nu + 2 + length(cc1[cc1==0])))
          
          Sc0 <- auxQy0*Si
          Sc0til <- auxQy02*Si
          
          LL1c <- LL1[cc1==1]
          LU1c <- LU1[cc1==1]
  
          aux1U <- MomTrunc::pmvnormt(lower = as.vector(LL1c),upper = as.vector(LU1c),mean = as.vector(muiic),sigma = as.matrix(Sc0til),nu = (nu + 2 + length(cc1[cc1==0])))
          aux2U <- MomTrunc::pmvnormt(lower = as.vector(LL1c),upper=as.vector(LU1c), mean = as.vector(muiic),sigma =as.matrix(Sc0),nu= (nu + length(cc1[cc1==0])))
          
          u <- as.numeric(aux1U/aux2U)*(1/auxQy0)
          Sc0til=round((Sc0til+t(Sc0til))/2,3)
          
          auxy <- MomTrunc::meanvarTMD(lower = as.vector(LL1c),upper=as.vector(LU1c),mu = as.vector(muiic), Sigma = as.matrix(Sc0til),nu=(nu + 2 + length(cc1[cc1==0])),dist = "t")
          w1aux <- auxy$mean
          w2aux <- auxy$EYY
          
          uy <- matrix(y1,nj[j],1)*u
          uy[cc1==1] <- w1aux*u
          
          uyy <- y1%*%t(y1)*u
          uyy[cc1==0,cc1==1] <- u*y1[cc1==0]%*%t(w1aux)
          uyy[cc1==1,cc1==0] <- u*w1aux%*%t(y1[cc1==0])
          uyy[cc1==1,cc1==1] <- u*w2aux
          uyy <- (uyy + t(uyy))/2
          
          ub <- (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uy - u*muii)
          ubb <- Lambda1 + (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uyy - uy%*%t(muii) - muii%*%t(uy) + u*muii%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)
          uyb <- (uyy - uy%*%t(muii))%*%(invGama%*%(z1*(1/sigmae))%*%Lambda1)
          
          auxb <- MomTrunc::onlymeanTMD(lower = as.vector(LL1c),upper=as.vector(LU1c),mu = as.vector(muiic), Sigma = as.matrix(Sc0),nu=(nu + length(cc1[cc1==0])),dist = "t")
          yh <- matrix(y1,nj[j],1)
          yh[cc1==1] <- auxb
          
          best <- (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(yh - muii)
          
        
          
          dp <- ((nu+nj[j])^2)*((gamma((nu+nj[j])/2)*gamma((nu+4+length(cc1[cc1==0]))/2))/(gamma((nu+length(cc1[cc1==0]))/2)*gamma((nu+4+nj[j])/2)))
          
          Sirc=(nu + Qy0)/(nu + 4 + length(cc1[cc1==0]))*Si
          Sirc=round(((t(Sirc)+Sirc)/2),2)
          auxEU <- MomTrunc::pmvnormt(lower = as.vector(LL1c),upper=as.vector(LU1c), mean = as.vector(muiic),sigma = as.matrix(Sirc),nu=(nu + 4 + length(cc1[cc1==0])))
          
          auxEw <- MomTrunc::meanvarTMD(lower = as.vector(LL1c),upper=as.vector(LU1c),mu = as.vector(muiic), Sigma = as.matrix(Sirc),nu=(nu + 4 + length(cc1[cc1==0])),dist = "t")
          
          Ew1aux <- auxEw$mean
          Ew2aux <- auxEw$EYY
          
          Eu2yy <- (dp/((nu + Qy0)^2))*(auxEU/aux2U)*y1%*%t(y1)
          Eu2yy[cc1==0,cc1==1] <- (dp/((nu + Qy0)^2))*(auxEU/aux2U)*y1[cc1==0]%*%t(Ew1aux)
          Eu2yy[cc1==1,cc1==0] <- (dp/((nu + Qy0)^2))*(auxEU/aux2U)* Ew1aux%*%t(y1[cc1==0])
          Eu2yy[cc1==1,cc1==1] <- (dp/((nu + Qy0)^2))*(auxEU/aux2U)*Ew2aux
          
          Eu2y <- (dp/((nu + Qy0)^2))*(auxEU/aux2U)*matrix(y1,nj[j],1)
          Eu2y[cc1==1] <- (dp/((nu + Qy0)^2))*(auxEU/aux2U)*Ew1aux
          
          Eu2 <- (dp/((nu + Qy0)^2))*(auxEU/aux2U)
          
          E2 <- Eu2yy - Eu2y%*%t(muii) - muii%*%t(Eu2y) + Eu2*(muii)%*%t(muii)
          E1 <- (uy - u*(muii))%*%t(uy - u*(muii))
          
               
        }    
        
      } 
      
      soma1 <- soma1 + ubb
      soma2 <- soma2 + (sum(diag(uyy%*%invGama)) - t(uy)%*%invGama%*%muii - t(muii)%*%invGama%*%uy - sum(diag(t(uyb)%*%invGama%*%z1)) - sum(diag(uyb%*%t(z1)%*%invGama))
                        + t(muii)%*%invGama%*%z1%*%ub  + t(ub)%*%t(z1)%*%invGama%*%muii + u*t(muii)%*%invGama%*%muii + sum(diag(ubb%*%t(z1)%*%invGama%*%z1))) 
      soma3 <- soma3 + (u*t(x1)%*%invGama%*%x1)
      soma4 <- soma4 + (t(x1)%*%invGama%*%(uy - z1%*%ub))
         
   
      soma7 <- soma7 +  (((nu+nj[j])/(nu+nj[j]+2))*t(x1)%*%SIGMAinv%*%x1 - ((nu+nj[j]+2)/(nu+nj[j]))*t(x1)%*%SIGMAinv%*%(E2)%*%SIGMAinv%*%x1 + (t(x1)%*%SIGMAinv%*%(E1)%*%SIGMAinv%*%x1))
        
      ui[j] <- u
      uyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),j] <- uy
      uyyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]<- uyy
      ubi[(((j-1)*q1)+1) : (j*q1), j] <- ub
      ubbi[(((j-1)*q1)+1) : (j*q1), (((j-1)*q1)+1) : (j*q1)]<- ubb
      uybi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(((j-1)*q1)+1) : (j*q1)]<- uyb
      yest[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))] <- z1%*%best + muii
      biest[(((j-1)*q1)+1) : (j*q1), j] <- best 
      yhi[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))] <- yh
      
    } 
    
    yhatorg <- apply(yhi,1,sum)
    yfit <- apply(yest,1,sum)  
    yfit[cc==1] <- yhatorg[cc==1]
    
   
    Infbetas[1:p,1:p] <- soma7
    
    Infbetas <- (Infbetas + t(Infbetas))/2         
    
    beta1 <- solve(soma3)%*%soma4
    gamma1 <- as.vector(c(beta1))
    sigmae <- (1/N)*(soma2)
    sigmae <- as.numeric(sigmae)
    D1 <- (1/m)*(soma1)
    iD1 <- solve(D1)
    
   
    if(nu.fixed==FALSE)
    {
      nu <- optimize(f = logliktslmec,interval = c(3,30),tol = 0.00001, maximum = TRUE,y=y,x=x,z=z,cc=cc,ttc=ttc,nj=nj,
                     LL=LL,LU=LU,betas=beta1,sigmae=sigmae,D1=D1,phi1=phi1,phi2=phi2,struc=struc)$maximum
    }
    
  
    if(struc=="DEC")                                                                     
    {
      phis <- optim(c(phi1,phi2), FCit,lower =c(0.01,0.01), upper=c(0.9,30),method = "L-BFGS-B", hessian=TRUE, beta1=beta1,sigmae=sigmae,ttc=ttc,ubi=ubi,ubbi=ubbi,uybi=uybi,uyyi=uyyi,uyi=uyi,ui=ui,x=x,z=z,nj=nj,struc=struc)$par
      phi1 <- phis[1]
      phi2 <- phis[2]
     teta1 <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],phi1,phi2,nu)
     
    } else if(struc=="DEC(AR)"){
      phi2 <- 1                                   
      phi1 <- optimize(f=FCiphi1t, lower= 0.0001, upper=0.9, tol = 0.001, phi2=phi2, beta1=beta1,sigmae=sigmae,
                       ttc=ttc,ubi=ubi,ubbi=ubbi,uybi=uybi,uyyi=uyyi,uyi=uyi,ui=ui,x=x,z=z,nj=nj,struc=struc)$minimum
     
      teta1 <- c(beta1, sigmae,D1[upper.tri(D1, diag = T)],phi1,nu)
    
    } else if(struc=="SYM"){
      phi2 <- 0
      phi1 <- optimize(f=FCiphi1t, lower= 0.0001, upper=0.9,tol = 0.001, phi2=phi2, beta1=beta1,sigmae=sigmae,
                       ttc=ttc,ubi=ubi,ubbi=ubbi,uybi=uybi,uyyi=uyyi,uyi=uyi,ui=ui,x=x,z=z,nj=nj,struc=struc)$minimum
      
      teta1 <- c(beta1, sigmae,D1[upper.tri(D1, diag = T)],phi1,nu)
      
    } else{
     
      teta1 <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],nu)
  
    }
    
    
    
    loglik <- logliktslmec(nu=nu,y=y,x=x,z=z,cc=cc,ttc=ttc,nj=nj,LL=LL,LU=LU,betas=beta1,sigmae=sigmae,D1=D1,phi1=phi1,phi2=phi2,struc=struc)
       loglikp1 <- loglik 
       
       
       
       if(count > 1){
         criterio <- sqrt(((loglikp1/loglikp)-1)%*%((loglikp1/loglikp)-1))
         setTkProgressBar(pb, count, label=paste("Iter ",count,"/",iter.max,"     -     ",floor((count)/(iter.max)*100),"% done",sep = ""))
       }    
       if(count==iter.max){criterio <- precision*0.0001}

     
       
    teta <- teta1
    loglikp <- loglikp1
    
    
  } 
  
  
  dd <- D1[upper.tri(D1, diag = T)]
  
  
  npar <- length(c(teta1))
  AICc <- -2*loglikp + 2*npar
  BICc <- -2*loglikp + log(N)*npar
  
  SE=round(sqrt(diag(solve(Infbetas))),3)
  intPar=round(qt(0.975,nu)*SE,3)
  
  
  tableB  = data.frame(round(beta1,3),SE,paste("<",round(beta1,3)-round(intPar,3),",",round(beta1,3)+round(intPar,3),">"))
  rownames(tableB) = paste("beta",1:p)
  colnames(tableB) = c("Est","SE","IConf(95%)")
  
  
  tableS  = data.frame(round(sigmae,3))
  rownames(tableS) = "Sigma^2"
  colnames(tableS) = c("Est")
  
  
  
  if(struc=="DEC"){
    phi=c(phi1,phi2)  
    tableP  = data.frame(round(phi,3))
    rownames(tableP) = paste("Phi",1:2)
    colnames(tableP) = c("Est")
  }
  if(struc=="DEC(AR)"){
    phi=phi1 
    tableP  = data.frame(round(phi,3))
    rownames(tableP) = paste("Phi",1)
    colnames(tableP) = c("Est")
  }
  if(struc=="SYM"){
    phi=phi1  
    tableP  = data.frame(round(phi,3))
    colnames(tableP) = c("Est")
  }
  if(struc=="UNC"){ phi=NULL; tableP =NULL }
  
   
  nnp=0
  for(al in 1:dim(D1)[1]) 
  {noa=paste(1:al,al,sep = "")
  nnp=c(nnp,noa)
  }
  nnp=nnp[-1]
  tableA  = data.frame(round(dd,3))
  rownames(tableA) = paste("Alpha",nnp)
  colnames(tableA) = c("Est")
  
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  
  obj.out <- list(beta1 = beta1, sigmae= sigmae, phi=phi, dd = dd,nu=nu, loglik=loglik,
                  AIC=AICc, BIC=BICc, iter = count,
                  ubi = ubi, ubbi = ubbi, uybi = uybi, uyi = uyi, uyyi = uyyi,ui=ui, MI=Infbetas,
                  Prev= NULL, time=time.taken, SE=SE,tableB=tableB,tableS=tableS,tableP=tableP,
                  tableA=tableA,yest=yfit, yog = yhatorg)
  
  if  (count == iter.max)
  {
    setTkProgressBar(pb, iter.max, label=paste("MaxIter reached ",count,"/",iter.max,"    -    100 % done",sep = ""))
    close(pb)
  }
  else
  {
    setTkProgressBar(pb, iter.max, label=paste("Convergence at Iter ",count,"/",iter.max,"    -    100 % done",sep = ""))
    close(pb)
  }
  
  class(obj.out) <- "DECTLMEC"
 
  return(obj.out)
}

EMCensArpT<- function(cc,y,x,z,ttc,nj,Arp,initial,cens.type,LL,LU,nu.fixed,iter.max,precision)
{
  
  start.time <- Sys.time()
  pb = tkProgressBar(title = "AR(p)-tLMEC by EM", min = 0,max = iter.max, width = 300)
  setTkProgressBar(pb, 0, label=paste("Iter ",0,"/",iter.max,"     -     ",0,"% done",sep = ""))
  
  
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  if(cens.type=="left"){
    LL=rep(-Inf,length(cc))
    LU=rep(Inf,length(cc))
    LU[cc==1]=y[cc==1]
    LL=as.vector(LL)
    LU=as.vector(LU)
  }
  
  if(cens.type=="right"){
    LL=rep(-Inf,length(cc))
    LL[cc==1]=y[cc==1]
    LU=rep(Inf,length(cc))
    LL=as.vector(LL)
    LU=as.vector(LU)
  }
  
  if(cens.type=="interval"){
    LL=LL
    LU=LU
    LL=as.vector(LL)
    LU=as.vector(LU)
  }
  

  m <- length(nj)[1]
  N <- sum(nj)
  q1 <- dim(z)[2]
  m2 <- m*q1
  p <- dim(x)[2]
  

  if(!is.null(initial)){
    beta1 <- matrix(initial$betas,p,1)
    sigmae <- initial$sigma2
    D1 <- initial$alphas
    iD1 <- solve(D1)
    iD1 <- (iD1 + t(iD1))/2
    nu <- initial$nu
    pii = as.numeric(pacf((y - x%*%beta1),lag.max=Arp,plot=F)$acf) 
    phis<-initial$phi
    
    if(Arp!="UNC"){
      pii = as.numeric(pacf((y - x%*%beta1),lag.max=Arp,plot=F)$acf)
      phi = phis}
    if(Arp=="UNC"){
      pii=0
      Arp=0
      pis = 0
      phi = 0}
  }
  
  if(is.null(initial)){ 
    beta1=solve(t(x)%*%x)%*%t(x)%*%y  
    sigmae= 0.25 
    D1=0.1*diag(dim(z)[2])
    nu=3
    pis = as.numeric(pacf((y - x%*%beta1),lag.max=Arp,plot=F)$acf)
    iD1<- solve(D1)
    if(Arp!="UNC"){
      pis = as.numeric(pacf((y - x%*%beta1),lag.max=Arp,plot=F)$acf)
      pii=pis
      phi = estphit(pis)}
    if(Arp=="UNC"){
      Arp=0
      pii=0
      pis = 0
      phi = 0}
  }

    qr <- length(D1[lower.tri(D1, diag = T)])
  W <- x
  gamma1 <- as.vector(c(beta1))
  

  
  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],pii,nu)


  criterio <- 1
  count <- 0
  
  loglik <- logliktArplmec(nu=nu,y=y,x=x,z=z,cc=cc,ttc=ttc,nj=nj,LL=LL,LU=LU,betas=beta1,sigmae=sigmae,D1=D1,pii)
  
  loglikp <- loglik 
  
  
  while(criterio > precision){
    
    count <- count + 1
    soma1 <- matrix(0,q1,q1)
    soma2 <- 0
    soma3 <- matrix(0,p,p)
    soma4 <- matrix(0,p,1)
     soma7 <- matrix(0,p,p)
    Infbetas <- matrix(0,p,p)
    
    ui <- rep(0,m) 
    uyi <- matrix(0,N,m) 
    uyyi <- matrix(0,N,N) 
    ubi <- matrix(0,m2,m)  
    ubbi <- matrix(0,m2,m2) 
    uybi <- matrix(0,N,m2)  
    yest <- matrix(0,N,1)    
    biest <- matrix(0,m2,m) 
    yhi <- matrix(0,N,1)    
   
    
    for (j in 1:m){
      cc1 <- cc[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      y1 <- y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      x1 <- matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
      z1 <- matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
      tt1 <- ttc[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      W1 <- x1
      
      LL1 <- LL[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      LU1 <- LU[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      
      
      
      muii <- W1%*%gamma1
      
      
      if(Arp==0){Gama=diag(1,nj[j])
      eGamma=Gama*sigmae}
      if(Arp!=0){
        eGamma<-MatArp(pii,tt1,sigmae)
        Gama <-eGamma/sigmae
      }
        
     
      
      invGama <- solve(Gama)
      SIGMA <- (sigmae*Gama + (z1)%*%D1%*%t(z1)) 
      SIGMA <-(SIGMA+t(SIGMA))/2
      SIGMAinv <- solve(SIGMA)
      Lambda1 <- solve(iD1 + (t(z1)%*%invGama%*%z1)*(1/sigmae))
      Lambda1 <- (Lambda1 + t(Lambda1))/2
      
      dm <- as.numeric(t(y1 - muii)%*%SIGMAinv%*%(y1-muii))
      cdm <- as.numeric((nu+nj[j])/(nu+dm))
      
      if(sum(cc1)==0)
      {
        u <- cdm
        uy <- matrix(y1,nj[j],1)*cdm
        uyy <- (y1%*%t(y1))*cdm
        
        ub <- (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uy - u*muii)
        ubb <- Lambda1 + (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uyy - uy%*%t(muii) - muii%*%t(uy) + u*muii%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)
        uyb <- (uyy - uy%*%t(muii))%*%(invGama%*%(z1*(1/sigmae))%*%Lambda1)
        
        yh <- matrix(y1,nj[j],1)
        best <- (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(yh - muii)
        
       
        
        Eu2yy <- (cdm^2)*(y1%*%t(y1))
        Eu2y <- (cdm^2)*(matrix(y1,nj[j],1))
        Eu2 <- cdm^2 
        
        E2 <- Eu2yy - Eu2y%*%t(muii) - muii%*%t(Eu2y) + Eu2*(muii)%*%t(muii)
        E1 <- (uy - u*(muii))%*%t(uy - u*(muii))
        
            }
      if(sum(cc1)>=1)
      {
        if(sum(cc1)==nj[j])
        {
          aux1U <- MomTrunc::pmvnormt(lower = as.vector(LL1),upper = as.vector(LU1),mean = as.vector(muii),sigma = as.matrix((nu/(nu + 2))*SIGMA),nu = (nu+2))
          aux2U <- MomTrunc::pmvnormt(lower = as.vector(LL1),upper=as.vector(LU1), mean = as.vector(muii),sigma = as.matrix(SIGMA),nu=nu)
          
          u <- as.numeric(aux1U/aux2U)
          
          auxy <- MomTrunc::meanvarTMD(lower = as.vector(LL1),upper=as.vector(LU1),mu = as.vector(muii), Sigma = as.matrix((nu/(nu + 2))*SIGMA),nu=(nu+2),dist = "t")
          
          uy <- u*auxy$mean
          uyy <- u*auxy$EYY
          
          ub <- (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uy - u*muii)
          ubb <- Lambda1 + (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uyy - uy%*%t(muii) - muii%*%t(uy) + u*muii%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)
          uyb <- (uyy - uy%*%t(muii))%*%(invGama%*%(z1*(1/sigmae))%*%Lambda1)
          
          
          auxb <- MomTrunc::onlymeanTMD(lower = as.vector(LL1),upper=as.vector(LU1),mu = as.vector(muii), Sigma = as.matrix(SIGMA),nu=(nu),dist = "t")
          yh <- auxb
          best <- (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(yh - muii)
          
      
          
          cp <- (((nu+nj[j])/nu)^2)*((gamma((nu+nj[j])/2)*gamma((nu+4)/2))/(gamma(nu/2)*gamma((nu+nj[j]+4)/2)))
          
          auxw <- MomTrunc::meanvarTMD(lower = as.vector(LL1),upper=as.vector(LU1),mu = as.vector(muii), Sigma = as.matrix((nu/(nu + 4))*SIGMA),nu=(nu+4),dist = "t")
          
         auxEU <- MomTrunc::pmvnormt(lower = as.vector(LL1),upper=as.vector(LU1), mean = as.vector(muii),sigma = as.matrix((nu/(nu + 4))*SIGMA),nu=(nu+4))
          
          
          Eu2yy <- cp*(auxEU/aux2U)*auxw$EYY
          Eu2y <- cp*(auxEU/aux2U)*auxw$mean
          Eu2 <- cp*(auxEU/aux2U)
          
          E2 <- Eu2yy - Eu2y%*%t(muii) - muii%*%t(Eu2y) + Eu2*(muii)%*%t(muii)
          E1 <- (uy - u*(muii))%*%t(uy - u*(muii))
          
            }
        else{
          
          muiic <-  W1[cc1==1,]%*%gamma1 + SIGMA[cc1==1,cc1==0]%*%solve(SIGMA[cc1==0,cc1==0])%*%(y1[cc1==0]-W1[cc1==0,]%*%gamma1)
          Si <- SIGMA[cc1==1,cc1==1]-SIGMA[cc1==1,cc1==0]%*%solve(SIGMA[cc1==0,cc1==0])%*%SIGMA[cc1==0,cc1==1]
          Si <- (Si+t(Si))/2
          
          Qy0 <- as.numeric(t(y1[cc1==0]-W1[cc1==0,]%*%gamma1)%*%solve(SIGMA[cc1==0,cc1==0])%*%(y1[cc1==0]-W1[cc1==0,]%*%gamma1))
          
          auxQy0 <- as.numeric((nu + Qy0)/(nu + length(cc1[cc1==0])))
          auxQy02 <- as.numeric((nu + Qy0)/(nu + 2 + length(cc1[cc1==0])))
          
          Sc0 <- auxQy0*Si
          Sc0til <- auxQy02*Si
          
          LL1c <- LL1[cc1==1]
          LU1c <- LU1[cc1==1]
         
          aux1U <- MomTrunc::pmvnormt(lower = as.vector(LL1c),upper = as.vector(LU1c),mean = as.vector(muiic),sigma = as.matrix(Sc0til),nu = (nu + 2 + length(cc1[cc1==0])))
          aux2U <- MomTrunc::pmvnormt(lower = as.vector(LL1c),upper=as.vector(LU1c), mean = as.vector(muiic),sigma =as.matrix(Sc0),nu= (nu + length(cc1[cc1==0])))
          
          u <- as.numeric(aux1U/aux2U)*(1/auxQy0)
          Sc0til=round((Sc0til+t(Sc0til))/2,3)
          auxy <- MomTrunc::meanvarTMD(lower = as.vector(LL1c),upper=as.vector(LU1c),mu = as.vector(muiic), Sigma = as.matrix(Sc0til),nu=(nu + 2 + length(cc1[cc1==0])),dist = "t")
          w1aux <- auxy$mean
          w2aux <- auxy$EYY
          
          uy <- matrix(y1,nj[j],1)*u
          uy[cc1==1] <- w1aux*u
          
          uyy <- y1%*%t(y1)*u
          uyy[cc1==0,cc1==1] <- u*y1[cc1==0]%*%t(w1aux)
          uyy[cc1==1,cc1==0] <- u*w1aux%*%t(y1[cc1==0])
          uyy[cc1==1,cc1==1] <- u*w2aux
          uyy <- (uyy + t(uyy))/2
          
          ub <- (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uy - u*muii)
          ubb <- Lambda1 + (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uyy - uy%*%t(muii) - muii%*%t(uy) + u*muii%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)
          uyb <- (uyy - uy%*%t(muii))%*%(invGama%*%(z1*(1/sigmae))%*%Lambda1)
          
          auxb <- MomTrunc::onlymeanTMD(lower = as.vector(LL1c),upper=as.vector(LU1c),mu = as.vector(muiic), Sigma = as.matrix(Sc0),nu=(nu + length(cc1[cc1==0])),dist = "t")
          yh <- matrix(y1,nj[j],1)
          yh[cc1==1] <- auxb
          
          best <- (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(yh - muii)
          
           dp <- ((nu+nj[j])^2)*((gamma((nu+nj[j])/2)*gamma((nu+4+length(cc1[cc1==0]))/2))/(gamma((nu+length(cc1[cc1==0]))/2)*gamma((nu+4+nj[j])/2)))
          
          Sirc=as.numeric((nu + Qy0)/(nu + 4 + length(cc1[cc1==0])))*Si
          Sirc=round(Sirc,3)
          auxEU <- MomTrunc::pmvnormt(lower = as.vector(LL1c),upper=as.vector(LU1c), mean = as.vector(muiic),sigma = as.matrix(Sirc),nu=(nu + 4 + length(cc1[cc1==0])))
          
          auxEw <- MomTrunc::meanvarTMD(lower = as.vector(LL1c),upper=as.vector(LU1c),mu = as.vector(muiic), Sigma = as.matrix(Sirc),nu=(nu + 4 + length(cc1[cc1==0])),dist = "t")
          
          Ew1aux <- auxEw$mean
          Ew2aux <- auxEw$EYY
          
          Eu2yy <- (dp/((nu + Qy0)^2))*(auxEU/aux2U)*y1%*%t(y1)
          Eu2yy[cc1==0,cc1==1] <- (dp/((nu + Qy0)^2))*(auxEU/aux2U)*y1[cc1==0]%*%t(Ew1aux)
          Eu2yy[cc1==1,cc1==0] <- (dp/((nu + Qy0)^2))*(auxEU/aux2U)* Ew1aux%*%t(y1[cc1==0])
          Eu2yy[cc1==1,cc1==1] <- (dp/((nu + Qy0)^2))*(auxEU/aux2U)*Ew2aux
          
          Eu2y <- (dp/((nu + Qy0)^2))*(auxEU/aux2U)*matrix(y1,nj[j],1)
          Eu2y[cc1==1] <- (dp/((nu + Qy0)^2))*(auxEU/aux2U)*Ew1aux
          
          Eu2 <- (dp/((nu + Qy0)^2))*(auxEU/aux2U)
          
          E2 <- Eu2yy - Eu2y%*%t(muii) - muii%*%t(Eu2y) + Eu2*(muii)%*%t(muii)
          E1 <- (uy - u*(muii))%*%t(uy - u*(muii))
           }       
        
      } 
      
      soma1 <- soma1 + ubb
      soma2 <- soma2 + (sum(diag(uyy%*%invGama)) - t(uy)%*%invGama%*%muii - t(muii)%*%invGama%*%uy - sum(diag(t(uyb)%*%invGama%*%z1)) - sum(diag(uyb%*%t(z1)%*%invGama))
                        + t(muii)%*%invGama%*%z1%*%ub  + t(ub)%*%t(z1)%*%invGama%*%muii + u*t(muii)%*%invGama%*%muii + sum(diag(ubb%*%t(z1)%*%invGama%*%z1))) 
      soma3 <- soma3 + (u*t(x1)%*%invGama%*%x1)
      soma4 <- soma4 + (t(x1)%*%invGama%*%(uy - z1%*%ub))
     soma7 <- soma7 +  (((nu+nj[j])/(nu+nj[j]+2))*t(x1)%*%SIGMAinv%*%x1 - ((nu+nj[j]+2)/(nu+nj[j]))*t(x1)%*%SIGMAinv%*%(E2)%*%SIGMAinv%*%x1 + (t(x1)%*%SIGMAinv%*%(E1)%*%SIGMAinv%*%x1))
      ui[j] <- u
      uyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),j] <- uy
      uyyi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]<- uyy
      ubi[(((j-1)*q1)+1) : (j*q1), j] <- ub
      ubbi[(((j-1)*q1)+1) : (j*q1), (((j-1)*q1)+1) : (j*q1)]<- ubb
      uybi[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),(((j-1)*q1)+1) : (j*q1)]<- uyb
      yest[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))] <- z1%*%best + muii
      biest[(((j-1)*q1)+1) : (j*q1), j] <- best 
      yhi[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))] <- yh
      
      
      
    } 
    
    yhatorg <- apply(yhi,1,sum) 
    yfit <- apply(yest,1,sum)   
    yfit[cc==1] <- yhatorg[cc==1]
    
      Infbetas[1:p,1:p] <- soma7
      
    Infbetas <- (Infbetas + t(Infbetas))/2         
    beta1 <- solve(soma3)%*%soma4
    gamma1 <- as.vector(c(beta1))
    sigmae <- (1/N)*(soma2)
    sigmae <- as.numeric(sigmae)
    D1 <- (1/m)*(soma1)
    iD1 <- solve(D1)
    
     
    if(nu.fixed==FALSE)
    {
      nu <- optimize(f = logliktArplmec, interval = c(3,50),tol = 0.00001, maximum = TRUE,y=y,x=x,z=z,cc=cc,ttc=ttc,nj=nj,LL=LL,LU=LU,betas=beta1,sigmae=sigmae,D1=D1,pii=pii)$maximum
     
    }
    
    
    if(Arp!=0){
      pii <- optim(pii, method = "L-BFGS-B", FCiArpt, lower =rep(-.999,Arp), upper =rep(.999,Arp), beta1=beta1,sigmae=sigmae,ttc=ttc,ubi=ubi,ubbi=ubbi,uybi=uybi,uyyi=uyyi,uyi=uyi,ui=ui,x=x,z=z,nj=nj,hessian=TRUE)$par
      phi=estphit(pii)
      }
    if(Arp==0){phi=0}
    
    
     
     
 
    
    loglik <- logliktArplmec(nu=nu,y=y,x=x,z=z,cc=cc,ttc=ttc,nj=nj,LL=LL,LU=LU,betas=beta1,sigmae=sigmae,D1=D1,pii)
      loglikp1 <- loglik
    
    
    if(count > 1){
      criterio <- sqrt(((loglikp1/loglikp)-1)%*%((loglikp1/loglikp)-1))
      setTkProgressBar(pb, count, label=paste("Iter ",count,"/",iter.max,"     -     ",floor((count)/(iter.max)*100),"% done",sep = ""))
          }    
      if(count==iter.max){criterio <- precision*0.0001}
    
    loglikp <- loglikp1

  } 
  
  dd <- D1[upper.tri(D1, diag = T)]
 
  npar <- length(c(teta))
  AICc <- -2*loglikp + 2*npar
  BICc <- -2*loglikp + log(N)*npar
  SE=round(sqrt(diag(solve(Infbetas))),3)
  intPar=round(qt(0.975,nu)*SE,3)
  
  
  tableB  = data.frame(round(beta1,3),SE,paste("<",round(beta1,3)-round(intPar,3),",",round(beta1,3)+round(intPar,3),">"))
  rownames(tableB) = paste("beta",1:p)
  colnames(tableB) = c("Est","SE","IConf(95%)")
  
  
  tableS  = data.frame(round(sigmae,3))
   rownames(tableS) = "Sigma^2"
  colnames(tableS) = c("Est")
  
  if(Arp!=0){
  tableP  = data.frame(round(phi,3))
  rownames(tableP) = paste("Phi",1:Arp)
  colnames(tableP) = c("Est")}
  
  if(Arp==0){
    tableP  = NULL;phi=NULL}
  
  
  
  nnp=0
  
  for(al in 1:dim(D1)[1]) 
  {noa=paste(1:al,al,sep = "")
  nnp=c(nnp,noa)
  }
  nnp=nnp[-1]
  tableA  = data.frame(round(dd,3))
  rownames(tableA) = paste("Alpha",nnp)
  colnames(tableA) = c("Est")
 
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
 
    
    obj.out <- list(beta1 = beta1, sigmae= sigmae, phi=phi, dd = dd,nu=nu, loglik=loglik,
                    AIC=AICc, BIC=BICc, iter = count,
                    ubi = ubi, ubbi = ubbi, uybi = uybi, uyi = uyi, uyyi = uyyi,ui=ui, MI=Infbetas,
                    Prev= NULL, time=time.taken, SE=SE,tableB=tableB,tableS=tableS,tableP=tableP,
                    tableA=tableA,yest=yfit, yog = yhatorg)
    
    if  (count == iter.max)
    {
      setTkProgressBar(pb, iter.max, label=paste("MaxIter reached ",count,"/",iter.max,"    -    100 % done",sep = ""))
      close(pb)
    }
    else
    {
      setTkProgressBar(pb, iter.max, label=paste("Convergence at Iter ",count,"/",iter.max,"    -    100 % done",sep = ""))
      close(pb)
    }
    
  class(obj.out) <- "ARpTLMEC"
  
  return(obj.out)
  
}