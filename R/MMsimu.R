MMsimu=function(m,x,z,nj,beta,sigmae,D1,phi,percCensu,cens.type){
 Arp=length(phi)
 rx="NO";rz="NO"
  if(is.null(x)){x<-matrix(runif(sum(nj)*length(beta),-1,1),sum(nj),length(beta));rx="yes"}
  if(is.null(z)){z<-matrix(runif(sum(nj)*dim(D1)[1],-1,1),sum(nj),dim(D1)[1]);rz="yes"}
  #m=length(nj)
  y<-matrix(0,sum(nj),1)
  #phi = estphit(pis)
  for(i in 1:m){
    b<-rmvnorm(1,rep(0,dim(D1)[1]), D1)
    errorp <- c(arima.sim(n = nj[i], list(order=c(Arp,0,0), ar = phi), sd = sqrt(sigmae)))
    y[(sum(nj[1:i-1])+1) : (sum(nj[1:i])),]=x[(sum(nj[1:i-1])+1) : (sum(nj[1:i])),]%*%beta+z[(sum(nj[1:i-1])+1) : (sum(nj[1:i])),]%*%t(b)+errorp
  }
    yy=y
    y_cc=y
    cc=rep(0,length(y_cc))
   if(percCensu!=0)
   {
    if(cens.type=="left") {
      aa=sort(y, decreasing = FALSE)
      bb=aa[1:(percCensu*sum(nj))]
      cutof<-bb[percCensu*sum(nj)]
      cc=matrix(1,sum(nj),1)*(y< cutof)
      y[cc==1]=cutof
      y_cc=y}
     
    if(cens.type=="right") {
       aa=sort(y, decreasing = TRUE)
       bb=aa[1:(percCensu*sum(nj))]
       cutof<-bb[percCensu*sum(nj)]
       cc=matrix(1,sum(nj),1)*(y> cutof)
       y[cc==1]=cutof
       y_cc=y}
     
    if(cens.type=="interval") {
      aa=sort(y, decreasing = FALSE)
      bbi=aa[1:(percCensu*sum(nj)*0.5)]
      aa=sort(y, decreasing = TRUE)
      bbs=aa[1:(percCensu*sum(nj)*0.5)]
      cutofi<-bbi[percCensu*sum(nj)*0.5]
      cutofs<-bbs[percCensu*sum(nj)*0.5]
      cci=matrix(1,sum(nj),1)*(y< cutofi)
      y[cci==1]=cutofi
      ccs=matrix(1,sum(nj),1)*(y>cutofs)
      y[ccs==1]=cutofs
      y_cc=y
      cc=cci+ccs
     }
   }
  if(rx=="yes"){ 
    if(rz=="yes"){
      return(list(cc=cc, y_cc=y_cc, x=x,z=z))
    }else{return(list(cc=cc, y_cc=y_cc, x=x)) }
    }else{return(list(cc=cc, y_cc=y_cc))}
 # return(list(cc=cc, y_cc=y_cc, yy=yy))
}


