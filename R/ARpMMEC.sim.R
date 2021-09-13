#' @title Generating Censored Autoregressive Dataset with Mixed Effects, for normal distribution.
#' @import mvtnorm
#' @import mnormt
#' @import numDeriv
#' @import tmvtnorm
#' @import tcltk
#' @import stats
#' @description This function simulates a censored response variable with autoregressive errors of order \code{p}, with mixed effect and a established censoring rate. This function returns the censoring vector and censored response vector. 
#' @param m Number of individuals
#' @param x Design matrix of the fixed effects of order \code{n x s}, corresponding to vector of fixed effects.  
#' @param z Design matrix of the random effects of order\code{n x b}, corresponding to vector of random effects. 
#' @param tt Vector \code{1 x n} with the time the measurements were made, where \code{n} is the total number of measurements for all individuals.
#' @param nj Vector \code{1 x m} with the number of observations for each subject,  where \code{m} is the total number of individuals.
#' @param beta Vector of values fixed effects.
#' @param sigmae It's the value for sigma. 
#' @param D Covariance Matrix for the random effects.
#' @param phi Vector of length \code{Arp}, of values for autoregressive parameters. 
#' @param p.cens Censoring level for the process. Default is \code{0}
#' @param cens.type \code{left} for left censoring, \code{right} for right censoring and \code{interval} for intervalar censoring. Default is \code{left}
#' @return returns list:
#' \item{cc}{Vector of censoring indicators.}
#' \item{y_cc}{Vector of responses censoring.}
#' @examples
#' \dontrun{
#'  p.cens   = 0.1
#'  m           = 50
#'  D = matrix(c(0.049,0.001,0.001,0.002),2,2)
#'  sigma2 = 0.30
#'  phi    = c(0.48,-0.2)
#'  beta   = c(1,2,1)
#'  nj=rep(6,m) 
#'  tt=rep(seq(1:6),m)
#'  x<-matrix(runif(sum(nj)*length(beta),-1,1),sum(nj),length(beta))
#'  z<-matrix(runif(sum(nj)*dim(D)[1],-1,1),sum(nj),dim(D)[1])
#'  data=ARpMMEC.sim(m,x,z,tt,nj,beta,sigma2,D,phi,p.cens)
#'  y<-data$y_cc
#'  cc<-data$cc
#' }
#' @export
ARpMMEC.sim=function(m,x=NULL,z=NULL,tt=NULL,nj,beta,sigmae,D,phi,p.cens= 0,cens.type="left")
  {
  
   if(m==sum(nj))                        stop("not compatible sizes between m and nj")

  if(!is.null(x)){
    if(!is.numeric(x))                    stop("x must be a numeric matrix. Check documentation!")
    if(sum(is.na(x))>0)                   stop("There are some NA values in x.")
    if(!is.matrix(x))                     stop("x must be a matrix. Check documentation!")
    if(det(t(x)%*%x)==0)                  stop("the columns of x must be linearly independent.")
    if(length(x)==0)                      stop("The parameter x must be provided.")
    if(dim(x)[1]!=dim(z)[1])              stop("not compatible sizes between x and z")
    if(dim(x)[1]!=sum(nj))                stop("not compatible sizes between x and nj")
    if(dim(x)[2]!=length(beta))           stop("not compatible sizes between x and beta")
    if(dim(x)[1]!=length(tt))             stop("not compatible sizes between x and tt")
    
  }
  
    if(!is.null(z)){
    if(!is.numeric(z))                    stop("z must be a numeric matrix. Check documentation!")
    if(sum(is.na(z))>0)                   stop("There are some NA values in z.")
    if(!is.matrix(z))                     stop("z must be a matrix. Check documentation!")
    if(length(z)==0)                      stop("The parameter z must be provided.")
    if(dim(z)[1]!=sum(nj))                stop("not compatible sizes between z and nj")
    if(dim(z)[2]!=dim(D)[2])             stop("not compatible sizes between z and D")
  }
 
   if(!is.numeric(nj))                   stop("nj must be a numeric vector. Check documentation!")
  if(!is.vector(nj))                    stop("nj must be a vector. Check documentation!")
  if(sum(is.na(nj))>0)                  stop("There are some NA values in nj")
  if(length(nj)==0)                     stop("The parameter nj must be provided.")
  
    if(!is.numeric(beta))             stop("beta must be a numeric vector. Check documentation!")
    if(!is.vector(beta))              stop("beta must be a vector. Check documentation!")
    if(!is.numeric(sigmae))           stop("sigmae must be a scalar.")
    if(length(sigmae)>1)              stop("beta must be a scalar.")
    if(!is.matrix(D))                stop("D must be a matrix.")
    if(D[upper.tri(D)]!=D[lower.tri(D)])stop("D must be a simetric matrix.")
    if(!is.numeric(phi))              stop("phi must be a numeric vector. Check documentation!")
    if(sum(abs(phi))>=1)                   stop("the sum of the phi must be less than 1. Check documentation!")
  
 
    if(cens.type!="left" & cens.type!="right" & cens.type!="interval")stop('cens.type must be left, right or interval. Check documentation!')
  
  
 
  if(p.cens>1| p.cens<0)       stop("the percCensu must be between 0 and 1 . Check documentation!")

  MMsimu(m=m,x=x,z=z,tt=tt,nj=nj,beta=beta,sigmae=sigmae,D=D,phi=phi,percCensu=p.cens,cens.type=cens.type)
      
      }
  
