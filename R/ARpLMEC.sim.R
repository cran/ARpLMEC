#' @title Generating Censored Autoregressive Dataset with Linear Mixed Effects.
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
#' @param nj Vector \code{1 x m} with the number of observations for each subject,  where \code{m} is the total number of individuals.
#' @param beta Vector of values fixed effects.
#' @param sigmae It's the value for sigma. 
#' @param D1 Covariance Matrix for the random effects.
#' @param phi Vector of length \code{Arp}, of values for autoregressive parameters. 
#' @param p.cens Censoring level for the process. Default is \code{0}
#' @param cens.type \code{left} for left censoring, \code{right} for right censoring and \code{interval} for interval censoring. Default is \code{left}
#' @return returns list:
#' \item{cc}{Vector of censoring indicators.}
#' \item{y_cc}{Vector of responses censoring.}
#' @references Schumacher FL, Lachos VH, Dey DK (2017). Censored regression models with autoregressive
#' errors: A likelihood-based perspective. Canadian Journal of Statistics. 
#' \url{https://doi.org/10.1002/cjs.11338}
#' @references Garay AM, Castro LM, Leskow J, Lachos VH (2017). Censored linear regression models for irregularly
#' observed longitudinal data using the multivariate-t distribution. Statistical Methods in
#' Medical Research.
#' \url{https://doi.org/10.1177/0962280214551191}
#' @examples
#' p.cens   = 0.1
#' m           = 50
#' D = matrix(c(0.049,0.001,0.001,0.002),2,2)
#' sigma2 = 0.30
#' phi    = c(0.48,-0.2)
#' beta   = c(1,2,1)
#' nj=rep(6,m) 
#' x<-matrix(runif(sum(nj)*length(beta),-1,1),sum(nj),length(beta))
#' z<-matrix(runif(sum(nj)*dim(D)[1],-1,1),sum(nj),dim(D)[1])
#' data=ARpLMEC.sim(m,x,z,nj,beta,sigma2,D,phi,p.cens)
#' y<-data$y_cc
#' cc<-data$cc
#' 
#' @export
ARpLMEC.sim=function(m,x=NULL,z=NULL,nj,beta,sigmae,D1,phi,p.cens= 0,cens.type="left")
  {
  
  #Control questions for the m parameter
  if(m==sum(nj))                        stop("not compatible sizes between m and nj")

  #Control questions for the x parameter
  if(!is.null(x)){
    if(!is.numeric(x))                    stop("x must be a numeric matrix. Check documentation!")
    if(sum(is.na(x))>0)                   stop("There are some NA values in x.")
    if(!is.matrix(x))                     stop("x must be a matrix. Check documentation!")
    if(det(t(x)%*%x)==0)                  stop("the columns of x must be linearly independent.")
    if(length(x)==0)                      stop("The parameter x must be provided.")
    if(dim(x)[1]!=dim(z)[1])              stop("not compatible sizes between x and z")
    if(dim(x)[1]!=sum(nj))                stop("not compatible sizes between x and nj")
    if(dim(x)[2]!=length(beta))           stop("not compatible sizes between x and beta")
  }
  
  #Control questions for the z parameter
  if(!is.null(z)){
    if(!is.numeric(z))                    stop("z must be a numeric matrix. Check documentation!")
    if(sum(is.na(z))>0)                   stop("There are some NA values in z.")
    if(!is.matrix(z))                     stop("z must be a matrix. Check documentation!")
    if(length(z)==0)                      stop("The parameter z must be provided.")
    if(dim(z)[1]!=sum(nj))                stop("not compatible sizes between z and nj")
    if(dim(z)[2]!=dim(D1)[2])             stop("not compatible sizes between z and D1")
  }
 
  #Control questions for the nj parameter
  if(!is.numeric(nj))                   stop("nj must be a numeric vector. Check documentation!")
  if(!is.vector(nj))                    stop("nj must be a vector. Check documentation!")
  if(sum(is.na(nj))>0)                  stop("There are some NA values in nj")
  if(length(nj)==0)                     stop("The parameter nj must be provided.")
  
  # #Control questions for the Arp parameter
  #  if(length(Arp)!=1)                    stop("Arp must be a value.")
  # if(is.numeric(Arp))
  # { if(Arp!=round(Arp)|Arp<=0)         stop("Arp must be UNC or a positive integer value.")}
  # 
  #Control questions for the Initial parameter
    if(!is.numeric(beta))             stop("beta must be a numeric vector. Check documentation!")
    if(!is.vector(beta))              stop("beta must be a vector. Check documentation!")
    if(!is.numeric(sigmae))           stop("sigmae must be a scalar.")
    if(length(sigmae)>1)              stop("beta must be a scalar.")
    if(!is.matrix(D1))                stop("D1 must be a matrix.")
    if(D1[upper.tri(D1)]!=D1[lower.tri(D1)])stop("D1 must be a simetric matrix.")
    if(!is.numeric(phi))              stop("phi must be a numeric vector. Check documentation!")
  #  if(length(phi)!=Arp)              stop("not compatible sizes between the value Arp and parameter phi. Check documentation!")
    if(sum(abs(phi))>=1)                   stop("the sum of the phi must be less than 1. Check documentation!")
  
  #Control questions for the cens.type parameter
    if(cens.type!="left" & cens.type!="right" & cens.type!="interval")stop('cens.type must be left, right or interval. Check documentation!')
  
  
  #Control questions for the percCensu parameter
  if(p.cens>1| p.cens<0)       stop("the percCensu must be between 0 and 1 . Check documentation!")

  MMsimu(m=m,x=x,z=z,nj=nj,beta=beta,sigmae=sigmae,D1=D1,phi=phi,percCensu=p.cens,cens.type=cens.type)
      
      }
  
