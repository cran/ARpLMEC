% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ARpMMEC.est.R
\name{ARpMMEC.est}
\alias{ARpMMEC.est}
\title{Censored Mixed-Effects Models with Autoregressive Correlation Structure and DEC for Normal and t-Student Errors}
\usage{
ARpMMEC.est(
  y,
  x,
  z,
  tt,
  cc,
  nj,
  struc = "UNC",
  order = 1,
  initial = NULL,
  nu.fixed = TRUE,
  typeModel = "Normal",
  cens.type = "left",
  LI = NULL,
  LS = NULL,
  MaxIter = 200,
  error = 1e-04,
  Prev = FALSE,
  step = NULL,
  isubj = NULL,
  xpre = NULL,
  zpre = NULL
)
}
\arguments{
\item{y}{Vector \code{1 x n} of censored responses, where \code{n} is the sum of the number of observations of each individual}

\item{x}{Design matrix of the fixed effects of order \code{n x s}, corresponding to vector of fixed effects.}

\item{z}{Design matrix of the random effects of order\code{n x b}, corresponding to vector of random effects.}

\item{tt}{Vector \code{1 x n} with the time the measurements were made, where \code{n} is the total number of measurements for all individuals. Default it's considered regular times.}

\item{cc}{Vector of censoring indicators of length \code{n}, where \code{n} is the total of observations. For each observation: \code{0} if non-censored, \code{1} if censored.}

\item{nj}{Vector \code{1 x m} with the number of observations for each subject,  where \code{m} is the total number of individuals.}

\item{struc}{\code{UNC},\code{ARp},\code{DEC},\code{SYM} or \code{DEC(AR)} for uncorrelated ,autoregressive, DEC(phi1,phi2), DEC(phi1,phi2=1), DEC(DEC(phi1,phi2=1)) structure, respectively}

\item{order}{Order of the autoregressive process. Must be a positive integer value.}

\item{initial}{List with the initial values in the next orden: betas,sigma2,alphas,phi and nu. If it is not indicated it will be provided automatically. Default is \code{NULL}}

\item{nu.fixed}{Logical. Should estimate the parameter "nu" for the t-student distribution?. If is False indicates the value in the list of initial values. Default is \code{FALSE}}

\item{typeModel}{\code{Normal} for Normal distribution and \code{Student} for t-Student distribution. Default is \code{Normal}}

\item{cens.type}{\code{left} for left censoring, \code{right} for right censoring and \code{interval} for intervalar censoring. Default is \code{left}}

\item{LI}{Vector censoring lower limit indicator of length \code{n}. For each observation: \code{0} if non-censored, \code{-inf} if censored. It is only indicated for when \code{cens.type} is \code{both}. Default is \code{NULL}}

\item{LS}{Vector censoring upper limit indicator of length \code{n}. For each observation: \code{0} if non-censored, \code{inf} if censored.It is only indicated for when \code{cens.type} is \code{both}. Default is \code{NULL}}

\item{MaxIter}{The maximum number of iterations of the EM algorithm. Default is \code{200}}

\item{error}{The convergence maximum error. Default is \code{0.0001}}

\item{Prev}{Indicator of the prediction process. Available at the moment only for the \code{typeModel=normal} case.  Default is \code{FALSE}}

\item{step}{Number of steps for prediction. Default is \code{NULL}}

\item{isubj}{Vector indicator of subject included in the prediction process. Default is \code{NULL}}

\item{xpre}{Design matrix of the fixed effects to be predicted. Default is \code{NULL}.}

\item{zpre}{Design matrix of the random effects to be predicted. Default is \code{NULL}.}
}
\value{
returns list of class \dQuote{ARpMMEC}:
\item{FixEffect}{Data frame with: estimate, standar errors and confidence intervals of the fixed effects.}
\item{Sigma2}{Data frame with: estimate, standar errors and confidence intervals  of the variance of the white noise process.}
\item{Phi}{Data frame with: estimate, standar errors and confidence intervals  of the autoregressive parameters.}
\item{RandEffect}{Data frame with: estimate, standar errors and confidence intervals  of the random effects.}
\item{nu}{the parameter "nu" for the t-student distribution}
\item{Est}{Vector of parameters estimate (fixed Effects, sigma2, phi, random effects).}
\item{SE}{Vector of the standard errors of (fixed Effects, sigma2, phi, random effects).}
\item{Residual}{Vector of the marginal residuals.}
\item{loglik}{Log-likelihood value.}
\item{AIC}{Akaike information criterion.}
\item{BIC}{Bayesian information criterion.}
\item{AICc}{Corrected Akaike information criterion.}
\item{iter}{Number of iterations until convergence.}
\item{Yfit}{Vector "y" fitted}
\item{MI}{Information matrix}
\item{Prev}{Predicted values (if xpre and zpre is not \code{NULL}).}
\item{time}{Processing time.}
\item{others}{The first and second moments of the random effect and vector Y}
}
\description{
This functino fits left, right or intervalar censored mixed-effects linear model, with autoregressive errors of order \code{p}, using the EM algorithm. It returns estimates, standard errors and prediction of future observations.
}
\examples{
\dontrun{
p.cens   = 0.1
m           = 10
D = matrix(c(0.049,0.001,0.001,0.002),2,2)
sigma2 = 0.30
phi    = 0.6
beta   = c(1,2,1)
nj=rep(4,10)
tt=rep(1:4,length(nj))
x<-matrix(runif(sum(nj)*length(beta),-1,1),sum(nj),length(beta))
z<-matrix(runif(sum(nj)*dim(D)[1],-1,1),sum(nj),dim(D)[1])
data=ARpMMEC.sim(m,x,z,tt,nj,beta,sigma2,D,phi,struc="ARp",typeModel="Normal",p.cens=p.cens)

teste1=ARpMMEC.est(data$y_cc,x,z,tt,data$cc,nj,struc="ARp",order=1,typeModel="Normal",MaxIter = 2)
teste2=ARpMMEC.est(data$y_cc,x,z,tt,data$cc,nj,struc="ARp",order=1,typeModel="Student",MaxIter = 2)

xx=matrix(runif(6*length(beta),-1,1),6,length(beta))
zz=matrix(runif(6*dim(D)[1],-1,1),6,dim(D)[1])
isubj=c(1,4,5)
teste3=ARpMMEC.est(data$y_cc,x,z,tt,data$cc,nj,struc="ARp",order=1,typeModel="Normal",
                  MaxIter = 2,Prev=TRUE,step=2,isubj=isubj,xpre=xx,zpre=zz)
teste3$Prev

}


}
\references{
Olivari, R. C., Garay, A. M., Lachos, V. H., & Matos, L. A. (2021). Mixed-effects 
models for censored data with autoregressive errors. Journal of Biopharmaceutical Statistics, 31(3), 273-294.
\doi{10.1080/10543406.2020.1852246}
}
