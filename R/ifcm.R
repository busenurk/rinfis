#' Intuitionistic Fuzzy C-Means Clustering
#'
#' Intuitionistic fuzzy c-means clustering algorithm is an extension of fuzzy c-means clustering algorithm.
#'
#'@param x The data matrix where columns correspond to variables and rows to observations.
#'@param c Number of clusters
#'@param m A number greater than 1 giving the degree of fuzzification.
#'@param alpha degree of hesitation
#'@param maxitr Maximum number of iterations.
#'@param epsilon Error value
#'@param fgen The type of intuitionistic fuzzy generator. \code{Yager} or \code{Sugeno} type fuzzy generators
#'can be selected. (Default="Yager")
#'@param lambda The constant
#'@param verbose If TRUE, make some output during learning.
#'@author Erol Egrioglu, Eren Bas and Busenur Kizilaslan \cr Maintainer: Busenur Kizilaslan
#'\email{busenur.sarica@@gmail.com}
#'@references
#'\itemize{
#'\item Chaira, T. (2011)
#'\emph{A novel intuitionistic fuzzy C means clustering algorithm and its application
#' to medical images, Applied Soft Computing, Vol. 11(2), 1711-1717},
#'\url{https://www.sciencedirect.com/science/article/pii/S1568494610001067}.
#'\item Chaira, T., Ray, A. K., Salvetti, O. (2006)
#'\emph{Intuitionistic fuzzy C means clustering in
#'medical images, Advances in Pattern Recognition, 226-230},
#'\url{https://www.worldscientific.com/doi/10.1142/9789812772381_0037}.}\cr

#'@details The data given by x is clustered by generalized versions of the fuzzy c-means algorithm,
#'If verbose is TRUE, each iteration displays its number and the value of the objective function.
#'The parameters m defines the degree of fuzzification. It is defined for real values greater
#'than 1 and the bigger it is the more fuzzy the membership values of the clustered data points are.
#'@importFrom stats runif
#'@importFrom dplyr case_when
#'@return An object of class "ifcm" which is a list with components:\cr\cr
#'\item{membership}{a matrix with the membership values of the data points to the clusters.}
#'\item{hesitation}{a matrix with the hesitation values of the data points to the clusters.}
#'\item{nonmembership}{a matrix with the nonmembership values of the data points to the clusters.}
#'\item{centers}{the final cluster centers.}
#'\item{hardcluster}{the hardclustering results}
#'@export ifcm
#'@examples
#'ifcm(mtcars$mpg,2)
ifcm <- function(x,c,m=2,alpha=0.85,maxitr=100,epsilon=0.03,fgen=c("Yager","Sugeno"),lambda=2, verbose=FALSE){
x<-as.matrix(x)
  n=nrow(x)
  pi=matrix(NA, ncol=c,nrow=n)
  ustar=matrix(NA, ncol=c, nrow=n)
  r=matrix(runif(n*c,min=0,max=1), ncol=c, nrow=n)
  u=matrix(NA, ncol=c, nrow=n)

  if(sum(is.na(x))>0)
  stop("Check the missing value.")

  if(missing(c))
    stop("Argument 'c' must be a number.")

  if(alpha < 0)
    stop("Argument 'alpha' must be positive.")

  if(maxitr < 1)
    stop("Argument 'maxitr' must be positive.")


  for (i in 1:n){
    for (j in 1:c){
      u[i,j] = r[i,j]/sum(r[i,])}}

      fgen=match.arg(fgen)

      if(fgen=="Yager"){ pi=1-u-(1-u^alpha)^(1/alpha)
      } else if(fgen=="Sugeno"){pi=1-u-((1-u)/(1+lambda*u))
      } else{stop("Argument 'fgen' must be selected as Yager or Sugeno")}

  for (i in 1:n){
    for (j in 1:c){
      ifelse(u[i,j]+2*pi[i,j]>1, ustar[i,j]<-u[i,j], ustar[i,j]<-u[i,j]+pi[i,j])
    }
  }

  Uold=ustar

  for (v in 1:maxitr){
    ustar=Uold
    #kume merkezleri
    vstar=matrix(NA, ncol=ncol(x), nrow=c)
    c1=ustar^m



    for (j in 1:c){
      c2=matrix(NA,ncol=ncol(x),nrow=n)
      for (i in 1:n){
        c2[i,]=c1[i,j]*x[i,]
      }

      vstar[j,]=colSums(c2)/sum(c1[,j])}

    d=matrix(NA,ncol=n, nrow=c)
    Uup=matrix(NA, ncol=c, nrow=n)


    for(j in 1:c){
      for (i in 1:n){
        d[j,i]=sqrt(sum((x[i,]-vstar[j,])^2))
        d.square=d^(2/(m-1))}}

    for (i in 1:n){
      for(j in 1:c){
        t1=d.square[j,i]/d.square[,i]

        #guncellenmis uyelikler
        Uup[i,j]=1/sum(t1)
      }
    }

    newpi=matrix(NA, ncol=c,nrow=n)
    newustar=matrix(NA, ncol=c, nrow=n)

    for (i in 1:n){
      for(j in 1:c){

        if (fgen=="Yager") {
          newpi[i,j]=1-Uup[i,j]-(1-Uup[i,j]^alpha)^(1/alpha)
        } else { newpi[i,j]=1-Uup[i,j]-((1-Uup[i,j])/(1+lambda*u[i,j]))}

        ifelse(Uup[i,j]+2*newpi[i,j]>1, newustar[i,j]<-Uup[i,j], newustar[i,j]<-Uup[i,j]+newpi[i,j])
      }
    }

    Unew=newustar


    stop.cri<-0

    stop.cri[v]=sqrt(sum((Uold-Unew)^2))

    if(verbose==TRUE)
    cat("\nIteration", v ,":",stop.cri[v])

    if(stop.cri[v]>epsilon)
    {Uold<-Unew}
    else
    {break}
  }

nm=1-Unew-newpi

hardc <- apply(Unew, 1, which.max)

cat("\n\nIntuitionistic Fuzzy c-means clustering with", c , "clusters:\n\n")

output<-list(membership=Unew, hesitation=newpi, nonmembership=nm, centers=vstar, hardcluster=hardc)
return(output)
}


