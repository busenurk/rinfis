% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ifcm.R
\name{ifcm}
\alias{ifcm}
\title{Intuitionistic Fuzzy C-Means Clustering}
\usage{
ifcm(
  x,
  c,
  m = 2,
  alpha = 0.85,
  maxitr = 100,
  epsilon = 0.03,
  fgen = c("Yager", "Sugeno"),
  lambda = 2,
  verbose = FALSE,
  standardize = TRUE
)
}
\arguments{
\item{x}{The data matrix where columns correspond to variables and rows to observations.}

\item{c}{Number of clusters}

\item{m}{A number greater than 1 giving the degree of fuzzification.}

\item{alpha}{degree of hesitation}

\item{maxitr}{Maximum number of iterations.}

\item{epsilon}{Error value}

\item{fgen}{The type of intuitionistic fuzzy generator. \code{Yager} or \code{Sugeno} type fuzzy generators
can be selected. (Default="Yager")}

\item{lambda}{The constant}

\item{verbose}{If TRUE, make some output during learning.}

\item{standardize}{If TRUE the data matrix is standardized. If variables are in the same units already,
you might not wish to standardize. Default is \code{standardize=TRUE}.}
}
\value{
An object of class "ifcm" which is a list with components:\cr\cr
\item{membership}{a matrix with the membership values of the data points to the clusters.}
\item{hesitation}{a matrix with the hesitation values of the data points to the clusters.}
\item{nonmembership}{a matrix with the nonmembership values of the data points to the clusters.}
\item{centers}{the final cluster centers.}
\item{hardcluster}{the hardclustering results}
}
\description{
Intuitionistic fuzzy c-means clustering algorithm is an extension of fuzzy c-means clustering algorithm.
}
\details{
The data given by x is clustered by generalized versions of the fuzzy c-means algorithm,
If verbose is TRUE, each iteration displays its number and the value of the objective function.
The parameters m defines the degree of fuzzification. It is defined for real values greater
than 1 and the bigger it is the more fuzzy the membership values of the clustered data points are.
}
\examples{
ifcm(mtcars$mpg,2)
}
\references{
\itemize{
\item Chaira, T. (2011)
\emph{A novel intuitionistic fuzzy C means clustering algorithm and its application
to medical images, Applied Soft Computing, Vol. 11(2), 1711-1717},
\url{https://www.sciencedirect.com/science/article/pii/S1568494610001067}.
\item Chaira, T., Ray, A. K., Salvetti, O. (2006)
\emph{Intuitionistic fuzzy C means clustering in
medical images, Advances in Pattern Recognition, 226-230},
\url{https://www.worldscientific.com/doi/10.1142/9789812772381_0037}.}\cr
}
\author{
Erol Egrioglu, Eren Bas and Busenur Kizilaslan \cr Maintainer: Busenur Kizilaslan
\email{busenur.sarica@gmail.com}
}
