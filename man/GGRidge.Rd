\name{GGRidge}
\alias{GGRidge}
\title{
GGRidge: Graphical group ridge
}
\description{
The function has two main tasks. The first task is to find the graph structure ( using GraphicalLasso) which is used to identify the
conditionally correlated groups of predictors and then search for the optimal values of the  tuning parameters to estimate the GGRidge regression coeffcients. For given data, an estimated regression coeffcients, mean squared error (MSE) and the shrinkage parameters for the estimated groups of predictors are returned.
}
\usage{
GGRidge(data, kg=5, sq=c(0.01,5,0.01), k=5, PE=50)
}

\arguments{
  \item{data}{
A matrix of dimension \code{n} x  \code{v)} where \code{v=(p+1)} consisting of the values of the predictors and the response variable  (last variable).
}
  
  \item{kg}{
Number of splits in k-fold cross-validation to find the graph structure. Default value is kg=5.
}
  \item{sq}{
The starting (maximal) end values and number: increment of the sequence of the tuning parameters for the generated penalty term. Defaults are 0.01, 5 & 0.01 respectively.
}
  \item{k}{
Number of splits in k-fold cross-validation to find the ridge penalty parameter. Default value is k=5.
}
  \item{PE}{
 The number of iterations used in computing the prediction mean square error for each group of predictors. The default is PE=50.
}
 
}

\value{
A list of three values:
 \item{Coefficients}{ The vector of estimated regression parameters.}
\item{MSE}{ Mean square error for the parameters.}
\item{lambda.opt}{ Optimal values of the penalty.}

}
\references{


Claus Dethlefsen and Soren Hojsgaard (2005): A Common Platform for Graphical Models in {R}: The {gRbase} Package,
\emph{Journal of Statistical Software}, \url{https://www.jstatsoft.org/v14/i17/}, 14(17).

Gabor Csardi and Tamas Nepusz (2006): The igraph software package for complex network research, \emph{Inter Journal}, \url{https://igraph.org}.

Saeed Aldahmani and Taoufik Zoubeidi (2020): Graphical group ridge, \emph{ Journal of Statistical Computation and Simulation}.

Matt Galloway (2018):CVglasso: Lasso Penalized Precision Matrix Estimation, \url{https://CRAN.R-project.org/package=CVglasso}.

Scheetz, T.E., Kim, K.Y.A., Swiderski, R.E., Philp, A.R., Braun, T.A., Knudtson, K.L., Dorrance, A.M., DiBona, G.F., Huang, J., Casavant, T.L. and Sheffield, V.C. (2006). Regulation of gene expression in the mammalian eye and its relevance to eye disease. \emph{ Proceedings of the National Academy of Sciences}.
}
\author{
Saeed Aldahmani and Taoufik Zoubeidi
}


\examples{

###################################

# Example 

###################################
library(gRbase)
library("CVglasso")
library("MASS")
#load the data 
data(SampleEyeData)
x <- scale(SampleEyeData[,-51], center = TRUE, scale = TRUE) 
y <- scale(SampleEyeData[,51], center = TRUE)

Data<-as.matrix(cbind(x,y))



# Call GGRidge(.) .

results <- GGRidge(data=Data,kg=3,sq=c(0.1,.5,0.1),k=3,PE=30)


# The returned objects are

results$Coefficients
results$MSE
results$lambda.opt


}

\keyword{ package }
