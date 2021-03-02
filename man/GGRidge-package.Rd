\name{GGRidge-package}
\alias{GGRidge-package}
\docType{package}
\title{GGRidge: Graphical Group Ridge}
\description{
 The Graphical Group Ridge 'GGRidge' classifies ridge regression predictors in disjoint
groups of conditionally correlated variables and derives different
penalties (shrinkage parameters) for these groups of predictors. It
combines the ridge regression method with the graphical model
for high-dimensional data (i.e. the number of predictors exceeds
the number of cases) or ill-conditioned data (e.g. in the presence
of multicollinearity among predictors). The package reduces the mean square errors and the extent of over-shrinking of predictors
as compared to the ridge method. 
}
\details{
\tabular{ll}{
Package: \tab GGRidge\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2021-1-28\cr
License: \tab GPL-2 \cr
}
 
}
\author{
Saeed Aldahmani and Taoufik Zoubeidi

Maintainer: Saeed Aldahmani <saldahmani@uaeu.ac.ae>

}
\references{

Claus Dethlefsen and Soren Hojsgaard (2005): A Common Platform for Graphical Models in {R}: The {gRbase} Package,
\emph{Journal of Statistical Software}, \url{https://www.jstatsoft.org/v14/i17/}, 14(17).

Gabor Csardi and Tamas Nepusz (2006): The igraph software package for complex network research, \emph{Inter Journal}, \url{https://igraph.org}.

Saeed Aldahmani and Taoufik Zoubeidi (2020): Graphical group ridge, \emph{ Journal of Statistical Computation and Simulation}.

Matt Galloway (2018):CVglasso: Lasso Penalized Precision Matrix Estimation, \url{https://CRAN.R-project.org/package=CVglasso}.

Scheetz, T.E., Kim, K.Y.A., Swiderski, R.E., Philp, A.R., Braun, T.A., Knudtson, K.L., Dorrance, A.M., DiBona, G.F., Huang, J., Casavant, T.L. and Sheffield, V.C. (2006). Regulation of gene expression in the mammalian eye and its relevance to eye disease. \emph{ Proceedings of the National Academy of Sciences}.
}
\keyword{ package }



