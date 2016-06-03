\name{DiseaseInfList}
\alias{DiseaseInfList}
\docType{data}
\title{
	All Known Disease Information.
}
\description{
	A dataframe contains known disease information, including disease phenotype, disease genes and disease lncRNAs.
}
\usage{data("DiseaseInfList")}
\format{
  A data frame with 53 observations on the following 6 variables.
  \describe{
    \item{\code{OMIMId}}{a numeric vector}
    \item{\code{DLncs}}{known disease lncRNAs}
    \item{\code{DGenes}}{known disease genes}
    \item{\code{DGenNum}}{a numeric vector}
    \item{\code{DLncNum}}{a numeric vector}
    \item{\code{OMIMName}}{disease phenotype name}
  }
}
\examples{
data(DiseaseInfList)
head(DiseaseInfList)
## maybe str(DiseaseInfList) ; 
}
\keyword{datasets}
