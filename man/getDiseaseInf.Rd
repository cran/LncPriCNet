\name{getDiseaseInf}
\alias{getDiseaseInf}
\title{Get All Disease Information based on phenotype id (OMIM Id) Provided by LncPriCNet Package}
\description{
  Get all disease information provided by LncPriCNet package by OMIMID. This function will provide you known disease genes and lncRNAs corresponding phenotype.
}

\usage{
     getDiseaseInf(OMIMID)
}

\arguments{  
  \item{OMIMID}{A character. OMIM ID of one disease or phenotype.}
}


\value{
  A list.
}


\author{ Qianlan Yao <yaoqianlan@yahoo.com> }
\seealso{\code{\link{DiseaseInfList}}}

\examples{
\dontrun{
#############    Get disease information   ################
##breast cancer 
diseaseInf<-getDiseaseInf("114480")
diseaseInf
## 
data(DiseaseInfList)
head(DiseaseInfList)


}
}
\keyword{ file }

