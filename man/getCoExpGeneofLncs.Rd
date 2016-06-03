\name{getCoExpGeneofLncs}
\alias{getCoExpGeneofLncs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get Co-expressed Genes of LncRNAs}
\description{
	Get co-expressed genes of one or some lncRNAs based on the multi-omics network.
}
\usage{
getCoExpGeneofLncs(lncRNAs, GLNet, scoreCutoff = 0.6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lncRNAs}{
	A character or a character vector. The names of lncRNAs.
}
  \item{GLNet}{
	An adjacency matrix of the  gene-lncRNA network with colnames representing lncRNAs and rownames representing genes. The value represents co-expression relation and it ranges from 0 to 1.
}
  \item{scoreCutoff}{
	A numeric donotes if the co-expression scores between genes and lncRNAs above this value, the genes would be retained. A value ranges from 0 to 1. The default value is 0.6.
}
}

\value{
	A vector of gene names. 
}

\author{
	Qianlan Yao <yaoqianlan@yahoo.com>
}

\seealso{
	\code{\link{getTopDiseaseLncRNAs}}
}
\examples{

\dontrun{
#############    getCoExpGeneofLncs   ################
##breast cancer 
data("GLNetExample")
diseaseInf<-getDiseaseInf("114480")
lncRNA<-diseaseInf$LncRNAs[1]
genes<-getCoExpGeneofLncs(lncRNA, GLNet=GLNetExample, scoreCutoff = 0.6) 
print(genes[1:5])

}
}
\keyword{file }
