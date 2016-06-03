\name{getTopDiseaseLncRNAs}
\alias{getTopDiseaseLncRNAs}
\title{Get the Disease Risk LncRNAs}
\description{
  prioritize the disease candidate lncRNAs by integrated multi-omics information.
}
\usage{
    getTopDiseaseLncRNAs(pheSeed = NULL, genSeed = NULL, lncSeed = NULL, 
			candidates = NULL, showTop = 30, gamma = 0.7, x = 1/3, y = 1/3, 
			a = 1/3, b = 1/3, GNet, PNet, LNet, GLNet, PGNet, PLNet)
}
\arguments{  
  \item{pheSeed}{A character vector, represents the disease users want to study. Users should input diseaseName.}
  \item{genSeed}{ A character vector, The gene seeds are the known disease genes of corresponding phenotype. Users should input a gene seed vector.}
  \item{lncSeed}{ A character vector, The lncRNA seeds are the known disease lncRNAs of corresponding phenotype. Users should input a lncRNA seed vector.}
  \item{candidates}{A character vector. Users should input lncRNA candidates, a lncRNA seed vector. }
  \item{showTop}{An integer. The number of top ranked candidate lncRNAs users want to show.}
  \item{gamma}{Restart probability in RWR method. A value ranges from 0 to 1. The default value is 0.7. }
  \item{x}{Jumping probability between gene network and phenotype network. A value ranges from 0 to 1. The default value is 1/3.}
  \item{y}{Jumping probability between gene network and lncRNA network. A value ranges from 0 to 1. The default value is 1/3.}
  \item{a}{a denote the importance of the gene network. A value ranges from 0 to 1. The default value is 1/3.}
  \item{b}{b denote the importance of the phenotype network. A value ranges from 0 to 1. The default value is 1/3.}
  \item{GNet}{An adjacency matrix of the gene network with colnames and rownames representing genes. The value in it ranges from 0 to 1. }
  \item{PNet}{An adjacency matrix of the phenotype network with colnames and rownames representing phenotypes. The value in it ranges from 0 to 1.}
  \item{LNet}{An adjacency matrix of the lncRNA network with colnames and rownames representing lncRNAs. The value in it ranges from 0 to 1.}
  \item{GLNet}{An adjacency matrix of the  gene-lncRNA network with colnames representing lncRNAs and rownames representing genes. The value in it ranges from 0 to 1. }
  \item{PGNet}{An adjacency matrix of the phenotype-gene network with colnames representing genes and rownames representing phenotypes. The value in it ranges from 0 to 1. }
  \item{PLNet}{An adjacency matrix of the phenotype-lncRNA network with colnames representing lncRNAs and rownames representing phenotypes. The value in it ranges from 0 to 1.}
  }
\details{
The function \code{getTopDiseaseLncRNAs} could get the top ranked candidate lncRNAs by prioritization of the disease candidate lncRNAs using known disease seed nodes to perform random walk on multi-omics network.
Note:
If users want to use seeds we provided, pheSeed, genSeed and lncSeed can be obtained by \code{getDiseaseInf}.

}

\value{
  A dataframe.
}


\author{ Qianlan Yao <yaoqianlan@yahoo.com> }
\seealso{
\code{\link{DiseaseInfList}},	\code{\link{getCoExpGeneofLncs}},	\code{\link{getDiseaseInf}}.

}

\examples{
\dontrun{
#############    Prioritize candidate lncRNAs   ################ 

## Here we use six fake matrix(networks).  
data("GNetExample")
data("GLNetExample")
data("PNetExample")
data("LNetExample")
data("PGNetExample")
data("PLNetExample")

## Get three types seeds from DiseaseInfList we provided. 
diseaseName<-"BREAST CANCER" ;
data("DiseaseInfList")
loci<-match(diseaseName,DiseaseInfList[["OMIMName"]]);
pheSeedExample<-DiseaseInfList[loci,"OMIMId"];
genSeedExample<-unlist(strsplit(as.character(DiseaseInfList[loci,"DGenes"]),";"));
lncSeedExample<-unlist(strsplit(as.character(DiseaseInfList[loci,"DLncs"]),";"));
PNodes<-colnames(PNetExample) ;
GNodes<-colnames(GNetExample);
LNodes<-colnames(LNetExample);
pheSeedResult<-pheSeedExample;
genSeedResult<-genSeedExample;
lncSeedResult<-lncSeedExample;

## Get candidates.
candidateResult<-setdiff(LNodes,lncSeedResult);
##Run getTopDiseaseLncRNAs
example<-getTopDiseaseLncRNAs(pheSeed=pheSeedResult,genSeed=genSeedResult,lncSeed=lncSeedResult,
       candidates=candidateResult,showTop=30,GNet=GNetExample,PNet=PNetExample,LNet=LNetExample,
	   GLNet=GLNetExample,PGNet=PGNetExample,PLNet=PLNetExample);
head(example)


}
}
\keyword{ file }

