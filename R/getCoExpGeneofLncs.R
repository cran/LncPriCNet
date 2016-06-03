getCoExpGeneofLncs <-
function(lncRNAs,GLNet,scoreCutoff=0.6){
	GLNet[,as.character(lncRNAs)]->a
	gene<-names(a[a>scoreCutoff])
	return(gene)
}
