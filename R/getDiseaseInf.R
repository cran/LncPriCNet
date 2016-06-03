getDiseaseInf <-
function(OMIMID){
	DiseaseInfList<-get("DiseaseInfList",envir=.GlobalEnv)
	loci<-which(DiseaseInfList[["OMIMId"]]%in%OMIMID)
	if(length(loci)>=1){
		diseaseInf<-DiseaseInfList[loci,]
		print(paste("There are ",diseaseInf[["DLncNum"]]," seed lncRNA and ",diseaseInf[["DGenNum"]]," seed genes related with this phenotype"))
		DiseaseInf<-list()
		DiseaseInf[[1]]<-as.character(diseaseInf[["OMIMName"]])
		DiseaseInf[[2]]<-unlist(strsplit(as.character(diseaseInf[["DLncs"]]),";"))
		DiseaseInf[[3]]<-unlist(strsplit(as.character(diseaseInf[["DGenes"]]),";"))
		names(DiseaseInf)<-c("OMIM Name","LncRNAs","Genes")
	}else{print("The OMIM Id you input is not in our disease list.");DiseaseInf<-NA}
	return(DiseaseInf)
}
