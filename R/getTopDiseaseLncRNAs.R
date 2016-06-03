getTopDiseaseLncRNAs <-
function(pheSeed=NULL,genSeed=NULL,lncSeed=NULL
     ,candidates=NULL,showTop=30,gamma=0.7,x=1/3,y=1/3,a=1/3,b=1/3,GNet,PNet,LNet,GLNet,PGNet,PLNet){
		if(x+y>=1){print("The value of x and y range from [0,1]. The sum of x and y must smaller than 1.");break;};
		if(a+b>=1){print("The value of a and b range from [0,1]. The sum of a and b must smaller than 1.");break;};
		PNodes<-colnames(PNet) #get("PNodes",envir=envi)
		GNodes<-colnames(GNet)#get("GNodes",envir=envi)
		LNodes<-colnames(LNet)#get("LNodes",envir=envi)
		seeds<-getSeed(pheSeed,genSeed,lncSeed,GNodes,LNodes,PNodes)
		newpheSeed<-seeds[[1]]
		newlncSeed<-seeds[[2]]
		newgenSeed<-seeds[[3]]
		pheSeedLength<-length(newpheSeed)
		genSeedLength<-length(newgenSeed)
		lncSeedLength<-length(newlncSeed)
		if(sum(pheSeedLength,genSeedLength,lncSeedLength)==0){print("The seeds you input are not in our network. Error! Please re-enter.");break}
		
		if(is.null(newpheSeed)==TRUE){
			pheSeedNodeInf1<-paste("The number of the phenotype seeds you input are: 0")
		}else{pheSeedNodeInf1<-paste("The number of the phenotype seeds you input are: " ,length(pheSeed))}
		pheSeedNodeInf2<-paste("The number of the phenotype seeds used in prioritizing the candidate lncRNA are: " ,length(newpheSeed))
		pheSeedNodeInf3<-paste("The phenotype seeds used in prioritizing the candidate lncRNA are: ",paste(newpheSeed,collapse=";"))
		
		if(is.null(newlncSeed)==TRUE){
			lncSeedNodeInf1<-paste("The number of the lncRNA seeds you input are: 0")			
		}else{lncSeedNodeInf1<-paste("The number of the lncRNA seeds you input are: " ,length(lncSeed))}
		lncSeedNodeInf2<-paste("The number of the lncRNA seeds used in prioritizing the candidate lncRNA are: " ,length(newlncSeed))
		lncSeedNodeInf3<-paste("The lncRNA seeds used in prioritizing the candidate lncRNA are: ",paste(newlncSeed,collapse=";"))
		
		if(is.null(newgenSeed)==TRUE){
			genSeedNodeInf1<-paste("The number of the gene seeds you input are: 0")
		}else{genSeedNodeInf1<-paste("The number of the gene seeds you input are: " ,length(genSeed))}
		genSeedNodeInf2<-paste("The number of the gene seeds used in prioritizing the candidate lncRNA are: " ,length(newgenSeed))
		genSeedNodeInf3<-paste("The gene seeds used in prioritizing the candidate lncRNA are: ",paste(newgenSeed,collapse=";"))		

		Candidates<-getCandidates(pheSeed=newpheSeed,genSeed=newgenSeed,lncSeed=newlncSeed,candidates=candidates,GNodes=GNodes,LNodes=LNodes,PNodes=PNodes)
		if(length(Candidates)==0){print("The candidates you input are not in our network. Error! Please re-enter.");break}
		if(is.null(Candidates)==FALSE){
			CandidatesInf1<-paste("The number of the candidate lncRNA you input are: " ,length(candidates))
		}else{CandidatesInf1<-paste("The number of the candidate lncRNA you input are: " ,0)}
		CandidatesInf2<-paste("The number of lncRNA candidates to be prioritized are " ,length(Candidates))
		DiseaseInf<-rbind(pheSeedNodeInf1,pheSeedNodeInf2,pheSeedNodeInf3,lncSeedNodeInf1,lncSeedNodeInf2,lncSeedNodeInf3,genSeedNodeInf1,genSeedNodeInf2,genSeedNodeInf3,CandidatesInf1,CandidatesInf2)
		DiseaseInf<-unname(DiseaseInf)
		print(DiseaseInf)
		#write.table(DiseaseInf,"DiseaseInf.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
		
		Attrist<-getPriList(gamma,x,y,a,b,pheSeed=newpheSeed,genSeed=newgenSeed,lncSeed=newlncSeed,GNet=GNet,PNet=PNet,LNet=LNet,GLNet=GLNet,PGNet=PGNet,PLNet=PLNet)
		topLncRNAsID<-c()
		topLncRNAsScore<-c()
		
			loci<-which(names(Attrist$AllLncRNAsProbs)%in%Candidates)
			if(length(loci)<=showTop)show<-length(loci)
			if(length(loci)>showTop)show<-showTop
			candidateInf<-Attrist$AllLncRNAsProbs[loci]
			orderedCandidateInf<-candidateInf[order(candidateInf,decreasing=TRUE)]
			topLncRNAsID<-names(orderedCandidateInf)
			topLncRNAsScore<-orderedCandidateInf
			
		Rank<-c(1:length(topLncRNAsID))
		#LncRNAInf<-get("LncRNAInf",envir=envi) ##pubchem
		#topLncRNAsName<-sapply(topLncRNAsID,function(x){as.character(LncRNAInf[which(LncRNAInf[[1]]%in%x),"LncRNAName"][[1]])})
		TopLncRNAInf<-data.frame(Rank,topLncRNAsID,topLncRNAsScore)
		#TopLncRNAInf[order(TopLncRNAInf[["Rank"]]),]->TopLncRNAInf1
		print(show)
		TopLncRNAInf2<-TopLncRNAInf[1:show,]
		#write.table(TopLncRNAInf,"TopLncRNAInf.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
		return(TopLncRNAInf2)
		
		

}
