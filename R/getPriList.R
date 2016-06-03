getPriList <-
function(gamma,x,y,a,b,pheSeed,genSeed,lncSeed,GNet,PNet,LNet,GLNet,PGNet,PLNet)
{

	Ng<-dim(GNet)[1]  #vcount(genNetGraph)
	Np<-dim(PNet)[1]  #vcount(pheNetGraph)
	Nl<-dim(LNet)[1]  #vcount(lncNetGraph)
	PNodes<-colnames(PNet) #get("PNodes",envir=envi)
	GNodes<-colnames(GNet)#get("GNodes",envir=envi)
	LNodes<-colnames(LNet)#get("LNodes",envir=envi)
	v0<-rep(0,length=Np) ;names(v0)<-PNodes
	u0<-rep(0,length=Ng) ; names(u0)<-GNodes
	w0<-rep(0,length=Nl) ;names(w0)<-LNodes
	if(length(pheSeed)>0){v0[which(names(v0)%in%as.character(pheSeed))]<-1}else{v0<-v0;print("Warning! Please input phenotype seeds.")}
	if(length(genSeed)>0){u0[which(names(u0)%in%as.character(genSeed))]<-1}else{u0<-u0;print("Warning! Please input gene seeds.")}
	if(length(lncSeed)>0){w0[which(names(w0)%in%lncSeed)]<-1}else{w0<-w0;print("Warning! Please input lncRNA seeds.")}
	###nomalize initial vector
	if(sum(u0)!=0){u0<-u0/sum(u0)}##new
	if(sum(v0)!=0){v0<-v0/sum(v0)}##new
	if(sum(w0)!=0){w0<-w0/sum(w0)}##new
	P0<-c(a*u0,b*v0,(1-a-b)*w0) #initial vector
	rm(u0);rm(v0);rm(w0)
	P0<-P0/sum(P0)
	GNet<-as.matrix(GNet);PNet<-as.matrix(PNet);LNet<-as.matrix(LNet);
	GLNet<-as.matrix(GLNet);PGNet<-as.matrix(PGNet);PLNet<-as.matrix(PLNet);
	nW<-getnW(x,y,GNet,PNet,LNet,GLNet,PGNet,PLNet) ##transition matrix
	#rm(GNet);rm(PNet);rm(LNet);rm(PGNet);rm(PLNet);rm(GLNet)
	##rwr

	PT<-rwr(nW,P0,gamma)
	rm(nW)
	PT<-as.vector(PT)
	names(PT)<-names(P0)

	#pg<-PT[1:Ng] ##the probabilities of all genes
	#pp<-PT[(Ng+1):(Ng+Np)] ##the probabilities of all phenotypes
	pm<-PT[(Ng+Np+1):(Ng+Np+Nl)] ##the probabilities of all lncRNA
	AllLncRNAsProbs<-pm[order(pm,decreasing=TRUE)]
	r1<-which(names(pm) %in% lncSeed)
	visProbs2<-pm[-r1]
	AllLProbsNoSeedNode<-visProbs2[order(visProbs2,decreasing=TRUE)] ##new   not include seeds

	Attrlist<-list()
	
	Attrlist[[1]]<-AllLncRNAsProbs
	Attrlist[[2]]<-AllLProbsNoSeedNode
	names(Attrlist)<-c("AllLncRNAsProbs","AllLProbsNoSeedNode")
	return(Attrlist)
}
