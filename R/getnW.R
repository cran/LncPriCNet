getnW <-
function(x,y,GNet,PNet,LNet,GLNet,PGNet,PLNet){
	z<-1-x-y
	G2L<-GLNet
	#rm(GLNet)
	L2G<-t(G2L)
	rowSums(G2L)->sumG2L
	GLCase<-which(sumG2L!=0)
	G2L[GLCase,]<-y*G2L[GLCase,]/sumG2L[GLCase]

	rowSums(L2G)->sumL2G
	LGCase<-which(sumL2G!=0)
	L2G[LGCase,]<-y*L2G[LGCase,]/sumL2G[LGCase]

	
	P2G<-PGNet
	#rm(PGNet)
	G2P<-t(P2G)
	rowSums(P2G)->sumP2G
	PGCase<-which(sumP2G!=0)
	P2G[PGCase,]<-x*P2G[PGCase,]/sumP2G[PGCase]
	
	rowSums(G2P)->sumG2P
	GPCase<-which(sumG2P!=0)
	G2P[GPCase,]<-x*G2P[GPCase,]/sumG2P[GPCase]

	G<-GNet
	#rm(GNet)
	rowSums(G)->sumG
	sumG[which(sumG==0)]<-1
	#rowSums(G2P)->sumG2P
	rowSums(G2L)->sumG2L
	Gcase1<-which(sumG2P!=0&sumG2L!=0)
	Gcase2<-which(sumG2P==0&sumG2L!=0)
	Gcase3<-which(sumG2P!=0&sumG2L==0)
	Gcase4<-which(sumG2P==0&sumG2L==0)
	G[Gcase1,]<-(1-x-y)*(G[Gcase1,]/sumG[Gcase1])
	G[Gcase2,]<-(1-y)*(G[Gcase2,]/sumG[Gcase2])
	G[Gcase3,]<-(1-x)*(G[Gcase3,]/sumG[Gcase3])
	G[Gcase4,]<-G[Gcase4,]/sumG[Gcase4]


	P2L<-PLNet
	#rm(PLNet)
	L2P<-t(P2L)
	rowSums(P2L)->sumP2L	
	PLCase<-which(sumP2L!=0)
	P2L[PLCase,]<-z*P2L[PLCase,]/sumP2L[PLCase]

	rowSums(L2P)->sumL2P
	LPCase<-which(sumL2P!=0)
	L2P[LPCase,]<-z*L2P[LPCase,]/sumL2P[LPCase]

				

	P<-PNet
	#rm(PNet)
	rowSums(P)->sumP
	sumP[which(sumP==0)]<-1
	#rowSums(P2G)->sumP2G
	#rowSums(P2L)->sumP2L
	Pcase1<-which(sumP2G!=0&sumP2L!=0)
	Pcase2<-which(sumP2G==0&sumP2L!=0)
	Pcase3<-which(sumP2G!=0&sumP2L==0)
	Pcase4<-which(sumP2G==0&sumP2L==0)
	P[Pcase1,]<-(1-z-x)*(P[Pcase1,]/sumP[Pcase1])
	P[Pcase2,]<-(1-z)*(P[Pcase2,]/sumP[Pcase2])
	P[Pcase3,]<-(1-x)*(P[Pcase3,]/sumP[Pcase3])
	P[Pcase4,]<-P[Pcase4,]/sumP[Pcase4]
	
	
	L<-LNet
	#rm(LNet)
	rowSums(L)->sumL
	sumL[which(sumL==0)]<-1
	#rowSums(L2G)->sumL2G
	#rowSums(L2P)->sumL2P
	Lcase1<-which(sumL2G!=0&sumL2P!=0)
	Lcase2<-which(sumL2G==0&sumL2P!=0)
	Lcase3<-which(sumL2G!=0&sumL2P==0)
	Lcase4<-which(sumL2G==0&sumL2P==0)
	L[Lcase1,]<-(1-z-y)*(L[Lcase1,]/sumL[Lcase1])
	L[Lcase2,]<-(1-z)*(L[Lcase2,]/sumL[Lcase2])
	L[Lcase3,]<-(1-y)*(L[Lcase3,]/sumL[Lcase3])
	L[Lcase4,]<-L[Lcase4,]/sumL[Lcase4]
	
###nomalize
	matrixName<-c(colnames(G),colnames(P),colnames(L))
	M1<-cbind(as.matrix(G),as.matrix(G2P),as.matrix(G2L))
	rm(G);rm(G2P);rm(G2L)
	M2<-cbind(as.matrix(P2G),as.matrix(P),as.matrix(P2L))
	rm(P2G);rm(P);rm(P2L)
	M3<-cbind(as.matrix(L2G),as.matrix(L2P),as.matrix(L))
	rm(L2G);rm(L2P);rm(L)
	MM<-rbind(M1,M2,M3)
	rm(M1);rm(M2);rm(M3)
	names(MM)<-matrixName
	rowSums(MM)->sumMM
	MM<-MM/sumMM
	nW<-t(MM)
	rm(MM)
	#print(nW[1:5,1:5])
return(nW)
#  save(nW,file="nW.Rdata")
}
