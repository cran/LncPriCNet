getSeed <-
function(pheSeed,genSeed,lncSeed,GNodes=NULL,LNodes=NULL,PNodes=NULL){ 
pheSeeds<-intersect(pheSeed,PNodes)
lncSeeds<-intersect(lncSeed,LNodes)
genSeeds<-intersect(genSeed,GNodes)
seeds<-list()
seeds[[1]]<-pheSeeds
seeds[[2]]<-lncSeeds
seeds[[3]]<-genSeeds
names(seeds)<-c("pheSeeds","lncSeeds","genSeeds")
return(seeds)
}
