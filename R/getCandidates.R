getCandidates <-
function(pheSeed,genSeed,lncSeed,candidates=NULL,GNodes=NULL,LNodes=NULL,PNodes=NULL){
seeds<-getSeed(pheSeed,genSeed,lncSeed,GNodes,LNodes,PNodes)
seeds[[2]]->lncSeeds
if( length(candidates)==0){Candidates<-setdiff(LNodes,lncSeeds)}else{Candidates<-intersect(LNodes,candidates)}
return(Candidates)
}
