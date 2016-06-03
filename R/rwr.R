rwr <-
function(W,P0,gamma){
	 PT<-as.matrix(P0)
	 k<-0
	 delta<-1
    while(delta>1e-10){
		PT1 = (1-gamma)*W%*%PT+gamma*P0;
        delta = sum(abs(PT1 - PT));
		PT = PT1;
        k = k + 1;
		#print(paste("k and delta:",k,delta))
    }


    PT<-t(PT)
    return(PT)
}
