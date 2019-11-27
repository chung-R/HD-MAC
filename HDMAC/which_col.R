#######################################################################
# strsplit for number of columns 
#######################################################################

which_col = function(x){
  x1<-gsub(" ", "", x, fixed = TRUE)
  x2<-strsplit(x1,",")
  x3<-x2[[1]]
  len<-length(x2[[1]])
  cc<-NULL
  for (i in 1:len){
    aa<-gregexpr(pattern ='-',x3[i])
    
    if (aa[[1]][1]== -1){
      as.numeric(x3[i])
      cc<-c(cc,as.numeric(x3[i]))
    }
    else {
      x4<-strsplit(x3[i],"-")
      cc<-c(cc,c(as.numeric(x4[[1]][1]):as.numeric(x4[[1]][2])))
    }
  }
  return(cc)
}

