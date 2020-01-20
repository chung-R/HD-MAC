#######################################################################
# strsplit for number of columns 
#######################################################################

which_col = function(x){
  x1<-gsub(" ", "", x, fixed = TRUE)
  if (strsplit(x1,";") == x1) {
    x2<-strsplit(x1,",")[[1]]
    len<-length(x2)
    cc<-NULL
    for (i in 1:len){
      aa<-gregexpr(pattern ='-',x2[i])
      if (aa[[1]][1]== -1){
        as.numeric(x2[i])
        cc<-c(cc,as.numeric(x2[i]))
      }
      else {
        x3<-strsplit(x2[i],"-")
        cc<-c(cc,c(as.numeric(x3[[1]][1]):as.numeric(x3[[1]][2])))
      }
    }
  }
  else{
    cc<-strsplit(x1,";")[[1]]
  }
  return(cc)
}

