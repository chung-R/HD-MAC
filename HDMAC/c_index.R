############ needed information
## u5 = data.frame(Death = delta,
##                 Death_surtime = time)
## cut_hazard1 

c.index = function(u5,cut_hazard1){
  c11=0;
  for(i in 1:(nrow(u5)-1)){
    for(j in (i+1):nrow(u5)){
      if((u5$Death[i]==1) && (u5$Death[j]==1))
        c11=c11+1
      if((u5$Death[i]==1) && (u5$Death[j]==0) && (u5$Death_surtime[i] < u5$Death_surtime[j]))
        c11=c11+1
      if((u5$Death[i]==0) && (u5$Death[j]==1) && (u5$Death_surtime[j] < u5$Death_surtime[i]))
        c11=c11+1
    }
  }
  c12=0;
  for(i in 1:(nrow(u5)-1)){
    for(j in (i+1):nrow(u5)){
      if((u5$Death[i]==1) && (u5$Death[j]==1) && (u5$Death_surtime[i] < u5$Death_surtime[j]) && (cut_hazard1[i] > cut_hazard1[j]))
        c12=c12+1
      if((u5$Death[i]==1) && (u5$Death[j]==1) && (u5$Death_surtime[i] < u5$Death_surtime[j]) && (cut_hazard1[i] == cut_hazard1[j]))
        c12=c12+0.5
      if((u5$Death[i]==1) && (u5$Death[j]==1) && (u5$Death_surtime[j] < u5$Death_surtime[i]) && (cut_hazard1[j] > cut_hazard1[i]))
        c12=c12+1
      if((u5$Death[i]==1) && (u5$Death[j]==1) && (u5$Death_surtime[j] < u5$Death_surtime[i]) && (cut_hazard1[j] == cut_hazard1[i]))
        c12=c12+0.5
      if((u5$Death[i]==1) && (u5$Death[j]==0) && (u5$Death_surtime[i] < u5$Death_surtime[j]) && (cut_hazard1[i] > cut_hazard1[j]))
        c12=c12+1
      if((u5$Death[i]==1) && (u5$Death[j]==0) && (u5$Death_surtime[i] < u5$Death_surtime[j]) && (cut_hazard1[i] == cut_hazard1[j]))
        c12=c12+0.5
      if((u5$Death[i]==0) && (u5$Death[j]==1) && (u5$Death_surtime[j] < u5$Death_surtime[i]) && (cut_hazard1[j] > cut_hazard1[i]))
        c12=c12+1
      if((u5$Death[i]==0) && (u5$Death[j]==1) && (u5$Death_surtime[j] < u5$Death_surtime[i]) && (cut_hazard1[j] == cut_hazard1[i]))
        c12=c12+0.5
    }
  }
  return(concordance_index = c12/c11)
}
