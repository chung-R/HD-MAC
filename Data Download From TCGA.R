#### library packages ####
install.packages("cgdsr")
library("cgdsr")
symbols = as.vector(read.csv("symbols.csv", header = T)[,1])

#### CGDS objects for download genetic data ####
# 1. CGDS
mycgds           = CGDS("http://www.cbioportal.org/")
# 2. Cancer Study
mycancerstudy    = getCancerStudies(mycgds)[22,1]
# 3. Genetic Variable
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[5,1] # mRNA Expression z-Scores (RNA Seq V2 RSEM)
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[8,1] # Mutations
# 4. Case List
mycaselist       = getCaseLists(mycgds,mycancerstudy)[7,1] # Tumor Samples with mRNA data (RNA Seq V2)
mycaselist       = getCaseLists(mycgds,mycancerstudy)[4,1] # All Seq Tumors
# 5. Clinical Data
myClinicalData   = getClinicalData(mycgds,mycaselist)
write.csv(myClinicalData, "Ova_TCGA_mRNA_Clinical_Data.csv")

#### Download Data ####
#### 1. Use gene symbols and getProfileData() function ####

mRNA.df  = lapply(gene_name[1:length(symbols)], function(name){
  df = getProfileData(mycgds, name, mygeneticprofile, mycaselist)
  return(df)})

#### 2. Final Result ####
# a. Combine the list
mRNA.comb.list = (do.call(cbind,mRNA.df))
# b. Output the data
write.csv(mRNA.comb.list, "Ova_TCGA_mRNA_RNA_Seq_Zscore.csv")

##################################################################################
#### CGDS objects for download genetic data ####
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[14,1] # Mutations
mycaselist       = getCaseLists(mycgds,mycancerstudy)[2,1] # All Seq Tumors
myClinicalData   = getClinicalData(mycgds,mycaselist)
write.csv(myClinicalData, "Ova_TCGA_mutations_Clinical_Data.csv")

muta.df  = lapply(gene_name[1:length(symbols)], function(name){
  df = getProfileData(mycgds, name, mygeneticprofile, mycaselist)
  return(df)})

muta.comb.list = (do.call(cbind,muta.df))
write.csv(muta.comb.list, "Ova_TCGA_mutations.csv")