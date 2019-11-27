library(data.table)
Blad_SUB_8clinical_4937mut_8024mRNA = as.data.frame(fread("Blad_SUB_8clinical_4937mut_8024mRNA.csv",
                                                          header = T,
                                                          sep = ","))
Ova_TCGA_OS_clinical_muta_cleaned_313_13_670 = as.data.frame(fread("Ova_TCGA_OS_clinical_muta_cleaned_313_13_670.csv",
                                                          header = T,
                                                          sep = ","))
