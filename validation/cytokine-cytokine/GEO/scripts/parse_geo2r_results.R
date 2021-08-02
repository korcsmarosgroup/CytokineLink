library(tidyverse)
setwd('~/OneDrive - Norwich BioScience Institutes/Documents/CytokineLink/analysis/cytcyt_validation/GEO/')

parseGEO2r<-function(infile,cytokine){
  CKL<-read_csv('data/cytokine2cytokine_simplified.csv', show_col_types = F) #read predicted interactions
  cytokines<- c(CKL$Source_cytokine_genename, CKL$Target_cytokine_genename) %>% unique()
  TNF<- CKL %>% filter(Source_cytokine_genename == paste0(cytokine)) #get blocked cytokine interactions
  experiment<-read_tsv(paste0('data/',infile), show_col_types = F) #expression data
  all_cyt_in_dataset<- experiment %>% separate_rows(Gene.symbol, sep = '///') %>% filter(Gene.symbol %in% cytokines) #all cytokines in experiment
  experiment_sign<-experiment %>% separate_rows(Gene.symbol, sep = '///') %>% filter(Gene.symbol %in% cytokines) %>% filter(logFC >= 2 | logFC <= -2) %>% filter(adj.P.Val <= 0.05)#filt for significance/FC
  TNF<- TNF %>% filter(Target_cytokine_genename %in% all_cyt_in_dataset$Gene.symbol) #cutting down targets to present ones
  overlap<- TNF %>% filter(Target_cytokine_genename %in% experiment_sign$Gene.symbol) #TNF targeted and sign
  
  #2x2 for chi square
  cklDEG<-dim(overlap)[1] #targeted by blocked cytokine AND significant
  cklAll<-dim(TNF)[1] - dim(overlap)[1] #targeted by blocked cytokine
  allDEG <- setdiff(experiment_sign$Gene.symbol,overlap$Target_cytokine_genename) 
  allDEG <- length(allDEG) #significant cytokine
  allc <- length(unique(all_cyt_in_dataset$Gene.symbol)) - cklDEG - allDEG -cklAll#all cytokine
  x <- matrix(c(cklDEG, cklAll, allDEG,allc), byrow = F, 2, 2)
  print(infile)
  print(x)
  print(chisq.test(x))
}
experiments<-list.files('data/', pattern = 'GSE*')
parseGEO2r('GSE92415.top.table.tsv','TNF') #golimumab
for (j in experiments){
  parseGEO2r(j, 'TNF')
}

#manual run for tocilizumab vs healthy control
CKL<-read_csv('data/cytokine2cytokine_simplified.csv', show_col_types = F) #read predicted interactions
cytokines<- c(CKL$Source_cytokine_genename, CKL$Target_cytokine_genename) %>% unique()
TNF<- CKL %>% filter(Source_cytokine_genename == 'IL6') #get blocked cytokine interactions
experiment<-read_tsv('data/GSE93777.top.table.tsv', show_col_types = F) #expression data
all_cyt_in_dataset<- experiment %>% separate_rows(Gene.symbol, sep = '///') %>% filter(Gene.symbol %in% cytokines) #all cytokines in experiment
experiment_sign<-experiment %>% separate_rows(Gene.symbol, sep = '///') %>% filter(Gene.symbol %in% cytokines) %>% filter(logFC >= 1.5| logFC <= -1.5) %>% filter(adj.P.Val <= 0.05) #filt for significance/FC
TNF<- TNF %>% filter(Target_cytokine_genename %in% all_cyt_in_dataset$Gene.symbol) #cutting down targets to present ones
overlap<- TNF %>% filter(Target_cytokine_genename %in% experiment_sign$Gene.symbol) #TNF targeted and sign

#2x2 for chi square
cklDEG<-dim(overlap)[1] #targeted by blocked cytokine AND significant
cklAll<-dim(TNF)[1] - dim(overlap)[1] #targeted by blocked cytokine
allDEG <- setdiff(experiment_sign$Gene.symbol,overlap$Target_cytokine_genename) 
allDEG <- length(allDEG) #significant cytokine
allc <- length(unique(all_cyt_in_dataset$Gene.symbol)) - cklDEG - allDEG -cklAll#all cytokine
x <- matrix(c(cklDEG, cklAll, allDEG,allc), byrow = F, 2, 2)
print(infile)
print(x)
print(chisq.test(x))
