
# code to generate and characterize human Tissue-specific BioInt libraries

load(file="hsap_apid1.RData")
load(file="tissexp_list.RData")

names(tissexp_list)
## ubiq genes must be counted excluding sexual tissues
tiss_noteval_inubiq=c("breast", "cervix, uterine", "endometrium", "fallopiantube", "ovary", "placenta" , "prostate", "seminal vesicle", "testis")
temp_tissexp_df=tissexp_list[-c(which(names(tissexp_list) %in% tiss_noteval_inubiq))]
ubiq=names(which(table(unlist(temp_tissexp_df)) == length(temp_tissexp_df)))


library("igraph")
library("GOfuncR")
library("philentropy")

hsap_FEA_FunMod_list=list()
hsap_processing_FunMod_list=list()

for ( i in 1:length(tissexp_list)){ ## HPA_symbol_mean1TMM
  #for ( i in 1:1){
  hsap_FEA_FunMod_output=hsap_FEA_FunMod_FUNC_v1(tissexp_list[[i]], hsap_apid1, FWER_thresh=0.05, num_randomiz=500, specie_db, h_jacc=0.1, output_version="little")
  hsap_processing_FunMod_output=hsap_processing_FunMod_FUNC_v1(hsap_FEA_FunMod_output, ubiq, maincomp_thresh=90, min_size=10, max_size=200)
  
  hsap_FEA_FunMod_list[[i]]=hsap_FEA_FunMod_output
  hsap_processing_FunMod_list[[i]]=hsap_processing_FunMod_output
  
}
names(hsap_FEA_FunMod_list)=names(tissexp_list)
names(hsap_processing_FunMod_list)=names(tissexp_list)

hsap_eval_FunMod_list=hsap_eval_FunMod_FUNC_v1(hsap_FEA_FunMod_list, hsap_processing_FunMod_list, filtered="yes")
