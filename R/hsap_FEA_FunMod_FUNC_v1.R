

# 1: Function to GENERATE BioInt Libraries
  # Reconstruct Tissue-Specific PPI network
  # Perform functional enrichment analysis
  # Merge Gene Ontology terms according to user-defined Jaccard coefficient

# Function to run ITERATIVELY FOR EACH OBJECT in the list of tissue-specific gene list

# @param tissuegenes Vector of tissue-specific genes
# @param ppidata  Data frame of PPI network. Two columns named "from" and "to"
# @param FWER_thresh Numeric value of Family wise error rate threshold for functional enrichment analysis
# @param num_randomiz Numeric value of number of randomizations for functional enrichment analysis
# @param specie_db Package name of the specie database used for Functional enrichment as indicated in GOfuncR R package
# @param h_jacc Numeric value of Jaccard coefficient threshold to merge enriched Gene Ontology terms
# @param output_version Character string to indicate if user wants to retrieve Jaccard merging complete output (Default no, BEWARE: heavy output)

# @return If output_version == "litte":: List including enrichedGOvec, tissue_ppidatat_df, sharedBP_PPIindex_list
  # @return enrichedGOvec Vector of enriched GO terms after Jaccards merging
  # @return tissue_ppidatat_df Tissue-specific PPI network employed for the contruction of TS-BioInt library
  # @return sharedBP_PPIindex_list List of vectors of edge indexes (rows) to extract PPIs from tissue_ppidatat_df for each enriched GO term

# @return If output_version != "litte": List including additionally sharedBP_PPIindex_mat and groupGO_geneass
  # @return sharedBP_PPIindex_mat PPIs for each enriched GO term (heavy output)
  # @return groupGO_geneass gene - goterm group associattions after Jaccards merging

library("GOfuncR")
library("philentropy")


hsap_FEA_FunMod_FUNC_v1=function(tissuegenes, ppidata, FWER_thresh=0.05, num_randomiz=1000, specie_db, h_jacc=0.1, output_version="little"){
  
  allproteins=c(ppidata$from,ppidata$to,tissuegenes)
  allproteins=unique(allproteins)
  
  allPPI_symbol=ppidata
  allproteins_symbol=c(allPPI_symbol$from,allPPI_symbol$to)
  allproteins_symbol=unique(allproteins_symbol)
  
  allnodesdf=data.frame("gene_ids"=allproteins_symbol, "sel"=0, stringsAsFactors = F)
  allnodesdf$sel[which(allnodesdf$gene_ids %in% tissuegenes)]=1
  allnodesdf=allnodesdf[which(!is.na(allnodesdf[,1])),]
  # 1 FEA
  tissue_net_hyper_res= GOfuncR::go_enrich(allnodesdf, n_randset=num_randomiz, orgDb=specie_db$packageName, domains="biological_process")[[1]]
  tissuenet_FEABP=subset(tissue_net_hyper_res,FWER_overrep<=FWER_thresh)
  
  # map again to original gene identifier
  anno_genes_originalID = GOfuncR::get_anno_genes(go_ids=tissuenet_FEABP$node_id,genes=subset(allnodesdf,sel=1)$gene_ids,
                                                  database=specie_db$packageName)
  anno_genes_originalID=cbind(anno_genes_originalID, "GOname"=0)
  anno_genes_originalID$GOname=tissuenet_FEABP$node_name[match(anno_genes_originalID$go_id, tissuenet_FEABP$node_id)]
  
  
  ######## jaccard to merge similar functions
  
  annogenes_m=matrix(0, ncol=length(unique(anno_genes_originalID$GOname)), nrow=length(unique(anno_genes_originalID$gene)))
  colnames(annogenes_m)=unique(anno_genes_originalID$GOname)
  rownames(annogenes_m)=unique(anno_genes_originalID$gene)
  
  for ( i in 1:ncol(annogenes_m)){
    AA=anno_genes_originalID$gene[which(anno_genes_originalID$GOname==colnames(annogenes_m)[i])]
    annogenes_m[which(rownames(annogenes_m) %in% AA),i]=1
  }
  
  JACsim=philentropy::distance(t(annogenes_m), method="jaccard")
  colnames(JACsim)=unique(anno_genes_originalID$GOname)
  rownames(JACsim)=unique(anno_genes_originalID$GOname)
  
  clusters <- hclust(as.dist(JACsim), method = "average")
  group <- cutree(clusters, h = h_jacc)
  groups_vec=(unique(group))
  
  annogenes_mgroups=data.frame()
  for ( i in groups_vec){
    group_nam=names(group)[which(group==i)][which.min(nchar(names(group)[which(group==i)]))]
    annogenes_mgroupssub=data.frame("gene"=unique(anno_genes_originalID$gene[which(anno_genes_originalID$GOname %in% names(which(group==i)))]), "GOname"=group_nam,"go_id"=anno_genes_originalID$go_id[which(anno_genes_originalID$GOname==group_nam)[1]])
    annogenes_mgroups=rbind(annogenes_mgroups, annogenes_mgroupssub)
  }
  
  
  ########
  
  # keep only tissue-specific interactions
  indexes=c(which(is.element(ppidata$from, tissuegenes)),which(is.element(ppidata$to, tissuegenes)))
  indexes=indexes[which(duplicated(indexes))]
  tissue_ppidatat_df=ppidata[indexes,1:2]
  
  ### alternative function - different output format to be faster
  
  enrichedGOvec=unique(annogenes_mgroups$go_id)
  sharedBP_PPIindex_mat=matrix(0, nrow=length(enrichedGOvec), ncol=nrow(tissue_ppidatat_df))
  rownames(sharedBP_PPIindex_mat)=enrichedGOvec
  colnames(sharedBP_PPIindex_mat)=as.character(1:nrow(tissue_ppidatat_df))
  
  sharedBP_PPIindex_list=list()
  for (i in 1:length(enrichedGOvec)){
    genesvec=unique(annogenes_mgroups$gene[which(annogenes_mgroups$go_id==enrichedGOvec[i])])
    ind=c(which(is.element(tissue_ppidatat_df$from,genesvec)),which(is.element(tissue_ppidatat_df$to,genesvec)))
    sharedind=ind[which(duplicated(ind))]
    sharedBP_PPIindex_mat[i,sharedind]=1
    sharedBP_PPIindex_list[[i]]=sharedind
  }  
  
  names(sharedBP_PPIindex_list)=unique(annogenes_mgroups$GOname)
  enrichedGOvec=unique(annogenes_mgroups$GOname)
  
  jacc=paste0("merged_Jaccard", gsub("\\.", "", as.character(h_jacc)) , " ")
  vec=paste0("originalFEA ", length(group), " ; " , jacc ,length(groups_vec))
  
  
  if (output_version!="little"){
    list(sharedBP_PPIindex_list=sharedBP_PPIindex_list,
         sharedBP_PPIindex_mat=sharedBP_PPIindex_mat,
         tissue_ppidatat_df=tissue_ppidatat_df, 
         enrichedGOvec=enrichedGOvec, 
         groupGO_geneass=annogenes_mgroups,
         vec)
  }
  
  if (output_version=="little"){
    list(sharedBP_PPIindex_list=sharedBP_PPIindex_list,
         sharedBP_PPIindex_mat=0,
         tissue_ppidatat_df=tissue_ppidatat_df, 
         enrichedGOvec=enrichedGOvec, 
         groupGO_geneass=0,
         0)
  }
  
  
}