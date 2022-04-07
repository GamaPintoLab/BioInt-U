


# Function to PROCESS BioInt Libraries

# 2: Characterize the topological properties of BioInt unit 
# Filter BioInt unit according to user-defined threshold

# Function to run ITERATIVELY FOR EACH OBJECT in the list of hsap_FEA_FunMod_FUNC_v1 output list

# @param FunMod_output Output list from hsap_FEA_FunMod_FUNC_v1 function
# @param ubiq_genes_vec  Vector of genes expressed in all tissues (ubiquitous/housekeeping genes)
# @param maincomp_thresh Numeric value. Threshold % of proteins included in the main component of the BioInt unit to define decide whether to include subnetworks aside the main component
# @param min_size Numeric value to define lower threshold size (minmum number of proteins by BioInt unit)
# @param max_size Numeric value to define upper threshold size (maxium number of proteins by BioInt unit)

# @return List including genes_by_funct_list, componenet_subg_infodf, componenet_subg_infodf_filtered, jaccard_uppertri and funmod_overlapping_list
# @return genes_by_funct_list List of genes by BioInt unit
# @return componenet_subg_infodf Data frame including all properties of all BioInt unit (not filtered)
# @return componenet_subg_infodf_filtered  Data frame including all properties of selected BioInt unit (filtered by size)
# @return jaccard_uppertri Matrix (upper triangle) of Jaccard analysis for selected BioInt unit 
# @return funmod_overlapping_list List format of jaccard_uppertri matrix

library("igraph")
library("philentropy")


hsap_processing_FunMod_FUNC_v1=function(FunMod_output,ubiq_genes_vec, maincomp_thresh=90, min_size=10, max_size=300){
  
  componenet_subg_infodf=data.frame("goid"=FunMod_output[[4]], "goname"=names(FunMod_output[[1]]), "numprot"=0, "numcomp"=0, "csize"=0, "maincomp"=0, 
                                    #"meandist"=0,
                                    stringsAsFactors = F)
  
  for ( i in 1:length(FunMod_output[[1]])){
    subg=igraph::graph_from_data_frame(FunMod_output[[3]][FunMod_output[[1]][[i]],], directed = F)
    componenet_subg_infodf$numprot[i]=igraph::vcount(subg)
    componenet_subg_infodf$numcomp[i]=igraph::components(subg)$no
    ordered_csizes=igraph::components(subg)$csize[order(igraph::components(subg)$csize, decreasing = T)]
    componenet_subg_infodf$csize[i]=paste(ordered_csizes, sep=",", collapse = ",")
    componenet_subg_infodf$maincomp[i]=ordered_csizes[1]*100/igraph::vcount(subg)
    # componenet_subg_infodf$meandist[i]=mean_distance(subg, directed = F)
  }
  
  test_componenet_subg_infodf=cbind(componenet_subg_infodf, "genes"=0,"symbol"=0)
  for (i in 1:nrow(test_componenet_subg_infodf)){
    subg=igraph::graph_from_data_frame(FunMod_output[[3]][FunMod_output[[1]][[i]],], directed = F)
    if (test_componenet_subg_infodf$maincomp[i]>=maincomp_thresh){
      subg_comp=igraph::components(subg)
      maincomp=order(subg_comp$csize, decreasing = T)[1]
      geness=names(which(subg_comp$membership==maincomp))
      test_componenet_subg_infodf$genes[i]=paste(geness, sep="/", collapse="/")
      test_componenet_subg_infodf$symbol[i]=test_componenet_subg_infodf$genes[i]
    }  
    if (test_componenet_subg_infodf$maincomp[i]<maincomp_thresh){
      geness=igraph::vertex_attr(subg)$name
      test_componenet_subg_infodf$genes[i]=paste(geness, sep="/", collapse="/")
      test_componenet_subg_infodf$symbol[i]=test_componenet_subg_infodf$genes[i]
    }  
  }
  componenet_subg_infodf=test_componenet_subg_infodf
  componenet_subg_infodf_filtered=subset(test_componenet_subg_infodf, numprot>=min_size & numprot<=max_size)
  
  
  templist=list() # list of genes by function
  for ( i in 1:nrow(componenet_subg_infodf_filtered)){
    templist[[i]]=unlist(strsplit(componenet_subg_infodf_filtered$genes[i], split="/"))
  }
  names(templist)=componenet_subg_infodf_filtered$goname
  tempdf=data.frame(table(unlist(templist)),stringsAsFactors = F)
  
  allgenes=unique(unlist(templist))
  genefunmod_matrix=matrix(0, ncol=length(allgenes), nrow=length(templist))
  colnames(genefunmod_matrix)=allgenes
  rownames(genefunmod_matrix)=names(templist)
  
  for ( i in 1:length(templist)){
    inndex=which(is.element(colnames(genefunmod_matrix), templist[[i]]))
    genefunmod_matrix[i,inndex]=1
  }
  
  # which fun mod overlap with other funmods
  if (nrow(genefunmod_matrix)>1){
    funmod_overlapping_list=list()
    for ( i in 1:nrow(genefunmod_matrix)){
      submat=genefunmod_matrix[-i,which(genefunmod_matrix[i,]==1)]
      funmod_overlapping_list[[i]]=names(which(apply(submat, 1, function(x) sum(x))>=1))
    }
    names(funmod_overlapping_list)=names(templist)
    
    
    componenet_subg_infodf_filtered2=cbind(componenet_subg_infodf_filtered,"maincomp_size"=0, "overlap"=0, "uniquenness"=0, "ubiq_perc"=0)
    for ( i in 1:nrow(componenet_subg_infodf_filtered)){
      thefb=templist[[which(names(templist) ==  componenet_subg_infodf_filtered$goname[i])]]
      componenet_subg_infodf_filtered2$overlap[i]=summary(tempdf$Freq[which(tempdf$Var1 %in% thefb)])[4]
      componenet_subg_infodf_filtered2$uniquenness[i]=
        (length(which(tempdf$Freq[which(tempdf$Var1 %in% thefb)]==1))*100)/length(thefb)
      componenet_subg_infodf_filtered2$ubiq_perc[i]=
        (length(which(is.element(thefb,ubiq_genes_vec)))*100)/length(thefb)
      componenet_subg_infodf_filtered2$maincomp_size[i]=length(thefb)
      
    }
    
    for (i in 1:nrow(componenet_subg_infodf_filtered2)){
      ind=which(names(funmod_overlapping_list)==componenet_subg_infodf_filtered2$goname[i])
      componenet_subg_infodf_filtered2$overlap[i]=length(funmod_overlapping_list[[ind]])
    }
    
    geneinmods=unique(unlist(templist))
    funmod_overlapping_mat=matrix(0,ncol=length(funmod_overlapping_list), nrow=length(geneinmods))
    colnames(funmod_overlapping_mat)=names(funmod_overlapping_list)
    rownames(funmod_overlapping_mat)=(geneinmods)
    
    for ( j in 1:ncol(funmod_overlapping_mat)){
      ind=which(rownames(funmod_overlapping_mat) %in% templist[[j]])
      funmod_overlapping_mat[ind,j]=1
    }
    
    all_genes_mat_dist=philentropy::distance(t(funmod_overlapping_mat), method="jaccard")
    colnames(all_genes_mat_dist)=colnames(funmod_overlapping_mat)
    rownames(all_genes_mat_dist)=colnames(funmod_overlapping_mat)
    
    get_upper_tri <- function(all_genes_mat_dist){
      all_genes_mat_dist[lower.tri(all_genes_mat_dist)]<- NA
      return(all_genes_mat_dist)
    }
    
    testuppertri=get_upper_tri(all_genes_mat_dist)
    colnames(testuppertri)=colnames(all_genes_mat_dist)
    rownames(testuppertri)=rownames(all_genes_mat_dist)
    
    jaccard_uppertri=testuppertri
    componenet_subg_infodf_filtered=componenet_subg_infodf_filtered2
    genes_by_funct_list=templist
    
  }
  
  if (nrow(genefunmod_matrix)<=1){
    jaccard_uppertri=0
    funmod_overlapping_list=0
    all_genes_mat_dist=0
  }
  
  list(genes_by_funct_list=genes_by_funct_list,  
       componenet_subg_infodf=componenet_subg_infodf, 
       componenet_subg_infodf_filtered=componenet_subg_infodf_filtered, 
       jaccard_uppertri=jaccard_uppertri,
       funmod_overlapping_list=funmod_overlapping_list) 

}
