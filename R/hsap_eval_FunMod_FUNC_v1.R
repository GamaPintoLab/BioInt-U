
# 3: Function to CHARACTERIZE BioInt Libraries
# Recover additional properties of the BioInt units 

# Function to run ONLY ONCE FOR THE COMPLETE LIST output from hsap_FEA_FunMod_FUNC_v1 function

# @param funmod_alltissues_list0 Output list from hsap_FEA_FunMod_FUNC_v1 function
# @param funmod_alltissues_list  Output list from hsap_processing_FunMod_FUNC_v1 function 
# @param filtered Character string to characterize all BioInt unit (filtered !="yes") or just the selected in hsap_processing_FunMod_FUNC_v1 function (filtered =="yes")

# @return List including tissue_funmod_summarydf, FUNMOD_PROP, REDUNDANCY_df, MAINCOMP_df and SHAPE_df
# @return tissue_funmod_summarydf Data frame for summary of properties of BioInt units by tissue
# @return FUNMOD_PROP Data frame including all properties of BioInt units in all the tissues evaluated
# @return REDUNDANCY_df Data frame including additional properties of BioInt units
# @return MAINCOMP_df Data frame including additional properties of BioInt units
# @return SHAPE_df Data frame including additional properties of BioInt units


library("igraph")



hsap_eval_FunMod_FUNC_v1=function(funmod_alltissues_list0, funmod_alltissues_list, filtered="yes"){
  
  OP=3 
  if (filtered!="yes"){
    OP=2
  }
  
  tissue_funmod_summarydf=data.frame("tissue"=names(unlist(lapply(lapply(lapply(funmod_alltissues_list, "[[", OP), "[",2), nrow))),
                                     "n_funmod"=unlist(lapply(lapply(lapply(funmod_alltissues_list, "[[", OP), "[",2), nrow)),
                                     "n_uniquegenes"=0, "n_ubiqgenes"=0,"perc_ubiq"=0,
                                     "min_size"=0,  "1stQ_size"=0, "median_size"=0, "mean_size"=0, "3rdQ_size"=0,"max_size"=0)
  
  for (i in 1:length(funmod_alltissues_list)){
    tissue_genes=unique(unlist(strsplit(unlist(funmod_alltissues_list[[i]][[OP]][7]), split="/")))
    tissue_funmod_summarydf[i,3]=length(tissue_genes) # -after 90% main comp filtering
    tissue_funmod_summarydf[i,4]=length(intersect(ubiq, tissue_genes))
    tissue_funmod_summarydf[i,5]=(length(intersect(ubiq, tissue_genes))*100)/length(tissue_genes)
    tissue_funmod_summarydf[i,6:11]=summary(unlist(funmod_alltissues_list[[i]][[OP]][9])) ## maincompsize size-after 90% main comp filtering summary
  }
  
  
  
  #####################
  
  # overlap: 
  # uniqueness: % of unique proteins on funmod
  REDUNDANCY_df=data.frame()
  for (i in 1:length(funmod_alltissues_list)){
    genefreq=table((unlist(funmod_alltissues_list[[i]]$genes_by_funct_list)))
    uniquegenes=names(which(genefreq==1))
    subREDUNDANCY_df=data.frame("tissue"=names(funmod_alltissues_list)[i], 
                                "funmod"=funmod_alltissues_list[[i]][[OP]]$goname,
                                "uniqueness"=as.vector(unlist(lapply(funmod_alltissues_list[[i]]$genes_by_funct_list, function(x,y) (length(which(is.element(uniquegenes,x)))*100)/length(x)))),
                                "overlap"=as.vector(unlist(lapply(funmod_alltissues_list[[i]]$funmod_overlapping_list, function(x) length(x)))),
                                "nsub"=unlist(funmod_alltissues_list[[i]][[OP]]$numcomp),
                                "finalsize"=unlist(lapply(unlist(funmod_alltissues_list[[i]][[OP]]$genes), function(x) length(unique(unlist(strsplit(x, split="/")))))),
                                stringsAsFactors = F)
    selFM=funmod_alltissues_list[[i]][[OP]]$goname
    # subMAINCOMP_df$sel[which(is.element(as.character(subMAINCOMP_df$funmod),selFM))]=1
    REDUNDANCY_df=rbind(REDUNDANCY_df,subREDUNDANCY_df)
  }
  
  
  
  MAINCOMP_df=data.frame()
  for (i in 1:length(funmod_alltissues_list)){
    subMAINCOMP_df=data.frame("tissue"=names(funmod_alltissues_list)[i], 
                              "funmod"=funmod_alltissues_list[[i]][[OP]]$goname,
                              "maincomp_perc"=unlist(funmod_alltissues_list[[i]][[OP]]$maincomp),
                              "nsub"=unlist(funmod_alltissues_list[[i]][[OP]]$numcomp),
                              "finalsize"=unlist(lapply(unlist(funmod_alltissues_list[[i]][[OP]]$genes), function(x) length(unique(unlist(strsplit(x, split="/")))))),
                              "sel"=0, stringsAsFactors = F)
    
    selFM=funmod_alltissues_list[[i]][[OP]]$goname
    # subMAINCOMP_df$sel[which(is.element(as.character(subMAINCOMP_df$funmod),selFM))]=1
    MAINCOMP_df=rbind(MAINCOMP_df,subMAINCOMP_df)
  }
  
  
  SHAPE_df=data.frame()
  for (i in 1:length(funmod_alltissues_list)){
    hsgraphALL=igraph::graph_from_data_frame(funmod_alltissues_list0[[i]]$tissue_ppidatat_df,directed=F)
    hsgraphALL=igraph::simplify(hsgraphALL)
    bet_graph=igraph::betweenness(hsgraphALL)
    deg_graph=igraph::degree(hsgraphALL)
    clust_graph=igraph::transitivity(hsgraphALL, type="local", isolates="zero")
    subSHAPE_df=data.frame("tissue"=names(funmod_alltissues_list)[i], 
                           "funmod"=funmod_alltissues_list[[i]][[OP]]$goname,
                           "nsub_sel"=0, 
                           "diam"=0, 
                           # "mean_transitivity"=0,
                           "sd_ecc"=0,
                           
                           "between"=0,
                           "degree"=0,
                           "clustering"=0, 
                           "finalsize"=unlist(lapply(unlist(funmod_alltissues_list[[i]][[OP]]$genes), function(x) length(unique(unlist(strsplit(x, split="/")))))),
                           stringsAsFactors = F)

    for ( j in 1:nrow(subSHAPE_df)){
      geness=unlist(strsplit(funmod_alltissues_list[[i]][[OP]]$genes[j], split="/"))
      indexes=which(is.element(igraph::vertex_attr(hsgraphALL)$name, geness))
      funmod=igraph::induced_subgraph(hsgraphALL, indexes)
      subSHAPE_df$nsub_sel[j]=igraph::components(funmod)$no
      subSHAPE_df$diam[j]=igraph::diameter(funmod)
      #subSHAPE_df$mean_transitivity[j]=1/mean(transitivity(funmod, type="local", isolates="zero")) # including NAs as zero to the mean calculation
      subSHAPE_df$sd_ecc[j]=sd(igraph::eccentricity(funmod))
      
      subSHAPE_df$between[j]=mean(bet_graph[indexes])
      subSHAPE_df$degree[j]=mean(deg_graph[indexes])
      subSHAPE_df$clustering[j]=mean(clust_graph[indexes])
      
    }
    SHAPE_df=rbind(SHAPE_df,subSHAPE_df)
  }
  
  
  FUNMOD_PROP=cbind(REDUNDANCY_df[,-7], "maincomp_perc"=MAINCOMP_df$maincomp_perc, SHAPE_df[,c(4:8)])
  
  list(tissue_funmod_summarydf=tissue_funmod_summarydf,
       FUNMOD_PROP=FUNMOD_PROP, 
       REDUNDANCY_df=REDUNDANCY_df[,-7], 
       MAINCOMP_df=MAINCOMP_df,
       SHAPE_df=SHAPE_df)
  
}
