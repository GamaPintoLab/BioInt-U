# **BioInt-U method**



The goal of BioInt-U method is to identify *Biological Interacting* units. 



First, it reconstructs tissue-specific protein-protein interaction networks. Each tissue network is functionally enriched in Gene Ontology Biological Processes. The GO-BP terms are used to identify **BioInt units**, defined as groups of proteins physically interacting and sharing the same enriched terms. GO-BP sharing an excesive fraction of proteins can be collapsed or simplified using  Jaccard coefficient.



Once BioInt libraries are generated, we also make available the code to functionally and topologically characterize the BioInt units whitin their tissue-specific interactomes.



The R functions can be adapted to other tissue expression and protein-protein interaction data from species available in GOFuncR Bioconductor package.
