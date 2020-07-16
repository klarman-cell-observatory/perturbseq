import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

def read_perturbations_csv(my_adata,
                           cell2guide_csv,
                           guide2gene_csv = None,
                           pref = 'perturb',
                           sep = '\t',
                           copy = False):

    """
    Read in which perturbations were present in which cell.

    Args:
    --------

    my_adata
        adata
    cell2guide_csv 
        csv file where each line is a cell and each column in a guide. It has a 1 if the guide is present in the cell and a 0 otherwise.
    guide2gene_csv
        (optional) a csv file mapping which gene is targeted by each guide. This is useful for analyses aggregating across all guides of a gene.
    pref (default: "perturb")
        (optional) a prefix to add to annotations. This is meant to allow the user to have potentially multiple modes of perturbations allowed.
    sep (default: "\\t")
        (optional) separator in the csv file. 
    copy (default: False)
        (optional) whether to return a copy of the annotation data
    
    
    Returns:
    -----------
    
    Adds the to adata (or copy thereof) the following fields:
        adata.obsm[pref+'.cell2guide']
        adata.obsm[pref+'.cell2gene'] 
        
    """

    import pandas as pd

    if copy: my_adata = my_adata.copy()
    
    #read cell2guide
    cell2guide=pd.read_csv(cell2guide_csv,sep=sep)
    cell2guide.index=cell2guide['cell']
        
    #check that the annotated cells have a good overlap with the adata
    cells_annotated=list(set(cell2guide['cell']))
    cells_annotated_in_adata=list(set(cells_annotated).intersection(set(my_adata.obs_names)))
    percent_cells_annotated_in_adata=100*round(len(cells_annotated_in_adata)/my_adata.n_obs,2)
    print('adata cells:',my_adata.n_obs)
    print('annotated cells:',len(cells_annotated),'of which',percent_cells_annotated_in_adata,'percent in adata')
    
    #assign the guides to the cells in adata
    guides=list(set(cell2guide.columns).difference(set(['cell'])))
    cell2guide_for_adata=pd.DataFrame(0.0,
                      index=my_adata.obs_names,
                      columns=guides)
    cell2guide_for_adata.loc[cells_annotated_in_adata,guides]=cell2guide.loc[cells_annotated_in_adata,guides]
    
    #if provided a guide->gene file, save that in the adata object, as another obsm
    if guide2gene_csv!=None:
        guide2gene=pd.read_csv(guide2gene_csv,sep=sep)
        genes=list(set(guide2gene['gene']))
        cell2gene_for_adata=pd.DataFrame(0.0,
                      index=my_adata.obs_names,
                      columns=genes)
        #for each gene, check if at least one of the guides is present
        for gene in genes:
            guides_per_gene_init=list(guide2gene.loc[guide2gene['gene']==gene,'guide'])
            guides_per_gene=list(set(guides_per_gene_init).intersection(set(guides)))
            if len(guides_per_gene)>0:
                cell2gene_for_adata.loc[cells_annotated_in_adata,gene]=cell2guide_for_adata.loc[cells_annotated_in_adata,guides_per_gene].sum(axis=1)
    
    my_adata.obsm[pref+'.cell2guide']=cell2guide_for_adata
    my_adata.obsm[pref+'.cell2gene']=cell2gene_for_adata
    if copy:
        return(my_adata)


