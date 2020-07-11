import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

def read_perturbations_csv(adata_here,cell2guide_csv,guide2gene_csv=None,pref='perturb',
                           sep='\t',copy=False):
    
    import pandas as pd

    if copy: adata_here = adata_here.copy()
    
    #read cell2guide
    cell2guide=pd.read_csv(cell2guide_csv,sep=sep)
    cell2guide.index=cell2guide['cell']
        
    #check that the annotated cells have a good overlap with the adata
    cells_annotated=list(set(cell2guide['cell']))
    cells_annotated_in_adata=list(set(cells_annotated).intersection(set(adata_here.obs_names)))
    percent_cells_annotated_in_adata=100*round(len(cells_annotated_in_adata)/adata_here.n_obs,2)
    print('adata cells:',adata_here.n_obs)
    print('annotated cells:',len(cells_annotated),'of which',percent_cells_annotated_in_adata,'percent in adata')
    
    #assign the guides to the cells in adata
    guides=list(set(cell2guide.columns).difference(set(['cell'])))
    cell2guide_for_adata=pd.DataFrame(0.0,
                      index=adata_here.obs_names,
                      columns=guides)
    cell2guide_for_adata.loc[cells_annotated_in_adata,guides]=cell2guide.loc[cells_annotated_in_adata,guides]
    
    #if provided a guide->gene file, save that in the adata object, as another obsm
    if guide2gene_csv!=None:
        guide2gene=pd.read_csv(guide2gene_csv,sep=sep)
        genes=list(set(guide2gene['gene']))
        cell2gene_for_adata=pd.DataFrame(0.0,
                      index=adata_here.obs_names,
                      columns=genes)
        #for each gene, check if at least one of the guides is present
        for gene in genes:
            guides_per_gene_init=list(guide2gene.loc[guide2gene['gene']==gene,'guide'])
            guides_per_gene=list(set(guides_per_gene_init).intersection(set(guides)))
            if len(guides_per_gene)>0:
                cell2gene_for_adata.loc[cells_annotated_in_adata,gene]=cell2guide_for_adata.loc[cells_annotated_in_adata,guides_per_gene].sum(axis=1)
    
    adata_here.obsm[pref+'.cell2guide']=cell2guide_for_adata
    adata_here.obsm[pref+'.cell2gene']=cell2gene_for_adata
    if copy:
        return(adata_here)


