import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import scanpy as sc

def read_perturbations_csv(adata_here,cell2guide_csv,guide2gene_csv=None,pref='',
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
    print('annotated cells:',len(cells_annotated),'or',percent_cells_annotated_in_adata,'percent of adata')
    
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
                
    #make obs for each perturbation count
    for guide in guides:
        adata_here.obs[pref+guide]=cell2guide_for_adata.loc[:,guide]*1.0
    if guide2gene_csv!=None:
        for gene in genes:
            adata_here.obs[pref+gene]=cell2gene_for_adata.loc[:,gene]*1.0
            
    #and one obs for guide, one for gene. one compact (says "multiple"), and one with the full name
    guide_sum=cell2guide_for_adata.sum(axis=1)
    if guide2gene_csv!=None:
        gene_sum=cell2gene_for_adata.sum(axis=1)
    guide_anno_long=[]
    guide_anno_compact=[]
    gene_anno_long=[]
    gene_anno_compact=[]
    for i in range(adata_here.n_obs):
        x=cell2guide_for_adata.iloc[i,]
        guides_cell=list(x.loc[x>0].index)
        guides_cell.sort()
        guides_cell_str=','.join(guides_cell)

        if len(guides_cell)==0:
            guide_anno_compact_here='unassigned'
            guide_anno_long_here='unassigned'
        if len(guides_cell)==1:
            guide_anno_compact_here=guides_cell_str
            guide_anno_long_here=guides_cell_str
        if len(guides_cell)>1:
            guide_anno_compact_here='multiple'
            guide_anno_long_here=guides_cell_str
        guide_anno_long.append(guide_anno_long_here)
        guide_anno_compact.append(guide_anno_compact_here)

        if guide2gene_csv!=None:
            y=cell2gene_for_adata.iloc[i,]
            genes_cell=list(y.loc[y>0].index)
            genes_cell.sort()
            genes_cell_str=','.join(genes_cell)

            if len(genes_cell)==0:
                gene_anno_compact_here='unassigned'
                gene_anno_long_here='unassigned'
            if len(genes_cell)==1:
                gene_anno_compact_here=genes_cell_str
                gene_anno_long_here=genes_cell_str
            if len(genes_cell)>1:
                gene_anno_compact_here='multiple'
                gene_anno_long_here=genes_cell_str
            gene_anno_long.append(gene_anno_long_here)
            gene_anno_compact.append(gene_anno_compact_here)
        
    adata_here.obs['guide']=guide_anno_long
    adata_here.obs['guide.compact']=guide_anno_compact

    if guide2gene_csv!=None:
        adata_here.obs['gene']=gene_anno_long
        adata_here.obs['gene.compact']=gene_anno_compact

    #make an entry for unassigned cells
    adata_here.obs['unassigned']=1.0*(adata_here.obs['guide']=='unassigned')

    #adata_here.obsm[pref+'cell2guide']=cell2guide_for_adata
    #adata_here.obsm[pref+'cell2gene']=cell2gene_for_adata
    if copy:
        return(adata_here)

def read_perturbations_csv_old(adata_here,cell2guide_csv,guide2gene_csv=None,pref='',
                           sep='\t',copy=False):
    
    """
    Read in which perturbations were present in which cell.                                                
                                                                                                              
    Args                                                                                                      
    ----                                                                                                      
    my_adata: AnnData                                                                                         
        adata                                                                                                 
    cell2guide_csv: str                                                                                       
        csv file where each line is a cell and each column in a guide. It has a 1 if the guide is present in the cell and a 0 otherwise.                                                                                   
    guide2gene_csv: str                                                                                       
        (optional) a csv file mapping which gene is targeted by each guide. This is useful for analyses aggregating across all guides of a gene.                                                                           
    pref: str                                                                                                 
        (default: "perturb"): (optional) a prefix to add to annotations. This is meant to allow the user to have potentially multiple modes of perturbations allowed.                                                      
    sep: str                                                                                                  
        (default: "\\t"): (optional) separator in the csv file.                                               
    copy: bool                                                                                                
        (default: False): (optional) whether to return a copy of the annotation data                          
                            
    Returns                                                                                                   
    -------                                                                                                   
    None                                                                                                      
        adds the to adata (or copy thereof) the following fields:                                             
                                                                                                              
            adata.obsm[pref+'.cell2guide']                                                                    
                                                                                                              
            adata.obsm[pref+'.cell2gene'] 
    """

    if copy: adata_here = adata_here.copy()
    
    #read cell2guide
    cell2guide=pd.read_csv(cell2guide_csv,sep=sep)
    cell2guide.index=cell2guide['cell']
        
    #check that the annotated cells have a good overlap with the adata
    cells_annotated=list(set(cell2guide['cell']))
    cells_annotated_in_adata=list(set(cells_annotated).intersection(set(adata_here.obs_names)))
    percent_cells_annotated_in_adata=100*round(len(cells_annotated_in_adata)/adata_here.n_obs,2)
    
    print('annotated cells:',len(cells_annotated))
    print('adata cells:',adata_here.n_obs,'|',percent_cells_annotated_in_adata,'percent annotated')
    if len(cells_annotated_in_adata)==0:
        print('======\nERROR: no cells in your adata were annotated. Check that the cell names match between adata and your cell2guide file')
        return(None)

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
                
            
    #and one obs for guide, one for gene. one compact (says "multiple"), and one with the full name
    guide_sum=cell2guide_for_adata.sum(axis=1)
    guide_anno_long=[]
    guide_anno_compact=[]
    unassigned_vector=[]
    for i in range(adata_here.n_obs):
        x=cell2guide_for_adata.iloc[i,]
        guides_cell=list(x.loc[x>0].index)
        guides_cell.sort()
        guides_cell_str=','.join(guides_cell)
        #print(guides_cell_str)
        if len(guides_cell)==0:
            guide_anno_compact_here='unassigned'
            guide_anno_long_here='unassigned'
        if len(guides_cell)==1:
            guide_anno_compact_here=guides_cell_str
            guide_anno_long_here=guides_cell_str
        if len(guides_cell)>1:
            guide_anno_compact_here='multiple'
            guide_anno_long_here=guides_cell_str
        guide_anno_long.append(guide_anno_long_here)
        guide_anno_compact.append(guide_anno_compact_here)
        if guide_anno_compact_here=='unassigned':
            unassigned_vector.append(1.0)
        else:
            unassigned_vector.append(0.0)

    adata_here.obs[pref+'guide']=guide_anno_long
    adata_here.obs[pref+'guide.compact']=guide_anno_compact
    
    
    if guide2gene_csv!=None:
        gene_sum=cell2gene_for_adata.sum(axis=1)
        gene_anno_long=[]
        gene_anno_compact=[]
        for i in range(adata_here.n_obs):
            x=cell2gene_for_adata.iloc[i,]
            genes_cell=list(x.loc[x>0].index)
            genes_cell.sort()
            genes_cell_str=','.join(genes_cell)
            #print(guides_cell_str)
            if len(genes_cell)==0:
                gene_anno_compact_here='unassigned'
                gene_anno_long_here='unassigned'
            if len(genes_cell)==1:
                gene_anno_compact_here=genes_cell_str
                gene_anno_long_here=genes_cell_str
            if len(genes_cell)>1:
                gene_anno_compact_here='multiple'
                gene_anno_long_here=genes_cell_str
            gene_anno_long.append(gene_anno_long_here)
            gene_anno_compact.append(gene_anno_compact_here)

        adata_here.obs[pref+'gene']=gene_anno_long
        adata_here.obs[pref+'gene.compact']=gene_anno_compact
    
    #make obs for each perturbation count                                                                                                                                
    for guide in guides:
        adata_here.obs[pref+guide]=1.0*(cell2guide_for_adata.loc[:,guide]>0)
    #also make one for unassigned
    adata_here.obs['unassigned']=unassigned_vector
    if guide2gene_csv!=None:
        for gene in genes:
            adata_here.obs[pref+'gene.'+gene]=1.0*(cell2gene_for_adata.loc[:,gene]>0)


    adata_here.obsm[pref+'cell2guide']=cell2guide_for_adata
    if guide2gene_csv!=None:
        adata_here.obsm[pref+'cell2gene']=cell2gene_for_adata
    if copy:
        return(adata_here)


def compute_TPT(gbcs_dataset):

    '''                                                                                                                                                                              
    input: pandas data frame with the columns "cbc", "umi", "gbc", "r2" where every row is a read                                                                                    
    output: pandas data frame with the columns "gbc", "cbc", "umi", "cbc-umi-r2-count", "cbc-umi-count", "TPT"                                                                       
    NOTE: for the input, multiple reads corresponding to the same cbc-umi combination should be listed as separate lines!                                                            
    '''

    import copy
    import re
    import time
    import pandas as pd
    import numpy as np

    print("======== annotating cbc-umi pairs, and cbc-umi-r2")
    gbcs_dataset['cbcumi']=[x+'-'+y for x,y in zip(gbcs_dataset['cbc'],gbcs_dataset['umi'])]
    cbcumir2=list([x+'+'+y for x,y in zip(gbcs_dataset['cbcumi'],gbcs_dataset['r2'])])
    gbcs_dataset['cbcumir2']=list([x+'_gbc_'+y for x,y in zip(cbcumir2,gbcs_dataset['gbc'])])

    print("======== counting the numbers of reads supporting cbc-umi-r2 and for denominator cbc-umi")
    cbcumi_group=gbcs_dataset.groupby('cbcumi').size()
    cbcumi_r2_group=gbcs_dataset.groupby('cbcumir2').size()
    cbcumi_from_grouped_reads=[x.split('+')[0] for x in cbcumi_r2_group.index]

    print("======== computing TPT")
    #divide every cbc-umi-r2 value by the cbc-umi values                                                                                                                             
    combo_counts=pd.DataFrame({'cbc-umi-r2':cbcumi_r2_group,
                              'cbc-umi-r2-name':cbcumi_r2_group.index,
                              'cbc-umi':cbcumi_from_grouped_reads})
    combo_counts['cbc-umi-total']=copy.deepcopy(list(cbcumi_group.loc[combo_counts['cbc-umi']]))
    combo_counts['TPT']=1.0*combo_counts['cbc-umi-r2']/combo_counts['cbc-umi-total']

    combo_counts['gbc']=list([x.split('_gbc_')[1] for x in list(combo_counts.loc[:,'cbc-umi-r2-name'])])
    combo_counts['umi']=list([x.split('-')[1] for x in list(combo_counts['cbc-umi'])])
    combo_counts['cbc']=list([x.split('-')[0] for x in list(combo_counts['cbc-umi'])])

    print("======== compiling the final result")
    to_return=pd.DataFrame({'gbc':combo_counts['gbc'],
                           'cbc':combo_counts['cbc'],
                           'umi':combo_counts['umi'],
                           'cbc-umi-r2-count':combo_counts['cbc-umi-r2'],
                           'cbc-umi-count':combo_counts['cbc-umi-total'],
                           'TPT':combo_counts['TPT']})
    to_return=to_return.reset_index(drop=True)
    return(to_return)
