import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns

#========= may have to split some of these methods into utils later
def get_cls_adata(adata_here,n_neighbors):
    sc.pp.neighbors(adata_here,n_neighbors=n_neighbors,use_rep='X')
    sc.tl.louvain(adata_here)
    sc.tl.umap(adata_here)
    sc.pl.umap(adata_here,color=['louvain'])
    return(adata_here.obs['louvain'])

#a method to automatically order the groups
def average_df_by_group(df,grouping):
    #groups on the rows
    
    assert df.shape[0]==len(grouping)
    
    groups=list(set(grouping))
    compact_df=pd.DataFrame(0,index=groups,columns=df.columns)
    for group in groups:
        rows_here=pd.DataFrame(df.loc[grouping==group,:])
        
        if rows_here.shape[0]==1:
            compact_df.loc[group,:]=np.array(rows_here).flatten()
        else:
            compact_df.loc[group,:]=rows_here.mean(axis=0)
    return(compact_df)

def get_cl_order(df,cl):
    #clusters are on the rows
    cluster_level_df=average_df_by_group(df,cl)
    cluster_level_df.index=list(cluster_level_df.index)
    
    #cluster in some way
    from scipy.cluster import hierarchy
    from scipy.spatial.distance import pdist
    np.random.seed(seed=7)
    cl_Z = hierarchy.linkage(cluster_level_df,
                         optimal_ordering=True,
                         method='average')
    cl_dn = hierarchy.dendrogram(cl_Z, labels=cluster_level_df.index)

    ax = plt.gca()
    labels = ax.get_xmajorticklabels()
    ordered_cls=[]
    for i in range(len(labels)):
        ordered_cls.append(labels[i].get_text())
    plt.show()
    
    #return the order of the rows of cluster_level_df
    return(ordered_cls)

def sort_clustered_df(df,cl_rows=None,cl_cols=None,cluster_within=True):
    row_order=df.index
    col_order=df.columns
    if cl_rows is not None:
        ordered_row_cls=get_cl_order(df,cl_rows)
        row_order=[]
        for r_idx in range(len(ordered_row_cls)):
            row_cl=ordered_row_cls[r_idx]
            rows_cluster_here=[]
            for rowname_idx in range(df.shape[0]):
                rowname=df.index[rowname_idx]
                if cl_rows[rowname_idx]==row_cl:
                    rows_cluster_here.append(rowname)
            if cluster_within:
                rows_cluster_here=get_cl_order(df.loc[rows_cluster_here,:],pd.Series(rows_cluster_here,
                                                                                    index=rows_cluster_here))
            for i in range(len(rows_cluster_here)):
                row_order.append(rows_cluster_here[i])
            

    if cl_cols is not None:
        ordered_col_cls=get_cl_order(df.T,cl_cols)
        col_order=[]
        for c_idx in range(len(ordered_col_cls)):
            col_cl=ordered_col_cls[c_idx]
            cols_cluster_here=[]
            for colname_idx in range(df.shape[1]):
                colname=df.columns[colname_idx]
                if cl_cols[colname_idx]==col_cl:
                    cols_cluster_here.append(colname)
            if cluster_within:
                cols_cluster_here=get_cl_order(df.T.loc[cols_cluster_here,:],pd.Series(cols_cluster_here,
                                                                                    index=cols_cluster_here))
            for i in range(len(cols_cluster_here)):
                col_order.append(cols_cluster_here[i])
    return(df.loc[row_order,col_order])

def rename_by_order(ordered_df,label):
    
    old2new={}
    program_counter=0
    for i in range(ordered_df.shape[0]):
        cl=label[list(ordered_df.index)[i]]
        if cl not in old2new:
            old2new[cl]=program_counter
            program_counter+=1
    new_label=[]
    for i in range(ordered_df.shape[0]):
        cl=label[list(ordered_df.index)[i]]
        new_label.append(old2new[cl])
    new_label_df=pd.Series(new_label,index=ordered_df.index)
    return(new_label_df)

def cat2color(category_vector,color_map=None,cmap='Set2',**kwargs):
    
    # from https://stackoverflow.com/questions/26139423/plot-different-color-for-different-categorical-levels-using-matplotlib
    if color_map==None:
        color_labels = category_vector.unique()

        # List of RGB triplets
        rgb_values = sns.color_palette(cmap, len(color_labels))

        # Map label to RGB
        color_map = dict(zip(color_labels, rgb_values))

    color_vector=category_vector.map(color_map)
    return(color_vector)

def perturbation_modules(adata_here,input_type='bulk',perturbation_name='guide.compact',
                         n_neighbors=10,cluster_within=True,copy=False,cmap='Set2'):
    
    if copy: adata_here = adata_here.copy()
    
    #make perturbations x gene adata
    import copy
    sc_bulk=copy.deepcopy(adata_here.uns[input_type+'.'+perturbation_name])
    sc_bulk=sc_bulk.sort_index(axis=0)
    sc_bulk=sc_bulk.sort_index(axis=1)
    bulk_adata_perturbations = sc.AnnData(sc_bulk)
    
    #get the clustering
    perturbations_cl=get_cls_adata(bulk_adata_perturbations,
                           n_neighbors=n_neighbors)
    
    #sort the resulting clusters and perturbations within
    ordered_df=sort_clustered_df(sc_bulk,
                             perturbations_cl,
                             None,
                             cluster_within=cluster_within)
    
    #get the new order of modules and rename them, and give them colors
    perturbations_cl_sorted=rename_by_order(ordered_df,perturbations_cl)
    perturbations_cl_color=cat2color(perturbations_cl_sorted,cmap=cmap)
    perturbations_cl_sorted=pd.DataFrame({'module':perturbations_cl_sorted},
                                            index=perturbations_cl_sorted.index)
    perturbations_cl_color=pd.DataFrame({'color':perturbations_cl_color},
                                            index=perturbations_cl_color.index)
    
    adata_here.uns[input_type+'.'+perturbation_name+'.perturbation_module']=perturbations_cl_sorted
    adata_here.uns[input_type+'.'+perturbation_name+'.perturbation_module_color']=perturbations_cl_color
    adata_here.uns[input_type+'.'+perturbation_name]=ordered_df
    
    if copy:
        return(adata_here)

def gene_programs(adata_here,input_type='bulk',perturbation_name='guide.compact',
                  n_neighbors=5,cluster_within=True,copy=False,cmap='Set2'):
    
    if copy: adata_here = adata_here.copy()
    
    #make gene x pert adata
    import copy
    sc_bulk=copy.deepcopy(adata_here.uns[input_type+'.'+perturbation_name])
    sc_bulk=sc_bulk.sort_index(axis=0)
    sc_bulk=sc_bulk.sort_index(axis=1)
    bulk_adata_genes = sc.AnnData(sc_bulk.T)
    
    #get the clustering
    genes_cl=get_cls_adata(bulk_adata_genes,
                           n_neighbors=n_neighbors)
    
    #sort the resulting clusters and genes within
    ordered_df=sort_clustered_df(sc_bulk.T,
                             genes_cl,
                             None,
                             cluster_within=cluster_within)
    
    #get the new order of programs and rename them, and give them colors
    genes_cl_sorted=rename_by_order(ordered_df,genes_cl)
    genes_cl_color=cat2color(genes_cl_sorted,cmap=cmap)
    genes_cl_sorted=pd.DataFrame({'program':genes_cl_sorted},
                                            index=genes_cl_sorted.index)
    genes_cl_color=pd.DataFrame({'color':genes_cl_color},
                                            index=genes_cl_color.index)

    adata_here.var[input_type+'.'+perturbation_name+'.program']=genes_cl_sorted.loc[adata_here.var_names]
    adata_here.var[input_type+'.'+perturbation_name+'.program_color']=genes_cl_color.loc[adata_here.var_names]
    adata_here.uns[input_type+'.'+perturbation_name]=ordered_df.T
    
    if copy:
        return(adata_here)
    
def gene_programs_and_perturbation_modules(adata_here,input_type='bulk',perturbation_name='guide.compact',
                                           n_neighbors_programs=5,n_neighbors_modules=5,cluster_within=True,copy=False,cmap_programs='Set2',cmap_modules='Set2'):
    
    if copy: adata_here = adata_here.copy()
        
    gene_programs(adata_here,input_type,perturbation_name,
             n_neighbors=n_neighbors_programs,
                  cluster_within=cluster_within,copy=copy,cmap=cmap_programs)
    perturbation_modules(adata_here,input_type,perturbation_name,
             n_neighbors=n_neighbors_modules,
                         cluster_within=cluster_within,copy=copy,cmap=cmap_modules)
    
    if copy:
        return(adata_here)

#==================================================================

def bulk(adata_here,grouping_variable,by_batch=False,return_matrix=False):
    
    """Compute an in silico bulk set of expression profiles, based on cell labels
 
    Parameters
    ----------
    adata_here : `scanpy Anndata`
    grouping_variable : `str`
        The name of the variable that specifies a label for each cell. This variable must be accessible as `adata_here.obs[grouping_variable]`
    by_batch : `bool`
        Whether to combine data from cells with the same label but from different batches.
        If this is set to True, adata_here must have a adata_here.obs["batch"]
    
    Returns
    -------
    profile_matrix_df : a pandas DataFrame of size (number of conditions) x (number of genes). 
                        The number of conditions is: number of unique labels in `adata_here.obs[grouping_variable]` if by_batch==False
                                                     number of unique labels times the number of batches if by batch==True
    """
    
    #construct the profiles
    profiles=list(set(adata_here.obs[grouping_variable]))
    adata_here.obs['profile']=adata_here.obs[grouping_variable]
    
    if by_batch:
        profile_list=[]
        #make a new variable that combines batch and variable into 1
        for cell_idx in range(len(adata_here.obs_names)):
            profile=adata_here.obs['batch'][cell_idx]+'_'+adata_here.obs[grouping_variable][cell_idx]
            if profile not in profile_list:
                profile_list.append(profile)
        adata_here.obs['profile']=profile_list
        profiles=profile_list
        
    genes=adata_here.var_names
    profile_matrix=np.zeros((len(profiles),len(genes)))
    for profile_idx in range(len(profiles)):
        profile=profiles[profile_idx]
        cells_with_profile=list(adata_here.obs_names[adata_here.obs['profile']==profile])
        data_profile=adata_here[cells_with_profile,:].X.toarray()
        profile_matrix[profile_idx,:]=data_profile.mean(axis=0)
    profile_matrix_df=pd.DataFrame(profile_matrix)
    profile_matrix_df.index=profiles
    profile_matrix_df.columns=genes
    adata_here.uns['bulk.'+grouping_variable]=profile_matrix_df
    if return_matrix:
        return(profile_matrix_df)

def enriched_features(adata,f1='louvain',f2='batch',fdr=0.05,
                      copy=False,add_min_pval=True,pval_correction='fdr_bh',ps=1e-10):
    
    if copy: adata=adata.copy()
    
    import scipy

    f1s=list(set(adata.obs[f1]))
    f2s=list(set(adata.obs[f2]))

    oddsratios=np.zeros((len(f1s),len(f2s)))
    pvals=np.zeros((len(f1s),len(f2s)))
    proportions=np.zeros((len(f1s),len(f2s)))

    for f1_idx in range(len(f1s)):
        f1_here=f1s[f1_idx]
        
        cells_in_f1=list(adata.obs_names[adata.obs[f1]==f1_here])
        
        for f2_idx in range(len(f2s)):
            f2_here=f2s[f2_idx]
            cells_in_f2=list(adata.obs_names[adata.obs[f2]==f2_here])

            total=list(adata.obs_names)
            overlap=list(set(cells_in_f1).intersection(set(cells_in_f2)))
            contingency_table=np.array([[len(overlap),
                                         len(cells_in_f1)-len(overlap)],
                                        [len(cells_in_f2)-len(overlap),
                                         0]])
            contingency_table[1,1]=len(total)-contingency_table[0,0]-contingency_table[1,0]-contingency_table[0,1]
            oddsratio, pvalue = scipy.stats.fisher_exact(contingency_table)

            if oddsratio==0.0:
                oddsratios[f1_idx,f2_idx]=np.log2(ps)
            else:
                oddsratios[f1_idx,f2_idx]=np.log2(oddsratio) 
            pvals[f1_idx,f2_idx]=pvalue 
            proportion_cells_in_f1_from_f2=1.0*len(overlap)/len(cells_in_f1)
            proportions[f1_idx,f2_idx]=proportion_cells_in_f1_from_f2

    oddsratios_df=pd.DataFrame(oddsratios)
    oddsratios_df.index=f1s
    oddsratios_df.columns=f2s

    pvals_df=pd.DataFrame(pvals)
    pvals_df.index=f1s
    pvals_df.columns=f2s
    if add_min_pval:
        min_pval=np.min(pvals[np.nonzero(pvals)])
    else:
        min_pval=0
    #adjust pvals
    from statsmodels.stats.multitest import multipletests
    pvals_df=np.reshape(multipletests(np.array(pvals_df).flatten(),method=pval_correction)[1],
                        pvals_df.shape)
    pvals_df=-np.log10(pvals_df+min_pval)*np.sign(oddsratios_df)

    proportions_df=pd.DataFrame(proportions)
    proportions_df.index=f1s
    proportions_df.columns=f2s

    pref='enrich_'+f1+'_vs_'+f2
    adata.uns[pref+'.oddsratios']=oddsratios_df
    adata.uns[pref+'.p_adj.negLog10.signed']=pvals_df
    adata.uns[pref+'.proportions']=proportions_df
    
    if copy:
        return(adata)



#this method taken from Dixit et al., 2016, https://github.com/asncd/MIMOSCA/blob/master/GBC_CBC_pairing/fit_moi.ipynb
def moi(adata_here,pref='perturb',level='guide',gridsize=100,maxk=10):
    
    import scipy
    from numpy import unravel_index
    
    print('Computing MOI and detection probability using code from Dixit et al., 2016')

    if pref+'.'+level+'.perturbations_per_cell' not in adata_here.obs:
        print('ERROR: missing adata.obs['+pref+'.'+level+'.perturbations_per_cell], please run perturb.pp.perturbations_per_cell first')
        exit
    if pref+'.'+level+'.cells_per_perturbation' not in adata_here.uns:
        print('ERROR: missing adata.obs['+pref+'.'+level+'.cells_per_perturbation], please run perturb.pp.cells_per_perturbation first')
        exit
        
    moi_dist=np.array(list(adata_here.obs[pref+'.'+level+'.perturbations_per_cell']))
    num_virus=adata_here.uns[pref+'.'+level+'.cells_per_perturbation'].shape[0]

    n,bins=np.histogram(moi_dist,range(int(maxk)+1))

    #maximum number of viruses possible (per cell)
    #maxk

    #total number of unique barcodes
    print('number of distinct perturbations',num_virus)

    #gridsize for performing lambda and alpha search
    nums=gridsize

    #specify start and finishing MOI to search over, it is set to 0.1 and 3 here
    mois=np.linspace(0.1,5,nums) #(0.1,2,nums)
    #specify start and finishing detection probability to search over, it is set to 0.1 and 0.99 here
    detects=np.linspace(0.1,0.99,nums)

    #initialize search array
    LL=np.zeros((nums,nums))


    #loop through square grid of different poission parameters and detection probabilities
    for i in range(nums):
        for m in range(nums):

            #current parameter guesses
            moi_guess=mois[i]
            detect_guess=detects[m]

            #initialize possion distribution with current guess    
            pdf=scipy.stats.poisson.pmf(k=range(maxk),mu=moi_guess)

            #Zero truncation and renormalization
            pdf[0]=0
            pdf=np.divide(pdf,np.sum(pdf))


            #get probabilities after convolving with binomial distribution
            zibpdf=np.zeros((maxk,1))
            for k in range(maxk):
                pf=0
                for j in np.arange(k,maxk):
                    pf+=pdf[j]*scipy.stats.binom.pmf(k,j,detect_guess)
                zibpdf[k]=pf

            #evaluate log likelihood after multiplying with observed values
            ll=1.0
            for k in range(maxk):#range(len(n)):
                ll+=n[k]*np.log(zibpdf[k])
            LL[i,m]=ll

    #Log likelihood vs. paramter space
    plt.contour(np.round(detects,2),np.round(mois,2),LL,400,cmap='magma')
    plt.colorbar()
    plt.xlabel('Detection Probability')
    plt.ylabel('MOI')

    #Find parameters that maximize the log likelihood
    final_tuple=unravel_index(LL.argmax(), LL.shape)
    moi_guess=int(100*mois[final_tuple[0]])/100
    detect_guess=int(100*detects[final_tuple[1]])/100
    print('MOI:',moi_guess)
    print('Detection probability:',detect_guess)

    adata_here.uns['MOI']=moi_guess
    adata_here.uns['Detection_probability']=detect_guess
    plt.scatter(detect_guess,moi_guess,color='black',s=50)
