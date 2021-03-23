import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from perturbseq.util import display_progress

def score_programs(adata_here,program_name='bulk.guide.program',
                  pref=None,copy=False):
    
    """Get program score for each cell, by averaging program genes

    """

    if copy: adata_here = adata_here.copy()
    
    if pref is None:
        pref=program_name
    programs=list(set(adata_here.var[program_name]))
    for pro in programs:
        print('scoring',pro)
        pro_genes=adata_here.var_names[adata_here.var[program_name]==pro]
        adata.obs[pref+str(pro)]=adata[:,pro_genes].X.mean(axis=1)
        
    if copy:
        return(adata_here)

def obs_mean(adata_here,grouping_variable,
             obs,outpref='obs_mean',copy=False):

    """Get the mean of an obs across pre-defined groups
    
    """
    
    if grouping_variable not in adata_here.obs.columns:
        print("ERROR: Variable '"+grouping_variable+"' not in adata.obs")
        return()
    for obs_variable in obs:
        if obs_variable not in adata_here.obs.columns:
            print("ERROR: Variable '"+obs_variable+"' not in adata.obs")
            return()
            
    
    if copy: adata_here = adata_here.copy()
        
    #construct the profiles
    import copy
    profiles=list(set(adata_here.obs[grouping_variable]))
    profile_obs=copy.deepcopy(adata_here.obs[grouping_variable])
    
    profile_matrix=np.zeros((len(profiles),len(obs)))
    for profile_idx in range(len(profiles)):
        profile=profiles[profile_idx]
        cells_with_profile=list(adata_here.obs_names[profile_obs==profile])
        data_profile=adata_here[cells_with_profile,:].obs.loc[:,obs]
        profile_matrix[profile_idx,:]=data_profile.mean(axis=0)
    profile_matrix_df=pd.DataFrame(profile_matrix)
    profile_matrix_df.index=profiles
    profile_matrix_df.columns=obs
    adata_here.uns[outpref]=profile_matrix_df
        
    if copy:
        return(adata_here)

def _get_perturbations(adata_here,
                     perturbations_obs='guide'):
    
    #get the list of perturbations
    perturbations_list=list(set(adata_here.obs[perturbations_obs]).intersection(set(list(adata_here.obs.columns))).difference(['unassigned']))
    
    perturbations_list.sort()
    return(perturbations_list)

def get_perturbations(adata_here,
                     perturbations_obs='guide',
                     copy=False):
    
    """Get a list of perturbations in the dataset

    It will be stored as adata.uns['PS.'+perturbations_obs+'.list']
    """

    if copy: adata_here = adata_here.copy()
        
    #check if perturbations_obs in adata
    try:
        assert perturbations_obs in adata_here.obs.columns
    except AssertionError:
            print('ERROR: "'+perturbations_obs+'" is not in adata.obs.')
            return
        
    #get perturbations
    if 'PS.'+perturbations_obs+'.list' in adata_here.uns:
        print('WARNING: Over-writing '+'"PS.'+perturbations_obs+'.list"')
    adata_here.uns['PS.'+perturbations_obs+'.list']=_get_perturbations(adata_here,
                                                                       perturbations_obs=perturbations_obs)
    
    if copy:
        return(adata_here)


def downsample_counts(adata_here,
                      downsampling_prob,
                      my_rng=np.random.RandomState(1234)):
    """Downsample counts per cell to a fraction

    """

    import time
    from scipy.sparse import csr_matrix
    from scipy.sparse import coo_matrix
    
    if downsampling_prob>1 or downsampling_prob<0:
        print('WARNING: You provided a downsampling probability outside the range [0,1].')
        print('No downsampling will be performed')
        return()
    
    if downsampling_prob==1.0:
        print('WARNING: You provided a downsampling probability equal to 1. No downsampling will be performed')
        return()
    
    #quickly loop through the sparse matrix of counts
    downsampled_vals=[]
    m=adata_here.X
    nonzeros=m.nonzero()
    nonzero_r=nonzeros[0]
    nonzero_c=nonzeros[1]
    start = time.time()
    
    m.eliminate_zeros()
    original_vals=m.data
    num_elts=len(original_vals)
    
    m_downsampled_data=[]
    
    print('median reads in original dataset',
         np.median(np.array(m.sum(axis=1))))

    elt=0
    while elt<num_elts:
        if elt%5000000==0:
            end = time.time()
            print(int(100*elt/num_elts),'% done')
            start = time.time()
            
        m_downsampled_data.append(my_rng.binomial(original_vals[elt],
                                                 downsampling_prob,1)[0])
        elt+=1
    print(int(100*elt/num_elts),'% done')
    
    downsampled=csr_matrix((m_downsampled_data, 
                            m.indices, 
                            m.indptr), 
                           dtype=float,shape=m.shape)
    
    print('median reads in downsampled dataset',
          np.median(np.array(downsampled.sum(axis=1))))
    return(downsampled)

#preparing for linear model
#==========================
def perturb_overlap_obs(perturbation_list,adata_here,list_name):

    """Get perturbations present as obs

    """

    found_perturbations=[]
    count_perturbations=0
    for perturbation in perturbation_list:
        if perturbation not in adata_here.obs:
            print('Warning: '+perturbation+' is not in the provided dataset and will be ignored')
        else:
            count_perturbations+=1
            found_perturbations.append(perturbation)
    print('Found '+str(count_perturbations)+'/'+str(len(perturbation_list))+' '+list_name)

    return(found_perturbations)

def obs_to_design_matrix(adata_here,obs_names,
                                          binarize=True,covariate=False):

    """Design matrix for linear model from obs

    """

    #say that we expect numerical data here (0/1)
                                                                      
    keep_obs_names=[]
    for i in range(len(obs_names)):
        obs=obs_names[i]
        if obs not in adata_here.obs:
            print('WARNING: '+obs+' not in the available individual annotations. Please add it and re-run')
        else:
            keep_obs_names.append(obs)

    annotations=keep_obs_names
    cells=list(adata_here.obs_names)
    design_matrix_df=adata_here.obs.loc[cells,annotations]
    
    if binarize:
        design_matrix_df=(design_matrix_df.astype(float)>0.0)*1.0
    design_matrix_df.index=cells
    design_matrix_df.columns=annotations

    #go through all the cells, and figure out what combinations of perturbations there are                              
    #add these to the design matrix                                           
                                          
    design_matrix_df_uniq=design_matrix_df.drop_duplicates()
    column_names=design_matrix_df_uniq.columns
    interaction_terms=[]
    for i in range(design_matrix_df_uniq.shape[0]):
        current_columns=[]
        for j in range(len(column_names)):
            if design_matrix_df_uniq.iloc[i,j]>0:
                current_columns.append(column_names[j])
        if len(current_columns)>1:
            current_columns.sort()
            current_columns_join=','.join(current_columns)
            interaction_terms.append(current_columns_join)
    if not covariate:
        #add columns with the interaction terms
        for interaction_term in interaction_terms:
            interaction_columns=interaction_term.split(',')
            values=design_matrix_df.loc[:,interaction_columns].prod(axis=1)
            import copy
            design_matrix_df[interaction_term]=copy.deepcopy(values)

    return(design_matrix_df)

def split_train_valid_test(adata_here,
                           training_proportion=0.6,
                           validation_proportion=0.2,
                           test_proportion=0.2,
                           rng=None,copy_adata=False):

    """Split cells into training, validation and test

    """
    
    assert training_proportion<=1.0
    assert validation_proportion<=1.0
    assert test_proportion<=1.0
    assert (training_proportion+validation_proportion+test_proportion)<=1.0

    num_examples=adata_here.n_obs

    if rng==None:
        idx_shuff=np.random.RandomState(seed=77).permutation(range(num_examples))
    else:
        idx_shuff=rng.permutation(range(num_examples))

    training_threshold=int(num_examples*training_proportion)
    validation_threshold=int(num_examples*(training_proportion+validation_proportion))

    training=range(training_threshold)
    validation=range(training_threshold,min(validation_threshold,num_examples))
    test=range(validation_threshold,num_examples)

    #make obs with train, validation, test                                                               
    train_test_df=pd.DataFrame({'cell':adata_here.obs_names,
                               'train_valid_test':'train'},index=adata_here.obs_names)
    train_test_df=train_test_df.iloc[idx_shuff,:]
    train_test_df.iloc[training,1]='train'
    train_test_df.iloc[validation,1]='valid'
    train_test_df.iloc[test,1]='test'
    print('splitting',train_test_df.loc[adata_here.obs_names,'train_valid_test'].value_counts())
    return(train_test_df.loc[adata_here.obs_names,'train_valid_test'])


def subset_singly_perturbed(adata_here,perturbations_obs='guide',keep_unassigned=True,
                           copy=False):

    """Keep only cells with one perturbation

    """
    
    import pandas as pd

    if copy: adata_here = adata_here.copy()
        
    #===============
    perturbs_per_cell(adata_here,perturbations_obs=perturbations_obs)
    keep=adata_here.obs_names[adata_here.obs['perturbs_per_cell.'+perturbations_obs]==1]
    if keep_unassigned:
        unassigned=adata_here.obs_names[adata_here.obs[perturbations_obs]=='unassigned']
        keep=list(set(keep).union(set(unassigned)))
    if not copy:
        adata_here._inplace_subset_obs(keep)
    #===============
        
    if copy:
        return(adata_here)

def compute_TPT(gbcs_dataset):
    
    """Compute transcript per transcript

    input: pandas data frame with the columns "cbc", "umi", "gbc", "r2" where every row is a read 
    output: pandas data frame with the columns "gbc", "cbc", "umi", "cbc-umi-r2-count", "cbc-umi-count", "TPT"
    NOTE: for the input, multiple reads corresponding to the same cbc-umi combination should be listed as separate lines!
    """

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

def subsample_cells(adata_here,num_cells,grouping_variable,
                   rng=np.random.RandomState(1234)):

    """Subsample cells/perturbation

    """

    import copy
    cells_keep=[]
    groups=list(set(adata_here.obs[grouping_variable]))
    for group in groups:
        group_cells=list(adata_here.obs_names[adata_here.obs[grouping_variable]==group])
        if len(group_cells)<num_cells:
            print('warning: fewer cells than needed for '+group+'. skipping subsampling')
        else:
            group_cells=rng.choice(group_cells,num_cells,replace=False)
        for cell in group_cells:
            cells_keep.append(cell)
    return(adata_here[cells_keep,:])

#todo: work on it, or delete
def perturb2obs(adata_here,pref='',copy=False):
    if pref+'cell2guide' not in adata_here.obsm:
        print('ERROR: '+pref+'cell2guide'+' was not found in adata.obsm. Please first run perturb.read_perturbations_csv')
        exit
    else:
        if copy: adata_here = adata_here.copy()

        #go through each cell and make 2 obs.
        #1. the guide combination in the cell
        #2. annotate cells with multiple guides just as multiple
        anno_total_guide=[]
        anno_singles_guide=[]
        cell2guide=adata_here.obsm[pref+'cell2guide']
        if pref+'cell2gene' in adata_here.obsm:
            cell2gene=adata_here.obsm[pref+'cell2gene']
            anno_total_gene=[]
            anno_singles_gene=[]
        for i in range(adata_here.n_obs):
            cell=adata_here.obs_names[i]
            
            #get the guides
            guides_in_cell=sorted(list(cell2guide.loc[cell,cell2guide.loc[cell,:]>0].index))
            if len(guides_in_cell)==0:
                anno_singles_here='unperturbed'
                anno_total_here='unperturbed'
            elif len(guides_in_cell)==1:
                anno_singles_here=guides_in_cell[0]
                anno_total_here=guides_in_cell[0]
            elif len(guides_in_cell)>1:
                anno_singles_here='multiple'
                anno_total_here=','.join(guides_in_cell)
            else:
                print('Negative number of guides for cell',cell)
                
            anno_total_guide.append(anno_total_here)
            anno_singles_guide.append(anno_singles_here)
            
            #get the genes
            genes_in_cell=sorted(list(cell2gene.loc[cell,cell2gene.loc[cell,:]>0].index))
            if len(genes_in_cell)==0:
                anno_singles_here_gene='unperturbed'
                anno_total_here_gene='unperturbed'
            elif len(genes_in_cell)==1:
                anno_singles_here_gene=genes_in_cell[0]
                anno_total_here_gene=genes_in_cell[0]
            elif len(genes_in_cell)>1:
                anno_singles_here_gene='multiple'
                anno_total_here_gene=','.join(genes_in_cell)
            else:
                print('Negative number of guides for cell',cell)
                
            anno_total_gene.append(anno_total_here_gene)
            anno_singles_gene.append(anno_singles_here_gene)
        
        
        adata_here.obs[pref+'guide.total']=anno_total_guide
        adata_here.obs[pref+'guide']=anno_singles_guide
        adata_here.obs[pref+'gene.total']=anno_total_gene
        adata_here.obs[pref+'gene']=anno_singles_gene
        if copy:
            return(adata_here)
                

def annotate_controls(adata_here,perturbations_obs,control_guides=[],pref='',copy=False):
    
    """Make obs with control guides

    Will be saved as adata.obs[pref+'control']
    """

    #if no controls are specified, nothing gets assigned
    
    if copy: adata_here = adata_here.copy()

    if pref+'guide' not in adata_here.obs:
        print('ERROR: '+pref+'guide'+' was not found in adata.obs. Please first perturb.io.run read_perturbations_csv')
        exit
        
    control_anno=[]
    for i in range(adata_here.n_obs):
        guide=adata_here.obs[perturbations_obs][i]
        if guide in control_guides:
            control_status='control'
        else:
            control_status='not control'
        control_anno.append(control_status)
        
    adata_here.obs[pref+'control']=control_anno

    if copy:
        return(adata_here)

def delete_guides_from_varnames(adata_here,perturbations_obs='guide',pref='',copy=False):

    """Delete perturbations from adata.var_names

    
    """

    if copy: adata_here = adata_here.copy()

    guides=list(set(adata_here.obs[perturbations_obs]).difference(['unassigned','multiple']))
    guides_in_adata_varnames=list(set(adata_here.var_names).intersection(set(guides)))
    print('filtering out',len(guides_in_adata_varnames),'guide names from the expression matrix')
    if len(guides_in_adata_varnames)>0 and not copy:
        remaining_varnames=list(set(adata_here.var_names).difference(set(guides_in_adata_varnames)))
        adata_here._inplace_subset_var(remaining_varnames)
    
    if copy:
        return(adata_here[:,remaining_varnames])

#========= perturbation stats
 
def perturbs_per_cell(adata_here,perturbations_obs='guide',copy=False):
    
    """Number of perturbations per cell

    Stored as adata.uns['perturbs_per_cell.'+perturbations_obs]
    """

    if copy: adata_here = adata_here.copy()
    
    #get perturbations
    perturbations=_get_perturbations(adata_here,
                                     perturbations_obs=perturbations_obs)
    #find their obs
    perturbations=perturb_overlap_obs(perturbations,adata_here,list_name='')
    #then count the nonzero obs
    pe_df=adata_here.obs.loc[:,perturbations]>0.0
    perturbations_per_cell_val=pe_df.sum(axis=1)
    adata_here.obs['perturbs_per_cell.'+perturbations_obs]=perturbations_per_cell_val
    
    if copy:
        return(adata_here)

def cells_per_perturb(adata_here,perturbations_obs='guide',copy=False):

    """Counts the number of cells for each perturbation.

    """
    
    if copy: adata_here = adata_here.copy()

    #get perturbations                           
    perturbations=_get_perturbations(adata_here,
                                     perturbations_obs=perturbations_obs)
    #find their obs                                                                                      
    perturbations=perturb_overlap_obs(perturbations,adata_here,list_name='')

    cell2perturbs=obs_to_design_matrix(adata_here, perturbations, binarize=True, covariate=False)
    cell2perturbs_single=cell2perturbs.loc[:,perturbations]
    cell2perturbs_counts=cell2perturbs.sum(axis=0)
    cell2perturbs_single_counts=cell2perturbs_single.sum(axis=0)
    adata_here.uns['cells_per_perturb.'+perturbations_obs+'.incl_multi_inf']=pd.DataFrame(cell2perturbs_counts,index=cell2perturbs_counts.index,columns=['Number of cells'])
    adata_here.uns['cells_per_perturb.'+perturbations_obs]=pd.DataFrame(cell2perturbs_single_counts,index=cell2perturbs_single_counts.index,columns=['Number of cells'])

    if copy:
        return(adata_here)

