import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from perturbseq.util import display_progress

def _get_perturbations(adata_here,
                     perturbations_obs='guide'):
    
    #get the list of perturbations
    perturbations_list=list(set(adata_here.obs[perturbations_obs]).intersection(set(list(adata_here.obs.columns))).difference(['unassigned']))
    
    perturbations_list.sort()
    return(perturbations_list)

def get_perturbations(adata_here,
                     perturbations_obs='guide',
                     copy=False):
    
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
def check_id_list(perturbation_list,adata_here,list_type):
    found_perturbations=[]
    count_perturbations=0
    for perturbation in perturbation_list:
        if perturbation not in adata_here.obs:
            print('Warning: '+perturbation+' is not in the provided dataset and will be ignored')
        else:
            count_perturbations+=1
            found_perturbations.append(perturbation)
    print('Found '+str(count_perturbations)+'/'+str(len(perturbation_list))+' '+list_type)

    return(found_perturbations)

def multiple_annotations_to_design_matrix(adata_here,annotation_names,
                                          binarize=True):

    #say that we expect numerical data here (0/1)                                                                                                                  
    keep_anno_names=[]
    for i in range(len(annotation_names)):
        anno=annotation_names[i]
        if anno not in adata_here.obs:
            print('WARNING: '+anno+' not in the available annotations. Please add it and re-run')
        else:
            keep_anno_names.append(anno)

    annotations=keep_anno_names
    cells=list(adata_here.obs_names)
    design_matrix=np.zeros((len(cells),len(annotations)))

    for cell_idx in range(len(adata_here.obs_names)):

        if cell_idx%1000==0:
            display_progress(cell_idx,len(adata_here.obs_names))

        for anno_idx in range(len(annotations)):
            anno=annotations[anno_idx]
            design_matrix[cell_idx,anno_idx]=adata_here.obs[anno][cell_idx]
    display_progress(1,1)
    print('\n')
    design_matrix_df=pd.DataFrame(design_matrix)
    if binarize:
        design_matrix_df=(design_matrix_df>0.0)*1.0
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
    #add columns with the interaction terms
    for interaction_term in interaction_terms:
        interaction_columns=interaction_term.split(',')
        values=design_matrix_df.loc[:,interaction_columns].prod(axis=1)
        import copy
        design_matrix_df[interaction_term]=copy.deepcopy(values)

    return(design_matrix_df)


def split_train_valid_test(adata_here,training_proportion=0.6,validation_proportion=0.2,test_proportion=0.2,rng=None,copy=False):
    assert training_proportion<=1.0
    assert validation_proportion<=1.0
    assert test_proportion<=1.0
    assert (training_proportion+validation_proportion+test_proportion)<=1.0
    
    if copy: adata_here = adata_here.copy()
    
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
                               'train_test':'train'},index=adata_here.obs_names)
    train_test_df=train_test_df.iloc[idx_shuff,:]
    train_test_df.iloc[training,1]='train'
    train_test_df.iloc[validation,1]='valid'
    train_test_df.iloc[test,1]='test'
    adata_here.obs['train_test']=train_test_df.loc[adata_here.obs_names,'train_test']

    if copy:
        return(adata_here)

def filter_multiplets(adata_here,pref='',level='guide',keep_unassigned=True,
                           copy=False):
    
    import pandas as pd

    if copy: adata_here = adata_here.copy()
        
    #===============
    if pref+'guide.perturbations_per_cell' not in adata_here.obs:
        perturb.pp.perturbations_per_cell(adata_here,level=level)
    keep=adata_here.obs_names[adata_here.obs[pref+'guide.perturbations_per_cell']==1]
    if keep_unassigned:
        unassigned=adata_here.obs_names[adata_here.obs[pref+'guide']=='unassigned']
        keep=list(set(keep).union(set(unassigned)))  
    adata_here._inplace_subset_obs(keep)
    
    #===============
        
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

def subsample_cells(adata_here,num_cells,grouping_variable,
                   rng=np.random.RandomState(1234)):

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
                

def annotate_controls(adata_here,control_guides=[],pref='',copy=False):
    #if no controls are specified, nothing gets assigned
    
    if copy: adata_here = adata_here.copy()

    if pref+'guide' not in adata_here.obs:
        print('ERROR: '+pref+'guide'+' was not found in adata.obs. Please first perturb.io.run read_perturbations_csv')
        exit
        
    control_anno=[]
    for i in range(adata_here.n_obs):
        guide=adata_here.obs[pref+'guide'][i]
        if guide in control_guides:
            control_status='control'
        else:
            control_status='not control'
        control_anno.append(control_status)
        
    adata_here.obs[pref+'control']=control_anno

    if copy:
        return(adata_here)

def remove_guides_from_gene_names(adata_here,pref='',copy=False):

    if copy: adata_here = adata_here.copy()

    guides=list(set(adata_here.obs['guide']).difference(['unassigned','multiple']))
    guides_in_adata_varnames=list(set(adata_here.var_names).intersection(set(guides)))
    print('filtering out',len(guides_in_adata_varnames),'guide names from the expression matrix')
    if len(guides_in_adata_varnames)>0:
        remaining_varnames=list(set(adata_here.var_names).difference(set(guides_in_adata_varnames)))
        adata_here._inplace_subset_var(remaining_varnames)
    else:
        adata_here=adata_here
    if copy:
        return(adata_here)

#========= perturbation stats
 
def perturbations_per_cell(adata_here,pref='',level='guide',copy=False):
    
    if copy: adata_here = adata_here.copy()
    
    if pref+'cell2'+level not in adata_here.obsm:
        print('ERROR: '+pref+'cell2'+level+' was not found in adata.obsm. Please first run perturb.read_perturbations_csv')
        exit
    pe_df=adata_here.obsm[pref+'cell2'+level]>0.0
    perturbations_per_cell_val=pe_df.sum(axis=1)
    adata_here.obs[pref+level+'.perturbations_per_cell']=perturbations_per_cell_val
    
    if copy:
        return(adata_here)

def cells_per_perturbation(adata_here,level='guide',pref='',copy=False):

    """Counts the number of cells for each perturbation.

    Args:
        adata_here: adata
        level: guide

    Returns:
        something: something
    """
    
    if copy: adata_here = adata_here.copy()

    if pref+'cell2'+level not in adata_here.obsm:
        print('ERROR: '+pref+'cell2'+level+' was not found in adata.obsm. Please first run perturb.read_perturbations_csv')
        exit
    pe_df=adata_here.obsm[pref+'cell2'+level]>0.0
    
    cells_per_perturbation_val=pe_df.sum(axis=0)
    adata_here.uns[pref+level+'.cells_per_perturbation']=pd.DataFrame(cells_per_perturbation_val,index=cells_per_perturbation_val.index) 
    adata_here.uns[pref+level+'.cells_per_perturbation'].columns=['Number of cells']
    #also restrict to singly infected cells
    perturbations_per_cell_val=pe_df.sum(axis=1)
    singles_df=pe_df.loc[perturbations_per_cell_val==1,:]
    cells_per_perturbation_val_singles=singles_df.sum(axis=0)
    adata_here.uns[pref+level+'.cells_per_perturbation.singly_infected']=pd.DataFrame(cells_per_perturbation_val_singles,index=cells_per_perturbation_val_singles.index)
    adata_here.uns[pref+level+'.cells_per_perturbation.singly_infected'].columns=['Number of cells']
    
    if copy:
        return(adata_here)

#this method taken from Dixit et al., 2016, https://github.com/asncd/MIMOSCA/blob/master/GBC_CBC_pairing/fit_moi.ipynb
def moi(adata_here,pref='',level='guide',gridsize=100,maxk=10):
    
    import scipy
    from numpy import unravel_index
    
    print('Computing MOI and detection probability using code from Dixit et al., 2016')

    if pref+level+'.perturbations_per_cell' not in adata.obs:
        print('ERROR: missing adata.obs['+pref+level+'.perturbations_per_cell], please run perturb.pp.perturbations_per_cell first')
        exit
    if pref+level+'.cells_per_perturbation' not in adata.uns:
        print('ERROR: missing adata.obs['+pref+level+'.cells_per_perturbation], please run perturb.pp.cells_per_perturbation first')
        exit
        
    moi_dist=np.array(list(adata_here.obs[pref+level+'.perturbations_per_cell']))
    num_virus=adata_here.uns[pref+level+'.cells_per_perturbation'].shape[0]

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
