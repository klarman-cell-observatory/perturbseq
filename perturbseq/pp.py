import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

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

def subsample_cells(adata,num_cells,grouping_variable):
    import random
    cells_keep=[]
    groups=list(set(adata.obs[grouping_variable]))
    for group in groups:
        group_cells=list(adata.obs_names[adata.obs[grouping_variable]==group])
        if len(group_cells)<num_cells:
            print('warning: fewer cells than needed for '+group+'. skipping subsampling')
        else:
            group_cells=random.sample(group_cells,num_cells)
        for cell in group_cells:
            cells_keep.append(cell)
    return(adata[cells_keep,:])

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
                

def annotate_controls(adata_here,control_guides=['unperturbed'],pref='',copy=False):
    #if no controls are specified, unperturbed is the control
    
    if copy: adata_here = adata_here.copy()

    if pref+'guide.total' not in adata_here.obs:
        print('ERROR: '+pref+'guide.total'+' was not found in adata.obs. Please first run perturb.perturb2obs')
        exit
        
    control_anno=[]
    for i in range(adata_here.n_obs):
        guide=adata_here.obs[pref+'guide.total'][i]
        if guide in control_guides:
            control_status='control'
        else:
            control_status='not control'
        control_anno.append(control_status)
        
    adata_here.obs[pref+'control']=control_anno

    if copy:
        return(adata_here)

def remove_guides_from_gene_names(adata_here,pref=''):
    if pref+'cell2guide' not in adata_here.obsm:
        print('ERROR: '+pref+'cell2guide'+' was not found in adata.obsm. Please first run perturb.read_perturbations_csv')
        exit
    guides=adata_here.obsm[pref+'cell2guide'].columns
    guides_in_adata_varnames=list(set(adata_here.var_names).intersection(set(guides)))
    print('filtering out',len(guides_in_adata_varnames),'guide names from the expression matrix')
    if len(guides_in_adata_varnames)>0:
        remaining_varnames=list(set(adata_here.var_names).difference(set(guides_in_adata_varnames)))
        adata_out=adata_here[:,remaining_varnames]
    else:
        adata_out=adata_here
    return(adata_out)

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
