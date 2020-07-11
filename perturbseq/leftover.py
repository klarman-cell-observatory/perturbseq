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


def perturb2obs(adata_here,pref='perturb',copy=False):
    if pref+'.cell2guide' not in adata_here.obsm:
        print('ERROR: '+pref+'.cell2guide'+' was not found in adata.obsm. Please first run perturb.read_perturbations_csv')
        exit
    else:
        if copy: adata_here = adata_here.copy()

        #go through each cell and make 2 obs.
        #1. the guide combination in the cell
        #2. annotate cells with multiple guides just as multiple
        anno_total_guide=[]
        anno_singles_guide=[]
        cell2guide=adata_here.obsm[pref+'.cell2guide']
        if pref+'.cell2gene' in adata_here.obsm:
            cell2gene=adata_here.obsm[pref+'.cell2gene']
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
        
        
        adata_here.obs[pref+'.guide.total']=anno_total_guide
        adata_here.obs[pref+'.guide']=anno_singles_guide
        adata_here.obs[pref+'.gene.total']=anno_total_gene
        adata_here.obs[pref+'.gene']=anno_singles_gene
        if copy:
            return(adata_here)
                

def annotate_controls(adata_here,control_guides=['unperturbed'],pref='perturb',copy=False):
    #if no controls are specified, unperturbed is the control
    
    if copy: adata_here = adata_here.copy()

    if pref+'.guide.total' not in adata_here.obs:
        print('ERROR: '+pref+'.guide.total'+' was not found in adata.obs. Please first run perturb.perturb2obs')
        exit
        
    control_anno=[]
    for i in range(adata_here.n_obs):
        guide=adata_here.obs[pref+'.guide.total'][i]
        if guide in control_guides:
            control_status='control'
        else:
            control_status='not control'
        control_anno.append(control_status)
        
    adata_here.obs[pref+'.control']=control_anno

    if copy:
        return(adata_here)

def remove_guides_from_gene_names(adata_here,pref='perturb'):
    if pref+'.cell2guide' not in adata_here.obsm:
        print('ERROR: '+pref+'.cell2guide'+' was not found in adata.obsm. Please first run perturb.read_perturbations_csv')
        exit
    guides=adata_here.obsm[pref+'.cell2guide'].columns
    guides_in_adata_varnames=list(set(adata_here.var_names).intersection(set(guides)))
    print('filtering out',len(guides_in_adata_varnames),'guide names from the expression matrix')
    if len(guides_in_adata_varnames)>0:
        remaining_varnames=list(set(adata_here.var_names).difference(set(guides_in_adata_varnames)))
        adata_out=adata_here[:,remaining_varnames]
    else:
        adata_out=adata_here
    return(adata_out)

#========= perturbation stats
 
def perturbations_per_cell(adata_here,pref='perturb',level='guide',copy=False):
    
    if copy: adata_here = adata_here.copy()
    
    if pref+'.cell2'+level not in adata_here.obsm:
        print('ERROR: '+pref+'.cell2'+level+' was not found in adata.obsm. Please first run perturb.read_perturbations_csv')
        exit
    pe_df=adata_here.obsm[pref+'.cell2'+level]>0.0
    perturbations_per_cell_val=pe_df.sum(axis=1)
    adata_here.obs[pref+'.'+level+'.perturbations_per_cell']=perturbations_per_cell_val
    
    if copy:
        return(adata_here)

def cells_per_perturbation(adata_here,level='guide',pref='perturb',copy=False):
    
    if copy: adata_here = adata_here.copy()

    if pref+'.cell2'+level not in adata_here.obsm:
        print('ERROR: '+pref+'.cell2'+level+' was not found in adata.obsm. Please first run perturb.read_perturbations_csv')
        exit
    pe_df=adata_here.obsm[pref+'.cell2'+level]>0.0
    
    cells_per_perturbation_val=pe_df.sum(axis=0)
    adata_here.uns[pref+'.'+level+'.cells_per_perturbation']=pd.DataFrame(cells_per_perturbation_val,index=cells_per_perturbation_val.index) 
    adata_here.uns[pref+'.'+level+'.cells_per_perturbation'].columns=['Number of cells']
    #also restrict to singly infected cells
    perturbations_per_cell_val=pe_df.sum(axis=1)
    singles_df=pe_df.loc[perturbations_per_cell_val==1,:]
    cells_per_perturbation_val_singles=singles_df.sum(axis=0)
    adata_here.uns[pref+'.'+level+'.cells_per_perturbation.singly_infected']=pd.DataFrame(cells_per_perturbation_val_singles,index=cells_per_perturbation_val_singles.index)
    adata_here.uns[pref+'.'+level+'.cells_per_perturbation.singly_infected'].columns=['Number of cells']
    
    if copy:
        return(adata_here)

#this method taken from Dixit et al., 2016, https://github.com/asncd/MIMOSCA/blob/master/GBC_CBC_pairing/fit_moi.ipynb
def moi(adata_here,pref='perturb',level='guide',gridsize=100,maxk=10):
    
    import scipy
    from numpy import unravel_index
    
    print('Computing MOI and detection probability using code from Dixit et al., 2016')

    if pref+'.'+level+'.perturbations_per_cell' not in adata.obs:
        print('ERROR: missing adata.obs['+pref+'.'+level+'.perturbations_per_cell], please run perturb.perturbations_per_cell first')
        exit
    if pref+'.'+level+'.cells_per_perturbation' not in adata.uns:
        print('ERROR: missing adata.obs['+pref+'.'+level+'.cells_per_perturbation], please run perturb.cells_per_perturbation first')
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
