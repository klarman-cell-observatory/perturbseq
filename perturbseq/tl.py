import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

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
