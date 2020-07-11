import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import copy
import numpy as np

def enriched_features(adata,f1,f2,fdr=None,x_or=5,x_p=20,**kwargs):
    
    pref='enrich_'+f1+'_vs_'+f2
    
    for item in ['.oddsratios','.p_adj.negLog10.signed']:
        if pref+item not in adata.uns:
            print("ERROR: result adata.uns['"+pref+item+"'] not found! \nPlease first run perturb.tl.enrich_features(adata,f1='"+f1+"',f2='"+f2+"')")
            return
    
    import seaborn as sns
    item='.oddsratios'
    g=sns.clustermap(adata.uns[pref+item],
                         vmin=-x_or,vmax=x_or,**kwargs)
    plt.title('log2 odds ratio\n('+f1+' vs '+f2+')',loc='center')
    ax = g.ax_heatmap
    ax.set_ylabel(f1)
    ax.set_xlabel(f2)

    
    item='.p_adj.negLog10.signed'
    g=sns.clustermap(adata.uns[pref+item],
                         vmin=-x_p,vmax=x_p,**kwargs)
    plt.title('signed -log10(p adj)\n('+f1+' vs '+f2+')',loc='center')
    ax = g.ax_heatmap
    ax.set_ylabel(f1)
    ax.set_xlabel(f2)

def perturbations_per_cell(adata_here,vmin=0,vmax=5,ax=None,pref='perturb',level='guide',copy=False):
    if ax==None:
        fig,plots=plt.subplots(1)
        ax=plots
    import copy
    vals=copy.deepcopy(adata_here.obs[pref+'.'+level+'.perturbations_per_cell'])
    vals[vals>vmax]=vmax
    vals[vals<vmin]=vmin
    
    ax.hist(vals,bins=range(int(vmin),int(vmax+1)),
                color='black')
    ax.set_xticks([vmin,vmax])
    ax.set_ylabel('Number of cells')
    ax.set_xlabel('Distinct '+level+'s\nper cell')
    ax.grid(False)
    plt.show()

def cells_per_perturbation(adata_here,ax=None,pref='perturb',level='guide',vmin=0,vmax=None):
    if ax==None:
        fig,plots=plt.subplots(1)
        ax=plots
    import copy
    for var in [pref+'.'+level+'.cells_per_perturbation',pref+'.'+level+'.cells_per_perturbation.singly_infected']:
        color='lightgray'
        label='all'
        if var==pref+'.'+level+'.cells_per_perturbation.singly_infected':
            color='black'
            label='single inf.'
        n, bins, patches=ax.hist(adata_here.uns[var]['Number of cells'],
                             color=color,
          cumulative=True,fill=False,histtype='step',label=label)
        patches[0].set_xy(patches[0].get_xy()[:-1])
    ax.grid(False)
    if vmin!=None and vmax!=None:
        ax.set_xlim(vmin,vmax)
    
    ax.set_xlabel('Cells per\n'+level)
    ax.set_ylabel('Number of\n'+level+'s')
    ax.legend(loc='lower right') 

def ranked_cells_per_perturbation(adata_here,figwidth=2,figheight=10,pref='perturb',level='guide',vmin=0,vmax=None):
    if vmin!=None and vmax!=None:
        fig,plots=plt.subplots(1,2)
    else:
        fig,plots=plt.subplots(1,2,sharex=True)
    
    fig.set_size_inches(figwidth,figheight)

    import copy
    names=[pref+'.'+level+'.cells_per_perturbation.singly_infected',pref+'.'+level+'.cells_per_perturbation']
    for i in range(len(names)):
        df=adata_here.uns[names[i]]['Number of cells']
        
        color='lightgray'
        label='all'
        if names[i]==pref+'.'+level+'.cells_per_perturbation.singly_infected':
            color='black'
            label='single inf.'
            sorted_order=df.sort_values().index
        df=df.loc[sorted_order]
        
        plots[i].barh(df.index,df,color=color,label=label)
        plots[i].set_ylim(-0.5,df.shape[0]-0.5)
        plots[i].grid(False)
        if names[i]==pref+'.'+level+'.cells_per_perturbation.singly_infected':
            plots[i].set_yticklabels(sorted_order,fontsize=figheight/df.shape[0]*60)
        else:
            plots[i].set_yticklabels([])
        plots[i].set_title(label,fontsize=figwidth*5)
        plots[i].set_xlabel('Number\nof cells',fontsize=figwidth*5)
        if vmin!=None and vmax!=None:
            plots[i].set_xlim(vmin,vmax)
