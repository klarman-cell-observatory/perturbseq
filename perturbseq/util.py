import sys
import pandas as pd
import numpy as np
import copy

def corr_mat(df,axis=0,corr_type='spearman'):
    
    from scipy.stats import spearmanr
    from scipy.stats import pearsonr
    
    if axis==0: #do it for rows
        df_here=np.array(copy.deepcopy(df))
    if axis==1:
        df_here=np.array(copy.deepcopy(df.T))
        
    corr_mat=np.zeros((df_here.shape[0],df_here.shape[0]))
    for i in range(df_here.shape[0]):
        #if i%10==0:
         #   display_progress(i,df_here.shape[0])
        for j in range(i,df_here.shape[0]):
            a=df_here[i,:]
            b=df_here[j,:]
            if np.std(a)==0 or np.std(b)==0:
                if corr_type in ['pearson','spearman']:
                    continue
            if corr_type=='pearson':
                val=pearsonr(a,b)[0]
            if corr_type=='spearman':
                val=spearmanr(a,b)[0]
            if corr_type=='diff':
                val=np.mean(np.abs(a-b))
            corr_mat[i,j]=val
            corr_mat[j,i]=val
    if axis==0: #do it for rows
        corr_mat=pd.DataFrame(corr_mat,index=df.index,columns=df.index)
    if axis==1: #do it for rows
        corr_mat=pd.DataFrame(corr_mat,index=df.columns,columns=df.columns)
    return(corr_mat)

def get_perturbations(adata_here,pref='',compact=False,level='guide',
                           copy=False):
    
    import pandas as pd

    if not compact:
        perturbs=list(set(adata_here.obs[pref+level]).difference(['unassigned']))
    elif compact:
        perturbs=list(set(adata_here.obs[pref+level+'.compact']).difference(['unassigned','multiple']))
    perturbs.sort()
    return(perturbs)

def display_progress(cur_num,max_num):
    sys.stdout.write('\r'+str(int(100*(1.0*cur_num/max_num)))+' %')
    sys.stdout.flush() 

def end_progress(num_spaces,newline=False):
    text="\r"
    for i in range(num_spaces):
        text=text+' '
    if newline:
        text=text+'\n'
    sys.stdout.write(text)
