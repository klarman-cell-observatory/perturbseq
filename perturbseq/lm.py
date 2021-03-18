import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
import copy
from perturbseq.pp import perturb_overlap_obs,obs_to_design_matrix,split_train_valid_test
#import mimosca 

def make_X_y_covariates(adata_here,
                        model_name='linear_model',
                        include_expression=True,
                        y_obs=[],
                        perturbation_list=None,
                        control_cells=[],
                        covariates_list=[],
                        use_raw=False,
                        keep_unassigned=False):

    if copy: adata_here = adata_here.copy()

    #X_df                                                                                                                         
    if perturbation_list is None:
        perturbation_list=list(set(adata.obs['guide']).difference(['multiple','unassigned']))

    perturbation_list=perturb_overlap_obs(perturbation_list,adata_here,'perturbations')
    X_df=obs_to_design_matrix(adata_here,perturbation_list)

    #set control cells to 0                                                                                                       
    #TODO: think about this more                                                                                                  
    control_cells_in_data=list(set(control_cells).intersection(set(adata_here.obs_names)))
    if len(control_cells_in_data)>0:
        X_df.loc[control_cells_in_data,:]=0

    #y_df         
    #expression data  
    if include_expression:
        if use_raw:
            expression=adata_here.raw.X.toarray()
            genes=adata_here.raw.var_names
        else:
            expression=adata_here.X
            genes=adata_here.var_names
        y_df=pd.DataFrame(expression)
        y_df.index=adata_here.obs_names
        y_df.columns=genes
    y_obs_list=perturb_overlap_obs(y_obs,adata_here,'obs')
    if len(y_obs_list)>0:
        y_obs_df=adata_here.obs[y_obs_list]
        if include_expression:
            y_df=pd.concat([y_df,y_obs_df],axis=1)
        else:
            y_df=y_obs_df

    #covariates_df                                                                                                                                
    covariates_list=perturb_overlap_obs(covariates_list,adata_here,'covariates')
    covariates_df=obs_to_design_matrix(adata_here,covariates_list)

    #whether to keep or remove unassigned                                                                                                                                                                                                                                                            
    if not keep_unassigned:
        keep=list(set(control_cells).union(set(list(X_df.index[X_df.sum(axis=1)>0]))))

    if keep_unassigned:
        keep=list(set(control_cells).union(set(list(X_df.index[X_df.sum(axis=1)>0]))).union(adata_here[adata_here.obs['unassigned']>0,:].obs_names\
))

    X_df=X_df.loc[keep,:]
    y_df=y_df.loc[keep,:]
    covariates_df=covariates_df.loc[X_df.index,:]

    adata_here.obs[model_name+'.X_df']=X_df
    adata_here.obs[model_name+'.y_df']=y_df
    adata_here.obs[model_name+'.covariates_df']=covariates_df

    if copy:
        return(adata_here)


def train_lm(X_df,
             y_df,
             lm,
             training=None):
    
    #convert to arrays
    X=np.array(X_df)
    y=np.array(y_df)
    
    if training is None:
        lm.fit(X,y)
    else:
        lm.fit(X[training,:],y[training,:])
    
    coef=pd.DataFrame(lm.coef_,
                        index=y_df.columns,
                        columns=X_df.columns)
    
    return(coef,lm)
    

def train(adata_here,lm,
        model_name='linear_model',
        perturbations_list=None, 
        include_expression=True, 
        y_obs=[],
        covariates_list=[],
        rng=None,
        adjust=False,
        adjust_vars=[],
        use_raw=True,
        keep_unassigned=False,
        training_proportion=0.8,
        validation_proportion=0.1,
        test_proportion=0.1,
        copy=False):
    
    if copy: adata_here = adata_here.copy()
        
    #perturbation list
    if perturbations_list is None:
        perturbations_list=list(set(','.join(list(set(adata_here.obs['guide']).difference(['unassigned']))).split(',')))
        print(perturbations_list)
        
    #split into train/test unless already done
    if 'train_test' not in adata_here.obs:
        split_train_valid_test(adata_here,
                               training_proportion=training_proportion,
                               validation_proportion=validation_proportion,
                               test_proportion=test_proportion,
                               rng=rng,
                               copy=copy)
        
    
    #first, make datasets
    ###X_df,y_df,covariates_df=
    make_X_y_covariates(adata_here,
                        model_name=model_name,
                        include_expression=include_expression,
                                y_obs=y_obs,
                               perturbation_list=perturbations_list,
                               covariates_list=covariates_list,
                               use_raw=use_raw,
                               keep_unassigned=keep_unassigned)
    
    #check data
    assert (X_df.index==y_df.index).all()
    if covariates_df.__class__.__name__!="NoneType":
        assert (X_df.index==covariates_df.index).all()
    
    if covariates_df.__class__.__name__!="NoneType":
        X_plus_covariates=np.concatenate([np.array(X_df),np.array(covariates_df)],axis=1)
        X_plus_covariates_columns=list(X_df.columns)+list(covariates_df.columns)
    else:
        X_plus_covariates=np.array(X_df)
        X_plus_covariates_columns=list(X_df.columns)
    #subset of the data with the desired perturbations
    adata_subset=adata_here[X_df.index,:]
    training=[i for i in range(adata_subset.n_obs) if adata_subset.obs['train_test'][i]=='train']
    
    #train model
    print('\nFitting model\n',lm)
    coef,trained_lm=train_lm(pd.DataFrame(X_plus_covariates,
                                          index=X_df.index,
                                          columns=X_plus_covariates_columns),
                             y_df,
                             lm,
                             training=training)
    
    #adjust, if requested
    if not adjust and len(adjust_vars)>0:
        print('WARNING: you provided variables to adjust with EM but also set adjust=False. No adjustment will be done.')
    adjust_vars_idx=[]
    if adjust==True:
        #setup the adjustment variables
        for idx in range(X_df.shape[1]):
            if X_df.columns[idx] in adjust_vars:
                adjust_vars_idx.append(idx)
        #get adjusted X (this has the covariates also, though they have not been adjusted)
        X_adjust=np.array(mimosca.bayes_cov_col(y_df,
                                                pd.DataFrame(X_plus_covariates,
                                                            index=X_df.index),
                                                adjust_vars_idx,
                                                trained_lm))        
        print('Re-fitting model on adjusted data\n',lm)
        coef,trained_lm=train_lm(pd.DataFrame(X_adjust,
                                              index=X_df.index,
                                              columns=X_plus_covariates_columns),
                                 y_df,
                                 lm,
                                 training=None)    
        
        
    #save items into adata
    adata_here.uns[model_name+'.coef']=coef
    adata_here.uns[model_name+'.y']=y_df
    adata_here.uns[model_name+'.X']=pd.DataFrame(X_plus_covariates,
                                                index=X_df.index,
                                                columns=X_plus_covariates_columns)
    if adjust:
        adata_here.uns[model_name+'.X_adjust']=pd.DataFrame(X_adjust,
                                                index=X_df.index,
                                                columns=X_plus_covariates_columns)
    
    if copy:
        return(adata_here)
