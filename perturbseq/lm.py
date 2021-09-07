import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
import copy
from perturbseq.pp import perturb_overlap_obs,obs_to_design_matrix,split_train_valid_test
#import mimosca 

#from MIMOSCA
#https://github.com/asncd/MIMOSCA/blob/8966411848c75e935abd5b86e77f6bb00acea2b4/mimosca.py
def bayes_cov_col(Y,X,cols,lm):
    """
    @Y    = Expression matrix, cells x x genes, expecting pandas dataframe
    @X    = Covariate matrix, cells x covariates, expecting pandas dataframe
    @cols = The subset of columns that the EM should be performed over, expecting list
    @lm   = linear model object
    """

    #EM iterateit
    Yhat=pd.DataFrame(lm.predict(X))
    Yhat.index=Y.index
    Yhat.columns=Y.columns
    SSE_all=np.square(Y.subtract(Yhat))
    X_adjust=X.copy()


    df_SSE   = []
    df_logit = []

    for curcov in cols:

        curcells=X[X[curcov]>0].index

        if len(curcells)>2:

            X_notcur=X.copy()
            X_notcur[curcov]=[0]*len(X_notcur)

            X_sub=X_notcur.loc[curcells]

            Y_sub=Y.loc[curcells]

            GENE_var=2.0*Y_sub.var(axis=0)
            vargenes=GENE_var[GENE_var>0].index

            Yhat_notcur=pd.DataFrame(lm.predict(X_sub))
            Yhat_notcur.index=Y_sub.index
            Yhat_notcur.columns=Y_sub.columns

            SSE_notcur=np.square(Y_sub.subtract(Yhat_notcur))
            SSE=SSE_all.loc[curcells].subtract(SSE_notcur)
            SSE_sum=SSE.sum(axis=1)

            SSE_transform=SSE.div(GENE_var+0.5)[vargenes].sum(axis=1)
            logitify=np.divide(1.0,1.0+np.exp(SSE_transform))#sum))

            df_SSE.append(SSE_sum)
            df_logit.append(logitify)

            X_adjust[curcov].loc[curcells]=logitify

    return X_adjust

def _train_lm(X_df,
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

def compute_X_y(adata_here,
             model_name='lm',
            include_expression=True,
            y_obs=[],
            perturbations_list=None,
            perturbations_obs='guide',
            covariates_list=[],
            control_names=[],
            use_raw=False,
            keep_unassigned=False):
    
    #===========================
    #setup list of perturbations
    #===========================                                                                    
    if perturbations_list is None:
        try:
            assert perturbations_obs in adata_here.obs.columns
        except AssertionError:
            print('ERROR: "'+perturbations_obs+'" is not in adata.obs. It is needed for determining the perturbations present in the data.')
            return
        perturbations_list=_get_perturbations(adata_here,
                                             perturbations_obs=perturbations_obs)

    #get the subset of perturbations from the list that are in adata
    perturbations_list=perturb_overlap_obs(perturbations_list,adata_here,'perturbations')
    
    #===
    # X
    #===
    X_df=obs_to_design_matrix(adata_here,perturbations_list)
    #set control cells to 0                                                                                                                                                                                                        
    ##control_cells_in_data=list(set(control_cells).intersection(set(adata_here.obs_names)))
    control_vars=[]
    for control_var in control_names:
        if control_var not in adata_here.obs.columns:
            print('WARNING: control variable "'+control_var+'" not in dataset. Ignoring.')
        else:
            control_vars.append(control_var)
    if len(control_vars)>0:
        control_cells_in_data=list(adata_here[adata_here.obs[control_vars].astype(float).sum(axis=1)>0,:].obs_names)
        if len(control_cells_in_data)>0:
            X_df.loc[control_cells_in_data,:]=0
    else:
        control_cells_in_data=[]

    #===
    # y 
    #===
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

    #==========
    #covariates
    #==========
    covariates_list=perturb_overlap_obs(covariates_list,adata_here,'covariates')
    covariates_df=obs_to_design_matrix(adata_here,covariates_list,binarize=False,covariate=True)

    #whether to keep or remove unassigned      
    assigned_cells=set(list(X_df.index[X_df.sum(axis=1)>0])) #with the current perturbations!
    if 'unassigned' not in adata_here.obs.columns:
        print('WARNING: unassigned cells are not annotated and will be ignored')
    else:
        unassigned_cells=adata_here[adata_here.obs['unassigned']>0,:].obs_names
    keep=list(set(control_cells_in_data).union(assigned_cells))
    if keep_unassigned and 'unassigned' in adata_here.obs.columns:
        keep=list(set(control_cells_in_data).union(assigned_cells).union(unassigned_cells))
    
    X_df=X_df.loc[keep,:]
    y_df=y_df.loc[keep,:]
    covariates_df=covariates_df.loc[keep,:]
    print(X_df.shape,covariates_df.shape,y_df.shape)
    
    return(X_df,y_df,covariates_df)

def split_train_valid_test(adata_here,
                           training_proportion=0.6,
                           validation_proportion=0.2,
                           test_proportion=0.2,
                           rng=None,copy_adata=False):
    
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



def train_lm(adata_here,lm,
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
        control_names=[],
        copy_adata=False):
    
    if copy_adata: adata_here = adata_here.copy()
        
    #split into train/test unless already done
    #if 'train_valid_test' not in adata_here.obs:
    if True:
        print('WARNING: Over-writing adata.obs["'+'PS.train_valid_test"]')
        adata_here.obs['PS.train_valid_test']=split_train_valid_test(adata_here,
                               training_proportion=training_proportion,
                               validation_proportion=validation_proportion,
                               test_proportion=test_proportion,
                               rng=rng)
        
    
    #first, make datasets
    X_df,y_df,covariates_df=compute_X_y(adata_here,
                                model_name=model_name,
                                include_expression=include_expression,
                                y_obs=y_obs,
                               perturbations_list=perturbations_list,
                               covariates_list=covariates_list,
                               use_raw=use_raw,
                                        control_names=control_names,
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
        
    #training set
    X_cells=list(X_df.index)
    print(adata_here.obs['PS.train_valid_test'].value_counts())
    training=[cellidx for cellidx in range(X_df.shape[0]) if adata_here.obs['PS.train_valid_test'].loc[X_cells[cellidx]]=='train']
    validation=[cellidx for cellidx in range(X_df.shape[0]) if adata_here.obs['PS.train_valid_test'].loc[X_cells[cellidx]]=='valid']
    test=[cellidx for cellidx in range(X_df.shape[0]) if adata_here.obs['PS.train_valid_test'].loc[X_cells[cellidx]]=='test']
    
    print('train',len(training))
    print('valid',len(validation))
    print('test',len(test))
    
    #train model
    print('\nFitting model\n',lm)
    coef,trained_lm=_train_lm(pd.DataFrame(X_plus_covariates,
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
        X_adjust=pd.DataFrame(np.array(bayes_cov_col(y_df,
                                                pd.DataFrame(X_plus_covariates,
                                                            index=X_df.index),
                                                adjust_vars_idx,
                                                trained_lm)),
                              index=X_df.index,
                              columns=X_plus_covariates_columns)
        print('Re-fitting model on adjusted data\n',lm)
        coef,trained_lm=_train_lm(X_adjust,
                                 y_df,
                                 lm,
                                 training=training) 
        
    #compute the performance of the model
    ##y_pred=trained_lm.predict(np.array(X_df))
    ##performance=compute_performance(y_df,y_pred,training,validation,test)
        
        
    #save items into adata
    if 'PS.'+model_name+'.X' in adata_here.uns:
        print('WARNING: Over-writing adata.uns["'+'PS.'+model_name+'.X"]')
    adata_here.uns['PS.'+model_name+'.X']=X_df
    if 'PS.'+model_name+'.y' in adata_here.uns:
        print('WARNING: Over-writing adata.uns["'+'PS.'+model_name+'.y"]')
    adata_here.uns['PS.'+model_name+'.y']=y_df
    if 'PS.'+model_name+'.covariates' in adata_here.uns:
        print('WARNING: Over-writing adata.uns["'+'PS.'+model_name+'.covariates"]')
    adata_here.uns['PS.'+model_name+'.covariates']=covariates_df
    if 'PS.'+model_name+'.coef' in adata_here.uns:
        print('WARNING: Over-writing adata.uns["'+'PS.'+model_name+'.coef"]')
    adata_here.uns['PS.'+model_name+'.coef']=coef
    ##if 'PS.'+model_name+'.performance' in adata_here.uns:
    ##    print('WARNING: Over-writing adata.uns["'+'PS.'+model_name+'.performance"]')
    ##adata_here.uns['PS.'+model_name+'.performance']=performance

    if adjust:
        if 'PS.'+model_name+'.X' in adata_here.uns:
            print('WARNING: Over-writing adata.uns["'+'PS.'+model_name+'.X_adjust"]')
        adata_here.uns[model_name+'.X_adjust']=pd.DataFrame(X_adjust.loc[X_df.index,X_df.columns],
                                                index=X_df.index,
                                                columns=X_plus_covariates_columns)
    
    if copy_adata:
        return(adata_here)
    

