import numpy as np
import pandas as pd
from random import choices

def grouped_standardization(df,population_feature='gene_symbol',control_value='nontargeting',group_columns=['plate','well'], index_columns = ['tile','cell_0'],cat_columns=['gene_symbol','cell_barcode_0'],target_features=None,drop_features=False):
    '''Standardizes the numerical columns of df by evaluating the robust z-score. The null model for each
    measurement is estimated as its empirical distribution for the null_gene. If group_column is specified, the 
    null model is evaluated separately for each category in group_column. (E.g., standardizing by well.)'''

    #Warning, this will fail is dataframe contains repeated values for cells

    df_out = df.copy().drop_duplicates(subset=group_columns+index_columns)

    if target_features is None:
        target_features = [col for col in df.columns if col not in group_columns+index_columns+cat_columns]
    
    if drop_features:
        df = df[group_columns+index_columns+cat_columns+target_features]

    unstandardized_features = [col for col in df.columns if col not in target_features]
    
    group_medians = df.query(population_feature+'==@control_value').groupby(group_columns)[target_features].median()
    group_mads = df.query(population_feature+'==@control_value').groupby(group_columns)[target_features].mad() # THis is not actually a mad. Grrr.

    df_out = pd.concat([df_out[unstandardized_features].set_index(group_columns+index_columns),df_out.set_index(group_columns+index_columns)[target_features].subtract(group_medians).divide(group_mads).multiply(0.6745)],axis=1)

    return df_out.reset_index()

def split_by_population(df,population_values,population_feature='cell_barcode_0'):
    return [df.query(population_feature+'==@population') for population in population_values]

def split_null(df,control_feature='gene_symbol',control_value='nontargeting'):
    return df.query(control_feature+'==@control_value')

def compute_variate(df,variate='median',target_features=None,
                    index_features=['gene_symbol','cell_barcode_0'],
                    n_samples=None,sample_depth=None,
                    min_count=None,group_feature=None,wildcards=None):
    '''

    n_samples: Number of times to sample variate per index_feature group. If index_feature group contains sgRNA barcodes,
    e.g., for null model sampling, it will create n_samples PER GUIDE. YOU DON"T NEED TO DO 100k samples!

    n_depth: Number of cells to sample for variate computation.
    min_count: determines the minimum number of cells per index_group. If a group has fewer cells, it is dropped from the data-frame.  
    default=None 

    Modes:
    Normal
    Bootstrapping
    Subsampling
    n_samples = None => Normal Mode
    n_samples not None and n_depth = None => Bootstrap at full representation
    n_samples and n_depth not None => Subsample with replacement for variate.

    Note: for the null samples, it will generate N_guides x n_samples total samples. Can reduce n_bootstrap from 1E5 to 1E5/N_guides
    '''

    if target_features is None:
        target_features = [col for col in df.columns if col not in index_features]

    if variate=='median':
        if n_samples is None: #Normal-mode
            df_out = df.groupby(index_features)[target_features].median().reset_index()
            if 'population_count' in df.columns:
                df_out['population_count']=df.groupby(index_features)['population_count'].sum().reset_index().drop(columns=index_features)
            else:
                df_out['population_count']=df.groupby(index_features).apply(lambda x: len(x)).reset_index().drop(columns=index_features)
            if min_count is not None:
                df_out = df_out.query('population_count >= @min_count')
            return df_out
        
        else: #Bootstrapping mode
            if sample_depth is not None: 
                df_out = pd.concat([df.groupby(index_features).sample(sample_depth,replace=True)
                                    .groupby(index_features)[target_features].median()
                                    .reset_index() for _ in range(n_samples)
                                    ],axis=0).assign(population_count=sample_depth)
                return df_out
            else:
                raise(NotImplementedError('Bootstrapping without subsampling is not yet implemented'))

    if variate=='auc':
        #The tricky part here is that I need to compute and pass in a value for the null auc. Where to do this?
        raise NotImplementedError('Have not implemented deltaAUC as a variate yet.')
    
def add_pseudogenes(df,gene_feature='gene_symbol',guide_feature='cell_barcode_0',pseudogene_population='nontargeting',sample_depth=None,guides_per_pseudogene=4):
    ''' Duplicates the non-targeting cells and groups them into "pseudo-genes" to serve as a non-targeting control for gene level volcano plots.
    '''
    guide_list = df.query(gene_feature+'==@pseudogene_population')[guide_feature].unique()
    num_guides = len(guide_list)
    gene_symbol = []
    for count in range(int(num_guides/guides_per_pseudogene+1)):
        gene_symbol.append(guides_per_pseudogene*[f'pseudogene{count}'])
    gene_symbol = np.concatenate(gene_symbol)
    
    df_pseudo = df.query(gene_feature+'==@pseudogene_population')
    conversion_dict = {guide:gene for guide,gene in zip(guide_list,gene_symbol[:num_guides])}
    df_pseudo[gene_feature] = df_pseudo[guide_feature].apply(lambda x: conversion_dict[x])
    return pd.concat([df,df_pseudo],axis=0)

def merge_composite_variate_bootstrap(df_list,how='median',wildcards=None,drop_columns=[]):
    from random import sample
    population = None
    if wildcards is not None:
        if 'population' in wildcards.keys():
            population = wildcards.population

    column_names = df_list[0].drop(columns=drop_columns).columns

    return pd.DataFrame(np.median([sample(df.drop(columns=drop_columns).values.tolist(),k=len(df)) for df in df_list],axis=0),columns=column_names).assign(population=population).assign(population_count=wildcards.depth_string)
    #return pd.DataFrame(np.median([df.drop(columns=drop_columns).values for df in df_list],axis=0),columns=column_names).assign(population=population).assign(population_count=wildcards.depth_string)

def evaluate_p_value(df_exp,df_boot,two_sided=True,wildcards=None,population_column='population'):

    df_exp = df_exp.drop(columns=population_column).drop(columns='population_count')#'population')
    df_boot = df_boot.drop(columns='population').drop(columns='population_count')

    population = None
    if wildcards is not None:
        if 'population' in wildcards.keys():
            population = wildcards.population

    def _get_p(series):
        if two_sided:
            return 2*pd.concat([(df_boot>series).sum(),(df_boot<series).sum()],axis=1).transpose().min()/len(df_boot)
        else:
            return (df_boot>series).sum().transpose()/len(df_boot)
        
    return pd.concat([_get_p(df) for _,df in df_exp.iterrows()],axis=1).transpose().rename(columns=lambda x: 'p_'+x).assign(population=population)

def collect_volcano(df_list_var,df_list_pval,left_population_features=['population'],right_population_features=['population']):
    '''Used for collecting multiple dataframes, and merging them into a single output. 
    df_list_var: a list of DataFrames with computed variates
    df_list_pval: a list of DataFames with computed p-vals 
    population_features: list of features in the two sets of DataFrames that will be used to align the data.
    '''
    print(len(df_list_var))
    print(len(df_list_pval))

    return pd.concat(df_list_var,axis=0).merge(pd.concat(df_list_pval,axis=0),left_on=left_population_features,right_on=right_population_features)
        