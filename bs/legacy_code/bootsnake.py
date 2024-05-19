### Theres a snake in my boots!
import inspect
from random import sample
from random import choices
import pandas as pd
import numpy as np
from bs.watersnake import Ouroboros

class BootSnake(Ouroboros):

    @staticmethod
    def _grouped_standardization(df,population_feature='gene_symbol',control_value='nontargeting',group_columns=['plate','well'], index_columns = ['tile','cell_0'],cat_columns=['gene_symbol','cell_barcode_0'],target_features=None,drop_features=False):
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
        group_mads = df.query(population_feature+'==@control_value').groupby(group_columns)[target_features].mad()

        df_out = pd.concat([df_out[unstandardized_features].set_index(group_columns+index_columns),df_out.set_index(group_columns+index_columns)[target_features].subtract(group_medians).divide(group_mads).multiply(0.6745)],axis=1)

        return df_out.reset_index()
    
    @staticmethod
    def _split_by_population(df,population_values,population_feature='cell_barcode_0'):
        return [df.query(population_feature+'==@population') for population in population_values]
    
    @staticmethod
    def _split_null(df,control_feature='gene_symbol',control_value='nontargeting'):
        return df.query(control_feature+'==@control_value')
    
    @staticmethod
    def _compute_variate(df,variate='median',target_features=None,n_samples=None,sample_depth=None,index_features=['gene_symbol','cell_barcode_0'],group_feature=None,wildcards=None):
        # What does this retrun for n_samples=None? How does a pd.series save?
        # Implement auc
        # WARNING: I don't think this is evaluating correctly. does the groupby sample mean that I am averaging MORE than the target number of cells whenever there is more than one guide present in df??
        
        population = None
        if wildcards is not None:
            if 'population' in wildcards.keys():
                population = wildcards.population

        if target_features is None:
            target_features = [col for col in df.columns if col not in index_features]
        if variate=='median':
            if n_samples is None:
                return df.groupby(index_features)[target_features].median().assign(population=population)
            else:
                if group_feature is not None:
                    populations = choices(df[group_feature].unique(), k=n_samples)
                    return pd.concat([df.query(group_feature+'==@target_population')[target_features]
                                    .sample(n=sample_depth,replace=True).median() 
                                    for target_population in populations]
                                    ,axis=1).transpose().assign(population=population)
                else:
                    return pd.concat([df[target_features]
                                      .sample(n=sample_depth,replace=True).median() 
                                      for _ in range(n_samples)]
                                      ,axis=1).transpose().assign(population=population)
                    
        if variate=='auc':
            #The tricky part here is that I need to compute and pass in a value for the null auc. Where to do this?
            raise NotImplementedError('Have not implemented deltaAUC as a variate yet.')

    @staticmethod
    def _evaluate_composite_variate(df_list,how='median',wildcards=None):
        population = None
        if wildcards is not None:
            if 'population' in wildcards.keys():
                population = wildcards.population

        column_names = df_list[0].columns
        return pd.DataFrame(np.median([df.values for df in df_list],axis=0),columns=column_names).assign(population=population)

    @staticmethod
    def _evaluate_p_value(df_exp,df_boot,two_sided=True,wildcards=None):

        df_exp = df_exp.drop(columns='population')
        df_boot = df_boot.drop(columns='population')

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

    @staticmethod
    def _collect_volcano(df_list_var,df_list_pval,population_features=['population']):
        return pd.concat(df_list_var,axis=0).merge(pd.concat(df_list_pval,axis=0),left_on=population_features,right_on=population_features)
         

    ### Class Utility Methods
    @staticmethod
    def load_methods():
        methods = inspect.getmembers(BootSnake)
        for name, f in methods:
            if name not in ('__doc__', '__module__') and name.startswith('_'):
                #Ouroboros.add_method('Snake', name[1:], Ouroboros.call_from_snakemake(f))
                BootSnake.add_method('BootSnake', name[1:], BootSnake.call_from_snakemake(f))

    @staticmethod
    def add_method(class_, name, f):
        f = staticmethod(f)
        exec('%s.%s = f' % (class_, name))
    
BootSnake.load_methods()