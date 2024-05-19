import pandas as pd

def compute_variate(df,variate='median',target_features=None,
                    index_features=['gene_symbol','cell_barcode_0'],
                    n_samples=None,sample_depth=None,
                    min_count=None,group_feature=None,wildcards=None):
    '''

    n_samples: Number of times to sample variate per index_feature group. If index_feature group contains sgRNA barcodes,
    e.g., for null model sampling, it will create n_samples PER GUIDE. YOU DON"T NEED TO DO 100k samples!

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
                df_out['population_count']=df.groupby(index_features)['population_count'].sum().reset_index()
            else:
                df_out['population_count']=df.groupby(index_features).apply(lambda x: len(x)).reset_index()
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
    
def _testerer():
    return None