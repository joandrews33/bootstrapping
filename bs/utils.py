import numpy as np
import pandas as pd

def get_count_list(df,feature):

    """Creates a list of tuples counting the value and the abundance of values for a categorical feature."""

    return list(zip(df[feature].value_counts().index,df[feature].value_counts().values))

def trim_count_list(guide_list,min_counts):

    """ A recursive binary search algorithm for removing guides from the guidelist with fewer than min_counts counts. 
    Assumes the guidelist is sorted in descending abundance."""

    if len(guide_list) < 3:
        
        return [(g,c) for g,c in guide_list if c>=min_counts]

    a = int(np.floor(len(guide_list)/2))
    v = guide_list[a][1]
    if v >= min_counts:
        return guide_list[:a+1] + trim_count_list(guide_list[a+1:],min_counts)
    else:
        return trim_count_list(guide_list[:a+1],min_counts)
    
def drop_df_duplicates(input_file, output_file, drop_subsets=[['plate','well','tile','cell_0'],['plate','well','site','cell_1']]):
    df=pd.read_csv(input_file)
    for subset in drop_subsets:
        df=df.drop_duplicates(subset=subset)
    df.to_csv(output_file)

def get_gene_guide_depths(df,gene_feature='gene_symbol',guide_feature='cell_barcode_0',mode='volcano',
                          min_counts=None, target_genes=None,blacklist_genes=['nontargeting'],
                          DropDuplicates=False,drop_subsets=[['plate','well','tile','cell_0'],['plate','well','site','cell_1']]):

    # {pop_0}_+{pop_1}_+{pop_2}-N{count}+{pop_2}-N{count}
    # {pop}_{count-string}
    # {count-string} = {pop}_{count-string}+{pop}_{count-string}

    # {pop}_{sub_pop}+{sub_pop} #How can I do this unambiguously? Bottom up?

    # {pop}-N{count} : base
    # {pop}_{count_strings}+{count_strings}:

    if DropDuplicates:
        for subset in drop_subsets:
            df = df.dropduplicates(subset=subset)

    if target_genes is not None and ('all' not in target_genes):
        df = df.query(gene_feature+' in @target_genes')

    df = df.query(gene_feature+' not in @blacklist_genes')

    if mode=='volcano':
        df_barcode_counts=(df[[gene_feature,guide_feature]]
                       .groupby([gene_feature,guide_feature]).apply(lambda x: len(x))
                       .reset_index().rename(columns={0:'count'})
                       #.set_index([gene_feature,guide_feature])
                       )
    elif mode=='power':
        df_barcode_counts=(df[[gene_feature,guide_feature]]
                       .groupby([gene_feature,guide_feature]).apply(lambda x: '{power_levels}')
                       .reset_index().rename(columns={0:'count'})
                       #.set_index([gene_feature,guide_feature])
                       )
    
    if min_counts is not None:
        df_barcode_counts.query('count>@min_counts')

    return df_barcode_counts[gene_feature], df_barcode_counts[guide_feature], df_barcode_counts['count']

def get_guide_depths():
    return None