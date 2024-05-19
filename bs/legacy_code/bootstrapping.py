from random import sample
from random import choices
import pandas as pd
import numpy as np
from functools import partial
import re

""" A collection of functions for performing bootstrapping significance analysis of features in a target DataFrame, using a null model DataFrame."""

## Util Functions

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

#### FILE FORMATS
'''
'combined.csv'  <- Input format
'combined.trans.csv' <- Transforming data so that each feature belongs to (-Inf,Inf)
'combined(.trans).std.csv <- (Transformed) data that has been z-standardized against the null population.
'{population}_N{depth}.var.csv' <- Variates computed from cells of {population} sampled at {depth}. Population can be e.g., null, ACTB, or GATTACA
'{genesymbol}_N{cellbarcode}{depth}-{cellbarcode}{depth}-{cellbarcode}{depth}-{cellbarcode}{depth}.group_var.csv'
'{population}_N{depth-group}.group_var.csv' <-Grouped variates. Typically medians of guide variates for a given gene. {depth-group} is an ordered list of sampling depths for subpopulations, e.g., 176-101-23-5
'{population}_N{depth-group}.pvals.csv' <- pvals and variates for given population computed against null_N{depth-group}.(group_)var.csv
'volcano.csv' <- Aggregation of all subpopulation pvals into a single dataframe. false-discovery-rate transformed pvals. Typically Benjameni-Hochberg
'''

#I don't like these function names. Come up with something better...

def get_target_file_strings(df,population_features=['gene_symbol','cell_barcode_0'],mode='volcano',target_genes=None,blacklist_genes='nontargeting'):
    '''target_genes and blacklist_genes are lists of gene_symbols that should be included/excluded. If target_genes is listed as None, all genes except the blacklist_genes will be included.
    '''

    if not isinstance(population_features,list):
        population_features = [population_features]

    if (target_genes is not None) and (not isinstance(target_genes,list)):
        target_genes = [target_genes]

    if (blacklist_genes is not None) and (not isinstance(blacklist_genes,list)):
        blacklist_genes = [blacklist_genes]
    

    if mode=='volcano':
        df_barcode_counts=(pd.DataFrame(df[population_features].groupby(population_features,as_index=False)
            .apply(lambda x: len(x))).rename(columns={None:'count'})
            .sort_values(by=[population_features[0],'count'],ascending=[True,False])
            #.set_index(population_features)
                       )
    elif mode=='power':
        df_barcode_counts=(df[population_features]
                       .groupby(population_features).apply(lambda x: '{power_levels}')
                       .reset_index().rename(columns={0:'count'})
                       #.set_index(population_features)
                       )
        
    if target_genes is None:
        target_genes = df_barcode_counts[population_features[0]].unique()
    
    if blacklist_genes is not None:
        target_genes = [target for target in target_genes if target not in blacklist_genes]
        
    if len(population_features)==1:
        print(df_barcode_counts[population_features+['count']])
        target_strings = [pop_val+'-N'+str(count) for pop_val,count in df_barcode_counts[population_features+['count']]]
    elif len(population_features)==2:
        target_strings = []
        for pop_val, df_group in df_barcode_counts.groupby(population_features[0]):
            if pop_val in target_genes:
                sub_pops = '+'.join([subpop for subpop in df_group[population_features[1]]])
                count_string = '+'.join([str(val) for val in df_group['count']])
                target_strings.append(pop_val+'_'+sub_pops+'_N'+count_string)
                #target_strings.append(pop_val+'_'+sub_pops+'-N'+count_string)
        
    return target_strings


def drop_df_duplicates(input_file, output_file, drop_subsets=[['plate','well','tile','cell_0'],['plate','well','site','cell_1']]):
    df=pd.read_csv(input_file)
    for subset in drop_subsets:
        df=df.drop_duplicates(subset=subset)
    df.to_csv(output_file)

def get_gene_guide_depths(df,gene_feature='gene_symbol',guide_feature='cell_barcode_0',mode='volcano',
                          min_counts=None,
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

def population_basic_file_format(suffix,directory='process',temp_tags=tuple()):
    file_pattern = f'{directory}/{{population}}.{suffix}'
    if suffix in temp_tags:
        from snakemake.io import temp
        file_pattern = temp(file_pattern)
    return file_pattern

def basic_file_format(suffix,directory='process',temp_tags=tuple()):
    file_pattern = f'{directory}/{suffix}'
    if suffix in temp_tags:
        from snakemake.io import temp
        file_pattern = temp(file_pattern)
    return file_pattern

def population_depth_file_format(suffix,directory='process',temp_tags=tuple()):
    """Format population file pattern, for example:
    population_file_format('pvals.csv') => 'process/{population}_N{depth}.pvals.tif'
    """
    file_pattern = f'{directory}/{{population}}-N{{depth}}.{suffix}'
    if suffix in temp_tags:
        from snakemake.io import temp
        file_pattern = temp(file_pattern)
    return file_pattern

def null_depth_file_format(suffix,directory='process',temp_tags=tuple()):
    """Format population file pattern, for example:
    population_file_format('pvals.csv') => 'process/{population}_N{depth}.pvals.tif'
    """
    file_pattern = f'{directory}/null-N{{depth}}.{suffix}'
    if suffix in temp_tags:
        from snakemake.io import temp
        file_pattern = temp(file_pattern)
    return file_pattern

def composite_population_file_format(suffix,directory='process',temp_tags=tuple()):

    file_pattern = f'{directory}/{{population}}_{{subpop_string}}_N{{depth_string}}.{suffix}'
    if suffix in temp_tags:
        from snakemake.io import temp
        file_pattern = temp(file_pattern)
    return file_pattern

def null_composite_file_format(suffix,directory='process',temp_tags=tuple()):

    file_pattern = f'{directory}/null_N{{depth_string}}.{suffix}'
    if suffix in temp_tags:
        from snakemake.io import temp
        file_pattern = temp(file_pattern)
    return file_pattern

def split_composite_file(input_filename):
    #TODO: Add ability to overwrite the input suffix. e.g., to allow group_var to be replaced with var
    '''
    '{directory}/{{population}}_{{subpopulation-string}}_N{{depth-string}}.{suffix}'
    '{directory}/null_N{{depth-string}}.{suffix}'
    {population}, {subpopulation-string} and {depth-string} cannot contain '.', '_', '+' or '/'.
    '''
    directory = '/'.join(input_filename.split('/')[:-1])
    if directory:
        directory = directory +'/'
    suffix = '.'.join(input_filename.split('/')[-1].split('.')[1:])
    if suffix:
        suffix = '.'+suffix
    filename = input_filename.split('/')[-1].split('.')[0]
    if 'null' in filename:
        pattern = r'null_N(?P<depth_string>[^.]+)'
        variables = re.match(pattern,filename).groupdict()
        depths = variables['depth_string'].split('+')
        return [directory + 'null-N'+depth+ suffix for depth in depths]
    else:
        pattern = r'(?P<parent_pop>.+?)_(?P<subpop_string>.+)_N(?P<depth_string>[^.]+)'
        variables = re.match(pattern,filename).groupdict()
        sub_pops = variables['subpop_string'].split('+')
        depths = variables['depth_string'].split('+')
        return [directory +sub_pop + '-N'+depth+ suffix for sub_pop,depth in zip(sub_pops,depths)]

def split_composite_wildcards(wildcards,suffix,directory):
    '''
    '{directory}/{{population}}_{{subpop_string}}_N{{depth_string}}.{suffix}'
    '''

    directory = directory + '/'

    suffix='.'+suffix
    depths = wildcards.depth_string.split('+')
    if 'subpop_string' in wildcards.keys(): #The composite file is from a targeting gene
        sub_pops = wildcards.subpop_string.split('+')
        return [directory +sub_pop + '-N'+depth+ suffix for sub_pop,depth in zip(sub_pops,depths)]
    else: #the composite file is from the null population
        return [directory + 'null-N'+depth+ suffix for depth in depths]