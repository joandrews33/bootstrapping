import re
import pandas as pd

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
                #sub_pops = '+'.join([subpop for subpop in df_group[population_features[1]]])
                count_string = '+'.join([str(val) for val in df_group['count']])
                target_strings.append(pop_val+'_N'+count_string)#target_strings.append(pop_val+'_'+sub_pops+'_N'+count_string)
                #target_strings.append(pop_val+'_'+sub_pops+'-N'+count_string)
        
    return target_strings


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

    #file_pattern = f'{directory}/{{population}}_{{subpop_string}}_N{{depth_string}}.{suffix}'
    file_pattern = f'{directory}/{{population}}_N{{depth_string}}.{suffix}'
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