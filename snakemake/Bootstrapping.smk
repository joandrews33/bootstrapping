from bs.bootstrapping import get_gene_guide_depths, get_guide_depths, drop_df_duplicates, get_target_file_strings
from bs.bootstrapping import basic_file_format, population_basic_file_format, population_depth_file_format, null_depth_file_format, composite_population_file_format, null_composite_file_format, split_composite_file, split_composite_wildcards
from os.path import isfile 
import pandas as pd
from bs.bootsnake import BootSnake
from functools import partial

input_file = config['INPUT_DIRECTORY']+'/'+config['INPUT_FILE']

### In case the input file hasn't already removed duplicated cells from the genotype-phenotype alignment, do it now.
if config['DROP_DUPLICATES']:
    output_file = config['PROCESS_DIRECTORY']+'/'+config['INPUT_FILE']
    if not isfile(output_file):
        drop_df_duplicates(input_file,output_file)
    input_file = output_file

### Extract the major 

if 'gene' in config['POPULATION_FEATURE']:
    GENES, GUIDES, DEPTHS = get_gene_guide_depths(pd.read_csv(input_file))
elif 'barcode' in config['POPULATION_FEATURE']:
    GUIDES, DEPTHS = get_guide_depths(pd.read_csv(input_file))

#GUIDES=['ACAACCCCA','AGTACGACT','CAGCACCGT','CCTGGAGGA']

if config['TRANSFORM_FEATURES']:
    std_ext = 'trans.std.csv'
else:
    std_ext = 'std.csv'

### FILE FORMATS
    
#basic_file {directory}/{suffix}
#population_basic_file {directory}/{{population}}.{suffix}
#population_depth_file {directory}/{{population}}-N{{depth}}.{suffix}
#null_depth_file {directory}/null-N{{depth}}.{suffix}
#composite_population_file {directory}/{{parent_pop}}_{{subpop-string}}_N{{depth-string}}.{suffix}
#null_composite_file {directory}/null_N{{depth-string}}.{suffix}

basic_file = partial(basic_file_format,directory=config['PROCESS_DIRECTORY'])
population_basic_file = partial(population_basic_file_format,directory=config['PROCESS_DIRECTORY'])
population_depth_file = partial(population_depth_file_format,directory=config['PROCESS_DIRECTORY'])
null_depth_file = partial(null_depth_file_format,directory=config['PROCESS_DIRECTORY'])
composite_population_file = partial(composite_population_file_format,directory=config['PROCESS_DIRECTORY'])
null_composite_file = partial(null_composite_file_format,directory=config['PROCESS_DIRECTORY'])
split_composite_file = partial(split_composite_wildcards,directory=config['PROCESS_DIRECTORY'])

#split_output = partial(split_composite_file_format,directory=config['PROCESS_DIRECTORY'],suffix='var.csv')

### Specifying which files to produce

if 'all' in config['TARGET_GENES']:
    TARGET_FILESTRINGS = get_target_file_strings(pd.read_csv(input_file))
else:
    TARGET_FILESTRINGS = get_target_file_strings(pd.read_csv(input_file),target_genes=config['TARGET_GENES'])


#target_files = [config['PROCESS_DIRECTORY'] + '/' + filestring + '.pvals.csv' for filestring in TARGET_FILESTRINGS]

rule all:
    input:
        basic_file('volcano.csv')
        #target_files

# Optional rule to transform feature values so that they behave better. (e.g., more normal, variance-stabilized, support on the whole real line etc.)
rule transform_features:
    input:
        basic_file('combined.csv')
    output:
        basic_file('combined.trans.csv')
    run:
        print('"transform_features" not yet implemented.')

rule standardize_features:
    input:
        basic_file('combined.csv')
    output:
        basic_file('combined.std.csv')
    run:
        BootSnake.grouped_standardization(output=output,df=input[0],target_features=config['TARGET_FEATURES'],drop_features=True)

rule standardize_transformed_features:
    input:
        basic_file('combined.trans.csv')
    output:
        basic_file('combined.trans.std.csv')
    run:
        BootSnake.grouped_standardization(output=output,df=input[0],target_features=config['TARGET_FEATURES'],drop_features=True)


#Making splits of the input dataframe so that it doesn't have to load large files full of irrelevant data for every bootstrapping.
rule split_null:
    input: 
        basic_file('combined'+'.'+std_ext)
    output:
        temp(basic_file('null'+'.'+std_ext))
    run:
        BootSnake.split_null(output=output,df=input[0]) #can add config for null population feature/value if not gene_symbol/nontargeting

rule split_by_guides:
    input: 
        basic_file('combined'+'.'+std_ext)
    output:
        temp([expand(population_basic_file(std_ext),population=GUIDES)]) #Can I do this all in one go, or do I have to do one job per GUIDE?
    run:
        BootSnake.split_by_population(output=output,df=input[0],population_values=GUIDES)#,population_feature=config['POPULATION_FEATURE'])

rule compute_feature_variates:
    input:
        population_basic_file(std_ext)
    output:
        temp(population_depth_file('var.csv'))
    run:
        if "null" in wildcards.population:
            BootSnake.compute_variate(output=output,df=input[0],variate=config['VARIATE'],target_features=config['TARGET_FEATURES'])#group_feature='cell_barcode_0' Should add a grouping feature here. #,wildcards=wildcards)
        else:
            BootSnake.compute_variate(output=output,df=input[0],variate=config['VARIATE'],target_features=config['TARGET_FEATURES'])#,wildcards=wildcards)

rule bootstrap_feature_variates:
    input:
        population_basic_file(std_ext)
    output:
        temp(population_depth_file('var.boot.csv'))
    run:
        if "null" in wildcards.population:
            print('Warning: Bootstrapping is mixing cells from different null guides together.')
        sample_depth=int(wildcards.depth)#'DEPTHWILDCARD' #How do I make this work?
        BootSnake.compute_variate(output=output,df=input[0],variate=config['VARIATE'],target_features=config['TARGET_FEATURES'],n_samples=config['BOOTSTRAP_SAMPLES'],sample_depth=sample_depth)#,wildcards=wildcards)

rule compute_composite_variates:
    input:
        partial(split_composite_file,suffix='var.csv')#split_output(composite_population_file('group_var.csv')) #Implement this split function ops.bootstrapping.split_composite_file almost works, except that the file tag needs to change from group_var to var.csv. At what point will the wildcards get evaluated to literals? Is this process here even theoretically possible?
    output:
        temp(composite_population_file('group_var.csv'))
    run:
        BootSnake.evaluate_composite_variate(
            output=output,
            df_list=input,
            wildcards=wildcards
        )

rule bootstrap_composite_variates:
    input:
        partial(split_composite_file,suffix='var.boot.csv')#split_output(composite_population_file('group_var.boot.csv')) #Implement this split function
    output:
        temp(composite_population_file('group_var.boot.csv'))
    run:
        print('"bootstrap_composite_variates" not yet implemented')

rule bootstrap_null_composite_variates:
    input:
        partial(split_composite_file,suffix='var.boot.csv')#split_output(composite_population_file('group_var.boot.csv')) #Implement this split function
    output:
        temp(null_composite_file('group_var.boot.csv'))
    run:
        BootSnake.evaluate_composite_variate(
            output=output,
            df_list=input,
            wildcards=wildcards
        )

rule calculate_significance:
    input:
        population_depth_file('var.csv'),
        null_depth_file('var.boot.csv')
    output:
        temp(population_depth_file('pvals.csv'))
    run:
        BootSnake.evaluate_p_value(output=output,
                                   df_exp=input[0],
                                   df_boot=input[1],
                                   two_sided=True,wildcards=wildcards)

rule calculate_composite_significance:
    input:
        composite_population_file('group_var.csv'),
        null_composite_file('group_var.boot.csv')
    output:
        temp(composite_population_file('pvals.csv'))
    run:
        BootSnake.evaluate_p_value(output=output,
                                   df_exp=input[0],
                                   df_boot=input[1],
                                   two_sided=True,wildcards=wildcards)

rule sample_significance:
    input:
        population_depth_file('var.boot.csv'),
        null_depth_file('var.boot.csv')
    output:
        temp(population_depth_file('pvals.power.csv'))
    run:
        BootSnake.evaluate_p_value(output=output,
                                   df_exp=input[0],
                                   df_boot=input[1],
                                   two_sided=True,wildcards=wildcards)

rule sample_composite_significance:
    input:
        composite_population_file('group_var.boot.csv'),
        null_composite_file('group_var.boot.csv')
    output:
        temp(composite_population_file('pvals.power.csv'))
    run:
        BootSnake.evaluate_p_value(output=output,
                                   df_exp=input[0],
                                   df_boot=input[1],
                                   two_sided=True)

rule collect_volcano_dataframe:
    input:
        [expand(config['PROCESS_DIRECTORY'] + '/{filestring}.{extension}',filestring=TARGET_FILESTRINGS,extension=['group_var.csv'])],
        [expand(config['PROCESS_DIRECTORY'] + '/{filestring}.{extension}',filestring=TARGET_FILESTRINGS,extension=['pvals.csv'])]
    output:
        basic_file('volcano.csv')
    run:
        n_input_files = len(input)
        print(n_input_files)
        variate_files = input[:int(n_input_files/2)]
        p_files = input[int(n_input_files/2):]
        print('here')
        print(variate_files)
        BootSnake.collect_volcano(output=output,
                                   df_list_var=variate_files,
                                   df_list_pval=p_files
                                   )

