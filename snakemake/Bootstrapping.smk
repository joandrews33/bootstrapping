### Input
import sys, os
from functools import partial
import pandas as pd

sys.path.append(os.path.dirname(sys.path[0]))
from bs.snakemake_wrapper import snakeify as sfy
import bs.bootstrapping as strap
from bs.filenames import basic_file_format, population_basic_file_format, population_depth_file_format, composite_population_file_format, null_depth_file_format, null_composite_file_format, split_composite_wildcards, get_target_file_strings
from bs.utils import get_gene_guide_depths, get_guide_depths

input_file = config['INPUT_DIRECTORY']+'/'+config['INPUT_FILE']

#GUIDES,GENES = None #Drop GUIDES that don't reach min_counts
min_counts = config['MIN_COUNTS']
target_genes = config['TARGET_GENES']

if 'gene' in config['POPULATION_FEATURE']:
    GENES, GUIDES, DEPTHS = get_gene_guide_depths(pd.read_csv(input_file),min_counts=min_counts,target_genes=target_genes)
elif 'barcode' in config['POPULATION_FEATURE']:
    GUIDES, DEPTHS = get_guide_depths(pd.read_csv(input_file),min_counts=min_counts,target_genes=target_genes)

if 'all' in target_genes:
    TARGET_FILESTRINGS = get_target_file_strings(pd.read_csv(input_file))# ADD Mincounts here.
else:
    TARGET_FILESTRINGS = get_target_file_strings(pd.read_csv(input_file),target_genes=target_genes)

#File Format Definitions
if config['TRANSFORM_FEATURES']:
    std_ext = 'trans.std.csv'
else:
    std_ext = 'std.csv'

basic_file = partial(basic_file_format,directory=config['PROCESS_DIRECTORY'])
population_basic_file = partial(population_basic_file_format,directory=config['PROCESS_DIRECTORY'])
population_depth_file = partial(population_depth_file_format,directory=config['PROCESS_DIRECTORY'])
composite_population_file = partial(composite_population_file_format,directory=config['PROCESS_DIRECTORY'])
null_depth_file = partial(null_depth_file_format,directory=config['PROCESS_DIRECTORY'])
null_composite_file = partial(null_composite_file_format,directory=config['PROCESS_DIRECTORY'])
split_composite_file = partial(split_composite_wildcards,directory=config['PROCESS_DIRECTORY'])

# Optional rule to transform feature values so that they behave better. (e.g., more normal, variance-stabilized, support on the whole real line etc.)
rule all:
    input:
        #basic_file('combined.group_var.csv')
        #population_basic_file('group_var.csv').format(population='TNPO1')
        basic_file('volcano.genes.csv')

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
        sfy(strap.grouped_standardization)(output=output,df=input[0],target_features=config['TARGET_FEATURES'],drop_features=True)

rule standardize_transformed_features:
    input:
        basic_file('combined.trans.csv')
    output:
        basic_file('combined.trans.std.csv')
    run:
        sfy(strap.grouped_standardization)(output=output,df=input[0],target_features=config['TARGET_FEATURES'],drop_features=True)

rule split_null:
    input: 
        basic_file('combined'+'.'+std_ext)
    output:
        temp(basic_file('null'+'.'+std_ext))
    run:
        sfy(strap.split_null)(output=output,df=input[0]) #can add config for null population feature/value if not gene_symbol/nontargeting

rule compute_feature_variates:
    input:
        basic_file('combined'+'.'+std_ext)
    output:
        basic_file('combined.var.csv')
    run:
        sfy(strap.compute_variate)(output=output,df=input[0],
                                   variate=config['VARIATE'],index_features=['gene_symbol','cell_barcode_0'],
                                   target_features=config['TARGET_FEATURES'])#Add min counts here if using. 

rule compute_composite_feature_variates:
    input:
        basic_file('combined.var.csv')
    output:
        basic_file('combined.group_var.csv')
    run:
        sfy(strap.compute_variate)(output=output,df=input[0],
                                   variate='median',index_features=['gene_symbol'],
                                   target_features=config['TARGET_FEATURES'],min_counts=min_counts)

rule split_by_guide:
    input:
        basic_file('combined.var.csv')
    output:
        temp([expand(population_basic_file('var.csv'),population=GUIDES)])
    run:
        sfy(strap.split_by_population)(output=output,df=input[0],population_values=GUIDES)

rule split_by_gene:
    input:
        basic_file('combined.group_var.csv')
    output:
        [expand(population_basic_file('group_var.csv'),population=set(GENES))]
        #temp([expand(population_basic_file('group_var.csv'),population=GENES)])
    run:
        sfy(strap.split_by_population)(output=output,df=input[0],population_values=set(GENES),population_feature='gene_symbol')

### Add steps for bootstrapping null here...
# Add option to either mix all cells, or sample by nt guide.
rule bootstrap_null_variate:
    input:
        basic_file('null'+'.'+std_ext)
    output:
        null_depth_file('var.boot.csv')#temp(null_depth_file('var.boot.csv'))
    run:
        sample_depth =int(wildcards.depth)
        sfy(strap.compute_variate)(output=output,df=input[0],
                                   variate=config['VARIATE'],index_features=['gene_symbol','cell_barcode_0'],
                                   target_features=config['TARGET_FEATURES'],
                                   n_samples=config['BOOTSTRAP_SAMPLES'],sample_depth=sample_depth)

rule bootstrap_null_composite_variates:
    input:
        partial(split_composite_file,suffix='var.boot.csv')#split_output(composite_population_file('group_var.boot.csv')) #Implement this split function
    output:
        composite_population_file('group_var.boot.csv')#null_composite_file('group_var.boot.csv'))
    run:
        sfy(strap.merge_composite_variate_bootstrap)(
            output=output,
            df_list=input,
            wildcards=wildcards,
            drop_columns=['gene_symbol','cell_barcode_0']
        )

rule calculate_significance:
    input:
        population_basic_file('var.csv'),
        null_depth_file('var.boot.csv')
    output:
        temp(population_depth_file('pvals.csv'))
    run:
        sfy(strap.evaluate_p_value)(output=output,
                                   df_exp=input[0],
                                   df_boot=input[1],
                                   two_sided=True,wildcards=wildcards)

rule calculate_composite_significance:
    input:
        population_basic_file('group_var.csv'),
        null_composite_file('group_var.boot.csv') #Define This format
    output:
        composite_population_file('pvals.csv')#temp(composite_population_file('pvals.csv')) #Define This format
    run:
        sfy(strap.evaluate_p_value)(output=output,
                                   df_exp=input[0],
                                   df_boot=input[1],
                                   two_sided=True,wildcards=wildcards,
                                   population_column='gene_symbol')

rule collect_volcano_dataframe:
    input:
        [expand(population_basic_file('group_var.csv'),population=set(GENES))],
        #[expand(config['PROCESS_DIRECTORY'] + '/{filestring}.{extension}',filestring=TARGET_FILESTRINGS,extension=['group_var.csv'])],
        [expand(config['PROCESS_DIRECTORY'] + '/{filestring}.{extension}',filestring=TARGET_FILESTRINGS,extension=['pvals.csv'])] #Define target filestrings
    output:
        basic_file('volcano.genes.csv')
    run:
        n_input_files = len(input)
        variate_files = input[:int(n_input_files/2)]
        p_files = input[int(n_input_files/2):]
        sfy(strap.collect_volcano)(output=output,
                                   df_list_var=variate_files,
                                   df_list_pval=p_files,
                                   left_population_features='gene_symbol'
                                   )