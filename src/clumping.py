# Authors: John Rouhana, Jiali Wang
# Segre lab, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA
# Date: April 2021

import numpy as np
import pandas as pd

import argparse
import os

import time
import re

def str2bool(s):
    if str2bool=='True':
        return True
    elif str2bool=='False':
        return False
    else:
        print("Unexpected input")

def main():
    parser = argparse.ArgumentParser(description='GWAS clumping for single cell analysis')
    parser.add_argument('--proxy_table', type=str, required=True,
                    help='Path to GWASVar2Gene proxy table')
    parser.add_argument('--traits_table', type=str, required=True,
                    help='Path to GWASVar2Gene traits table')
    parser.add_argument('--sQTL_table', type=str, required=True,
                    help='Path to GWASVar2Gene sQTL_table')
    parser.add_argument('--eQTL_table', type=str, required=True,
                    help='Path to GWASVar2Gene eQTL_table')
    parser.add_argument('--opentarget', type=str, default=None,
                    help='Path to open targets file')
    parser.add_argument('-ancestry', type=str, default=None,
                    help='Case insensitive string regex indicating what ancestry to look for in GWAS. All GWAS not including string are excluded.')
    parser.add_argument('--significant_traits', type=str, required=True,
                    help='Tab delimited text file, with the first column the generic names of the traits of interest and the following columns the trait names that map to the generic names. Trait names (2-n columns) should match the trait column in traits_table exactly.')
    parser.add_argument('-treat_sig_traits_as_one', type=str, default=None,
                    help='If a string is given, all significant traits will be treated as a single trait under this shared trait name.')
    parser.add_argument('-out_directory', type=str, default='',
                    help='Directory where output files will be written. Defaults to working directory.')

    args = parser.parse_args()


    #Retrieve user-defined parameters
    proxy_file  = args.proxy_table
    sQTL_tbl=args.sQTL_table
    eQTL_tbl=args.eQTL_table
    traits_file = args.traits_table
    open_target=args.opentarget
    sig_traits_in = args.significant_traits
    ancestry = args.ancestry
    out_dir = args.out_directory
    tst = args.treat_sig_traits_as_one

    #add / if necessary
    if (out_dir != '') and (not out_dir.endswith('/')):
        out_dir = out_dir+'/'
    
    #create output folder if it doesn't exist
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)    


    #read in variants, split by significant/null first
#     variants_df = pd.read_csv(traits_file, sep='\t', dtype={'pubmed_id':str}, low_memory=False)
    variants_df = pd.read_csv(traits_file, sep='\t', low_memory=False)
    variants_df['GWAS_p_value']=variants_df['GWAS_p_value'].astype(float)
    
    #Verify genome wide significance
    variants_df = variants_df.loc[variants_df.GWAS_p_value < 5e-8].copy()

    #read in traits to use
    with open(sig_traits_in) as f:
        sig_traits = [x.split('\t') for x in f.read().splitlines()]
       
    #Adjust GWAS names in table to reflect desired annotations
    for trait_list in sig_traits:
        for trait in trait_list[1:]:
            variants_df.loc[variants_df.trait == trait, 'trait'] = trait_list[0]
    
    sig_traits = [x[0] for x in sig_traits]

    #Narrow variants_df now to only include ancestry of interest
    if ancestry is not None:
        variants_df = variants_df.loc[variants_df['Sample_Ancestries'].notna()].copy()
        variants_df = variants_df.loc[variants_df['Sample_Ancestries'].apply(lambda x: re_match(x, ancestry))].copy() 
        variants_df.reset_index(drop=True, inplace=True)

    #Separate by significant/not significant traits
    sig_variants_df = variants_df.loc[variants_df.trait.isin(sig_traits)].copy()
    remainder_variants_df = variants_df.loc[~(variants_df.index.isin(sig_variants_df.index))].copy()
    #Remove significant GWAS variants from null set
    null_variants_df = remainder_variants_df.loc[~(remainder_variants_df.GWAS_variant.isin(sig_variants_df.GWAS_variant))].copy()

    #Refresh indices
    sig_variants_df.reset_index(drop=True, inplace=True)
    remainder_variants_df.reset_index(drop=True, inplace=True)
    null_variants_df.reset_index(drop=True, inplace=True)

    #read in GWASVar2Gene table for LD data
    proxy_df = pd.read_csv(proxy_file, sep='\t', usecols=['GWAS_variant', 'LD_variant'])
    #read in eGene file
    tmp = []
    
    df=pd.read_csv(sQTL_tbl,sep='\t',usecols=['sVariant','gene_id','gene_name','GWAS_variant']).drop_duplicates()
    df = df.rename(columns={'sVariant': 'eVariant'})
    tmp.append(df)
    
    df=pd.read_csv(eQTL_tbl,sep='\t',usecols=['eVariant','gene_id','gene_name','GWAS_variant']).drop_duplicates()
    tmp.append(df)
           
    eGene_df = pd.concat(tmp)
    eGene_df = eGene_df.drop_duplicates().copy()
    eGene_df.reset_index(drop=True, inplace=True)
        
    sub_eGene_df = eGene_df[['GWAS_variant', 'gene_name']].drop_duplicates().copy()

    ## Add open target genes
    if open_target is not None:
        open_target_df=pd.read_csv(open_target, dtype=str,usecols=['variant_id', 'bestLocus2Genes'])[['variant_id','bestLocus2Genes']].rename(columns={'variant_id':'GWAS_variant','bestLocus2Genes':'gene_name'}).drop_duplicates().dropna()
        open_target_df['GWAS_variant'] = 'chr'+open_target_df['GWAS_variant']+'_b38'
        sub_eGene_df=pd.concat([sub_eGene_df,open_target_df]).drop_duplicates().dropna()
    

    #Run clumping for significant and null
    print("Clumping significant traits...")
    sig_sub_proxy_clumps, sig_sub_eGene_clumps, sig_all_clumps = clumping_processes(proxy_df, sig_variants_df, sub_eGene_df, all_traits_are_one = tst)
    
    #Write clumped DFs
    for value in sig_sub_proxy_clumps.keys():
        if value in sig_sub_eGene_clumps.keys():
            sig_all_df = pd.DataFrame({value:[val for val in sig_all_clumps[value] if val]})
            #Turn from sets to strings
            sig_all_df['genes'] = sig_all_df[value].apply(lambda x: ','.join([y for y in x if not re.match('chr[1-9]*', y)]))
            sig_all_df[value] = sig_all_df[value].apply(lambda x: ','.join([y for y in x if re.match('chr[1-9]*', y)]))
            
            ## add source of genes to sig_all_df
            source_file={}
            
            df=pd.read_csv(sQTL_tbl,sep='\t',usecols=['sVariant','gene_id','gene_name','GWAS_variant']).drop_duplicates()
            source_file['s']=df[['GWAS_variant','gene_name']].groupby('gene_name').aggregate(set).squeeze()
            
            df=pd.read_csv(eQTL_tbl,sep='\t',usecols=['eVariant','gene_id','gene_name','GWAS_variant']).drop_duplicates()
            source_file['e']=df[['GWAS_variant','gene_name']].groupby('gene_name').aggregate(set).squeeze()
            
            if open_target is not None:
                source_file['o']=open_target_df.groupby('gene_name').aggregate(set).squeeze()

            sig_all_df['source']=""
            for i in range(sig_all_df.shape[0]):
                source_mark_all_genes=[]
                for gene in sig_all_df.iloc[i]['genes'].split(','):
                    source_mark_this_gene=''
                    for source in source_file.keys():
                        if source_file[source].get(gene, set([])).intersection(sig_all_df.iloc[i][value].split(',')):
                            source_mark_this_gene+=source
                    source_mark_all_genes.append(source_mark_this_gene)
                sig_all_df['source'].iloc[i]=','.join(source_mark_all_genes)
        
            sig_all_df.to_csv(out_dir+value+'_consolidated_clumping.tsv', sep='\t', index=False)
       
    print("Clumping null traits...")
    null_sub_proxy_clumps, null_sub_eGene_clumps, null_all_clumps = clumping_processes(proxy_df, null_variants_df, sub_eGene_df, all_traits_are_one = 'Null')


    value = next(iter(null_sub_proxy_clumps)) 
    null_all_df = pd.DataFrame({value:[val for val in null_all_clumps['Null'] if val]})
    #Turn from sets to strings
    null_all_df['genes'] = null_all_df['Null'].apply(lambda x: ','.join([y for y in x if not re.match('chr[1-9]*', y)]))
    null_all_df['Null'] = null_all_df['Null'].apply(lambda x: ','.join([y for y in x if re.match('chr[1-9]*', y)]))
    null_all_df.to_csv(out_dir+'Null'+'_consolidated_clumping.tsv', sep='\t', index=False)



def clumping_processes(proxy_df, variants_df, sub_eGene_df, all_traits_are_one = None):
    """
    Manager for the GWAS variant clumping.
    proxy_df: GWASVar2Gene table with proxy variant data
    variants_df: GWASVar2Gene traits table (or appropriate subset)
    eGene_df: GWASVar2Gene eGene table 
    all_traits_are_one: whether to treat all traits in variants_df as a single trait.
    If None, traits are treated individually. If string, that string is used to annotate
    all traits as one trait.
    """
    #Find proxy df subset that has one of our GWAS variants in LD
    proxy_df_ld = proxy_df.loc[proxy_df.LD_variant.isin(variants_df.GWAS_variant)].copy()
    sub_proxy_df = proxy_df_ld[['GWAS_variant', 'LD_variant']].drop_duplicates().copy()
    sub_proxy_df.drop_duplicates(inplace=True)
    if type(all_traits_are_one) is str:
        variants_df['real_trait'] = variants_df.trait
        variants_df['trait'] = all_traits_are_one
    ld_clumps = {}

    #clump by trait
    for trait in variants_df.trait.drop_duplicates().values:
        print(trait)
        sub_variants_df = variants_df.loc[variants_df.trait==trait].copy()
        #cluster by proxy variants
        sub_proxy_df_ld = sub_proxy_df.loc[sub_proxy_df.LD_variant.isin(sub_variants_df.GWAS_variant)].copy()
        #Change data types to fit consolidate function
        sub_proxy_df_ld.drop_duplicates(inplace=True)
        sub_proxy_df_ld.dropna(inplace=True)
        #verify df not empty
        if not sub_proxy_df_ld.empty:
            sub_proxy_df_ld['consolidate_this'] = sub_proxy_df_ld.apply(lambda x: [x.LD_variant,x.GWAS_variant], axis=1)
        else:
            continue
        ld_clumps[trait] = consolidate([set(x) for x in sub_proxy_df_ld[['consolidate_this']].values.flat])
  
    #Get GWAS_variant
    sub_eGene_df = sub_eGene_df[['GWAS_variant', 'gene_name']].copy()
    sub_eGene_df.drop_duplicates(inplace=True)
    traits_eGene_df = sub_eGene_df.merge(variants_df[['GWAS_variant', 'trait']].copy(), on='GWAS_variant').copy().drop_duplicates()
    #cluster by eGene
    egene_clumps = {}
    for trait in variants_df.trait.drop_duplicates().values:

        #get trait relevant lines
        sub_egene_df = traits_eGene_df.loc[traits_eGene_df.trait == trait][['GWAS_variant', 'gene_name']].drop_duplicates().dropna().copy()
        #change data to fit consolidate function
        #verify df not empty
        if not sub_egene_df.empty:
            sub_egene_df['consolidate_this'] = sub_egene_df.apply(lambda x: [x.GWAS_variant, x.gene_name], axis=1)
        else:
            continue
        egene_clumps[trait] = consolidate([set(x) for x in sub_egene_df[['consolidate_this']].values.flat])

    #Prepare for set consolidation between two groupings

    #merge LD and eGene clumps by trait
    final_clumps = {} 
    for trait in variants_df.trait.drop_duplicates().values:
        #Verify both values have key
        if ((trait in ld_clumps) & (trait in egene_clumps)):
            final_clumps[trait] = consolidate(ld_clumps[trait]+egene_clumps[trait])
        elif ((trait in ld_clumps) & (trait not in egene_clumps)):
            final_clumps[trait] = ld_clumps[trait]
        elif ((trait not in ld_clumps) & (trait in egene_clumps)):
            final_clumps[trait] = egene_clumps[trait]
        else:
            continue
    return(ld_clumps, egene_clumps, final_clumps)



def consolidate(sets):
    """
    function to consolidate sets into super-sets
    taken from here:
    https://rosettacode.org/wiki/Set_consolidation#Python:_Iterative
    """
    setlist = [s for s in sets if s]
    for i, s1 in enumerate(setlist):
        if s1:
            for s2 in setlist[i+1:]:
                intersection = s1.intersection(s2)
                if intersection:
                    s2.update(s1)
                    s1.clear()
                    s1 = s2
    return [s for s in setlist if s]
 
def flatten(container):
    """
    function to flatten list of lists into 1 list
    """
    for i in container:
        if isinstance(i, (list,tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i

def get_files_in_folder(folder):
    """
    folder: folder from which you want to retrieve files
    
    returns: list of files inside a folder
    """
    only_files = [os.path.join(folder, f) for f in os.listdir(folder)
                  if os.path.isfile(os.path.join(folder, f))]
    return(only_files)

def re_match(string, pattern, case_insensitive = True):
    """
    returns T/F indicating whether there is a match within string
    string: string to search for pattern
    pattern: pattern to look for
    case_insensitive: whether the function should disregard case
    """
    if case_insensitive:
        if (re.search(pattern, string, re.IGNORECASE)):
            return(True)
        else:
            return(False)
    else:
        if (re.search(pattern, string)):
            return(True)
        else:
            return(False)

if __name__ == "__main__":
    start_time = time.time()
    print("Starting...")
    main()
    print("--- Done in %s seconds ---" % (time.time() - start_time))




