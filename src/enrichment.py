# Author: Jiali Wang
# Segre lab, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA
# Date: April 2021

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys 
import os
from glob import glob
from statsmodels.stats.multitest import fdrcorrection
import pandas as pd
import numpy as np
import time
import argparse

def main():
    #Define user-defined parameters
    parser = argparse.ArgumentParser(description='GWAS clumping for single cell enrichment analysis')
    
    parser.add_argument('--sc_path', type=str, required=True,
                    help='Path to the comma delimited single cell differential expression table.')
    parser.add_argument('--clumps_file', type=str, required=True,
                    help='Path to the clump file of trait of interest from GWAS clumping output.')
    parser.add_argument('--background_file', type=str, required=True,
                    help='Path to the clump file of background from GWAS clumping output.')
    parser.add_argument('--trait', type=str, required=True,
                    help='Name of the trait of interest.')
    parser.add_argument('--tissue', type=str, required=True,
                    help='Name of the tissue. Should match the values in the tissue column.')
    parser.add_argument('--out_directory', type=str, required=True,
                    help='Directory where output files will be written. Default to working directory.')
    parser.add_argument('--num_samples', type=int, default=1000000,
                    help='Number of MC sampling to take each run.')
    parser.add_argument('--cutoff1', type=float, default=0.5,
                    help='cutoff value of log2FC')
    parser.add_argument('--cutoff2', type=float, default=0.1,
                    help='cutoff value of pval FDR (pvals_fdr)')
    parser.add_argument('--cutoff3', type=float, default=0.05,
                    help='cutoff value of % cells expressing gene of given cell type (pct.1)')
    parser.add_argument('--method', type=str, required=True,default='Bayesian',
                    help='enrichment method, Bayesian or Permutation')
        
    args = parser.parse_args()

    def f_genes(x):
        """
        Compute number of genes in sc file and fraction log2FC>cutoff1 & pval_fdr<cutoff2 in GWAS locus and cell type x
        """
        genes_df=sc_df.reindex(x.split(',')).dropna()
        out=genes_df.mean()
        out['num_genes_in_sc']=genes_df.shape[0]
        return out

    
    # read data
    sc_df0 = pd.read_csv(args.sc_path, sep=',')
    # prepare sc data for heart tissue, rows are genes, columns are cell types, values are indicators of log2FC>cutoff1 & pval_fdr<cutoff2
    sc_df0=sc_df0.loc[sc_df0['tissue']==args.tissue]

    if "pct.1" in sc_df0.columns:
        sc_df=(sc_df0.pivot(index='gene', columns='celltype', values='log2FC')>args.cutoff1) & (sc_df0.pivot(index='gene', columns='celltype', values='pvals_fdr')<args.cutoff2) & (sc_df0.pivot(index='gene', columns='celltype', values='pct.1')>=args.cutoff3)
    else:
        sc_df=(sc_df0.pivot(index='gene', columns='celltype', values='log2FC')>args.cutoff1) & (sc_df0.pivot(index='gene', columns='celltype', values='pvals_fdr')<args.cutoff2)
    
    clumps_df=pd.read_csv(args.clumps_file, sep='\t',index_col=0).dropna()
    numps_df=pd.read_csv(args.background_file, sep='\t',index_col=0).dropna()
    
    ## Compute number of genes in sc file and fraction of log2FC>cutoff1 & pval_fdr<cutoff2 in GWAS locus and cell type in clumps file and numps file
    clumps_df=clumps_df.join(clumps_df['genes'].apply(f_genes)).dropna()
    numps_df=numps_df.join(numps_df['genes'].apply(f_genes)).dropna()
    
    null_95_percentile=numps_df.quantile(0.95) 
    avg_thr_stat=clumps_df.drop(['genes', 'num_genes_in_sc','source'], axis=1).mean() #observed statistics from GWAS catalog
    
    
    ## table 2
    if args.method == 'Bayesian':
        res=[]
        for cell in clumps_df.columns.drop(['genes', 'num_genes_in_sc','source']):  
            a=np.sum(numps_df.loc[:,cell]>=null_95_percentile[cell])
            b=numps_df.shape[0]-a
            c=np.sum(clumps_df.loc[:,cell]>=null_95_percentile[cell])
            d=clumps_df.shape[0]-c
            theta1=np.random.beta(c+1, d+1, size=args.num_samples)
            theta2=np.random.beta(a+1, b+1, size=args.num_samples)

            fold_enrichment=theta1/theta2
            res_fold_enrichment=pd.Series(fold_enrichment).describe([0.025, 0.5, 0.975]).T[['2.5%','50%','97.5%']]
            res_fold_enrichment.index='fold_enrichment_'+res_fold_enrichment.index
            res_fold_enrichment['Enrichment p-value']=(fold_enrichment<1).mean()
            res_fold_enrichment.name=cell

            res.append(res_fold_enrichment)
        table2=pd.DataFrame(res)  
        
    elif args.method == 'Permutation':
        null_95_percentile1=np.expand_dims(np.expand_dims(null_95_percentile.drop('num_genes_in_sc').values, axis=0), axis=0)
        numps_df1=numps_df.drop(['genes', 'num_genes_in_sc'], axis=1).values  #drop 'genes' and 'num_genes_in_sc' columns
        clumps_df1=clumps_df.drop(['genes', 'num_genes_in_sc','source'], axis=1).values
        res_batch=[]
        for batchs in range(100):
            boot_idx=np.random.randint(low=0,high=numps_df1.shape[0],size=int(clumps_df1.shape[0]*args.num_samples/100))
            temp=numps_df1[boot_idx,:]
            res_batch.append(np.sum(np.reshape(temp,[clumps_df1.shape[0],int(args.num_samples/100),numps_df1.shape[1]])>=null_95_percentile1,axis=0))

        null_thre_df=pd.DataFrame(np.concatenate(res_batch,axis=0),columns=numps_df.columns.drop(['genes', 'num_genes_in_sc'])) #concatenate results from 100 batches

        c=(clumps_df.drop(['genes', 'num_genes_in_sc','source'], axis=1)>=null_95_percentile.drop('num_genes_in_sc')).sum()
        table2=((c+1)/(null_thre_df+1)).describe([0.025,0.5,0.975]).T[['2.5%','50%','97.5%']]
        table2.columns='fold_enrichment_'+table2.columns
        table2['Enrichment p-value']=np.mean(null_thre_df>=c)

        
        
    table2['BH adj. Enrichment p-value (tissue-wide)']=fdrcorrection(table2['Enrichment p-value'], alpha=0.05)[1]
    table2['Average GWAS locus set statistic']=avg_thr_stat
    table2['Fraction of GWAS loci >95% percentile of null loci']=np.sum(clumps_df>=null_95_percentile)/clumps_df.shape[0]
    
    # add leading edge genes: genes in GWAS loci with score above 95th percentile of null loci
    for cell in numps_df.columns.drop(['genes', 'num_genes_in_sc']):
        qt95=numps_df[cell].quantile(0.95)
        GWAS_locus_95=clumps_df.loc[clumps_df[cell]>=qt95,['genes','source']].reset_index()

        GWAS_loci_df=[]
        gene_log2FC=[]
        for i in range(GWAS_locus_95.shape[0]):
            gene_list       =GWAS_locus_95.iloc[i,:]['genes'].split(',')
            gene_list_source=GWAS_locus_95.iloc[i,:]['source'].split(',')
            
            # gene_list get from sc_df, nan is not a good gene, replace with False
            gene_list=sc_df.reindex(gene_list)[cell].replace(float('nan'),False)
            good_genes       =gene_list[gene_list].index.tolist() 
            good_genes_source=np.array(gene_list_source)[gene_list].tolist()
            gene_log2FC.append(np.mean(sc_df0.loc[(sc_df0['gene'].isin(good_genes))&(sc_df0['celltype']==cell)&(sc_df0['tissue']==args.tissue),'log2FC']))

            if len(good_genes)!=0:
                GWAS_loci_df.append([GWAS_locus_95.iloc[i,0],','.join(good_genes),','.join(good_genes_source)])

        gene_log2FC=np.array(gene_log2FC)
        gene_log2FC=gene_log2FC[np.logical_not(np.isnan(gene_log2FC))].tolist()        
        GWAS_loci_df=pd.DataFrame(GWAS_loci_df,columns=['loci','genes','source'])
        GWAS_loci_df['ave_gene_log2FC']=gene_log2FC
        GWAS_loci_df=GWAS_loci_df.merge(clumps_df[cell],how='left',left_on='loci',right_index=True)
        GWAS_loci_df=GWAS_loci_df.sort_values(by=[cell,'ave_gene_log2FC'],ascending=False)
        table2.loc[cell,'Leading edge genes']=";".join(GWAS_loci_df['genes'].tolist())
        table2.loc[cell,'Leading edge GWAS loci']=";".join(GWAS_loci_df['loci'].tolist())
        table2.loc[cell,'Leading edge genes source']=";".join(GWAS_loci_df['source'].tolist())
       
    table2['GWAS']=args.trait #remove '_consolidated_clumping.tsv' suffix
    table2.index.name='Cell type'
    table2['Tissue']=args.tissue
    
    table2.rename({'fold_enrichment_2.5%': 'Fold-enrichment lower 95% CI',
          'fold_enrichment_50%':'Fold-enrichment','fold_enrichment_97.5%':'Fold-enrichment upper 95% CI'}, axis=1, inplace=True)
    
    table2.to_csv(args.out_directory+'table_'+args.trait+'_'+args.tissue+'.csv')
    
if __name__ == "__main__":
    start_time = time.time()
    print("Starting...")
    main()
    print("--- Done in %s seconds ---" % (time.time() - start_time))


