#!/bin/bash

## Example bash script to perform enrichment analysis for each cell type, trait and tissue combination for GTEx snRNA-seq tissues (Eraslan et al, BioRxiv 2021).  

# See README.md for required python libraries
# pip install statsmodels

## Run jobs sequentially
for clumping_file in /example/GTEx_clumping_files/*/*consolidated_clumping.tsv 
do
for tissue_file in /example/GTEx_clumping_files/*
do  

    mkdir ../out_ECLIPSER_GTEx/

    folder=$(echo $tissue_file |  rev | cut -f1- -d '/' | rev)  
    trait=$(echo $clumping_file | sed 's/ /_/g' | rev | cut -f1 -d'/' | rev | sed 's/_consolidated_clumping.tsv//g')
    tissue=$(echo $folder | rev | cut -d '/' -f1 | rev)

    test=$(echo $clumping_file | rev | cut -f1 -d'/' | rev)
    if [[ "$test" != "Null_consolidated_clumping.tsv" ]]
    then
    python -u /src/enrichment.py \
      --sc_path /data/GTEx_snRNAseq_DGE_broad_cell_types.csv \
      --clumps_file $clumping_file  \
      --background_file $folder/Null_consolidated_clumping.tsv \
      --out_directory ../out_ECLIPSER_GTEx/ \
      --trait $trait \
      --tissue $tissue \
      --method Bayesian
    fi 
done 
done


## Run jobs in parallel using LSF job scheduler
for clumping_file in /example/GTEx_clumping_files/*/*consolidated_clumping.tsv 
do
for tissue_file in /example/GTEx_clumping_files/*
do  

    mkdir ../out_ECLIPSER_GTEx/

    folder=$(echo $tissue_file |  rev | cut -f1- -d '/' | rev)  
    trait=$(echo $clumping_file | sed 's/ /_/g' | rev | cut -f1 -d'/' | rev | sed 's/_consolidated_clumping.tsv//g')
    tissue=$(echo $folder | rev | cut -d '/' -f1 | rev)

    test=$(echo $clumping_file | rev | cut -f1 -d'/' | rev)
    if [[ "$test" != "Null_consolidated_clumping.tsv" ]]
    then
    bsub -J test -q short -n 1 -M 1000 -R rusage[mem=1000] -e err -o out "python -u /src/enrichment.py \
      --sc_path /data/GTEx_snRNAseq_DGE_broad_cell_types.csv \
      --clumps_file $clumping_file  \
      --background_file $folder/Null_consolidated_clumping.tsv \
      --out_directory ../out_ECLIPSER_GTEx/ \
      --trait $trait \
      --tissue $tissue \
      --method Bayesian"
    fi 
done 
done
