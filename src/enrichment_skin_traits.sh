#!/bin/bash

# See README.md for required python libraries 
# pip install statsmodels

## Example bash script to perform enrichment analysis for each cell type, trait and tissue pair for eight skin diseases and traits using 
## differential gene expression from Reynold et al., Science 2021

## Run jobs sequentially
# Loop over all tissues (just skin in this example)

mkdir ../out_ECLIPSER_skin

for clumping_file in /example/Skin_clumping_files/*/*consolidated_clumping.tsv 
do
for tissue_file in /example/Skin_clumping_files/*
do  

    folder=$(echo $tissue_file |  rev | cut -f1- -d '/' | rev)  
    trait=$(echo $clumping_file | sed 's/ /_/g' | rev | cut -f1 -d'/' | rev | sed 's/_consolidated_clumping.tsv//g')
    tissue=$(echo $folder | rev | cut -d '/' -f1 | rev)

echo $folder
echo $trait
echo $tissue

    test=$(echo $clumping_file | rev | cut -f1 -d'/' | rev)
    if [[ "$test" != "Null_consolidated_clumping.tsv" ]]
    then
    python -u  /src/enrichment.py \
      --sc_path /data/skin_normal_final_clustering_DGE_Reynolds_etal_Science2021.csv \
      --clumps_file $clumping_file \
      --background_file $folder/Null_consolidated_clumping.tsv \
      --out_directory ../out_ECLIPSER_skin/ \
      --trait $trait \
      --tissue $tissue \
      --method Bayesian \
      --cutoff1 0.5
      --cutoff2 0.1
      --cutoff3 0.05
    fi 
done 
done
