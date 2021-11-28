#!/bin/bash

# Example: Run clumping on GWAS variants for eight skin diseases and traits taken from Open Targets Genetics

# See README.md for required python libraries (e.g. numpy/1.7.1-goolf-1.4.10-Python-2.7.3, pandas/2.5.0)
# The GWASvar2gene input files used for the clumping of the skin trait GWAS variants and gene to variant mapping can be downloaded from Zenodo and saved in the data/ folder.

# Run jobs in parallel on 8 skin traits using SLURM job scheduler

sbatch --mem=0 --nodes=2 --ntasks-per-node=1 --job-name ECLIPSER_skin_traits -e err_file -o out_file --wrap="python -u /src/clumping.py \
--proxy_table /data/GWASvar2gene_open_targets_proxy_table.tsv \
--traits_table /data/GWASvar2gene_open_targets_trait_table.tsv \
--eQTL_table /data/GWASvar2gene_open_targets_eQTL_table.tsv \
--sQTL_table /data/GWASvar2gene_open_targets_sQTL_table.tsv \
--opentarget /data/opentargets_variant_and_gene_table.csv \
--significant_traits /data/Reynolds_etal_skin_traits/skin_traits.txt \
-ancestry EUR \
-out_directory /example/Skin_clumping_files/"


# Run job on unix command line:
 
 python -u /src/clumping.py \
   --proxy_table /data/GWASVar2gene_open_targets_proxy_table.tsv \
   --traits_table /data/GWASvar2gene_open_targets_trait_table.tsv \
   --eQTL_table /data/GWASvar2gene_open_targets_eQTL_table.tsv \
   --sQTL_table /data/GWASvar2gene_open_targets_sQTL_table.tsv \
   --opentarget /data/opentargets_variant_and_gene_table.csv \
   --significant_traits /data/Reynolds_etal_skin_traits/skin_traits.txt \
   -ancestry EUR \
   -out_directory \"$out_path\" 
done
