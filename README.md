## ECLIPSER
This repository contains the software implementation of ECLIPSER (Enrichment of Causal Loci and Identification of Pathogenic cells in Single Cell Expression and Regulation data).

ECLIPSER tests whether the expression of genes mapped to GWAS loci of complex diseases or traits, based on eQTLs/sQTLs and other functional data, are enriched in specific cell types in one or more tissues.

Authors: Jiali Wang, John Rouhana written in the Segre lab under the guidance of Ayellet Segre, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA

Contributors to code: Andrew Hamel, Brian Cole, Massachusetts Eye and Ear, Harvard Medical School

Date: August 2021

## Repository structure
`src`: the directory contains scripts for the software pipeline and for plotting results with relevant instructions

`data`: the directory contains input files required to run ECLIPSER

`examples`: sample data for running ECLIPSER for broad cell type annotations in GTEx tissues and 21 complex traits. 

## Dependencies
ECLIPSER was written in Python (> 3.6) and requires the following libraries and modules:

- pandas 1.1.3
- numpy 1.18.1
- matplotlib 3.3.1
- seanborn 0.11.0 
- statsmodels 0.11.1
- plotnine 0.8.0

## Publications
Gokcen Eraslan*, Eugene Drokhlyansky*, Shankara Anand‡, Evgenij Fiskin‡, Ayshwarya Subramanian‡, Michal Slyper‡, Jiali Wang‡, Nicholas Van Wittenberghe, John M. Rouhana, Julia Waldman, Orr Ashenberg, Monkol Lek, Danielle Dionne, Thet Su Win, Michael S. Cuoco, Olena Kuksenko, Alexander M. Tsankov, Philip A. Branton, Jamie L. Marshall, Anna Greka, Gad Getz,  Ayellet V. Segrè#,  François Aguet#,  Orit Rozenblatt-Rosen#,  Kristin G. Ardlie#,  Aviv Regev#. Single-nucleus cross-tissue molecular reference maps towards understanding disease gene function. Science. 2022. 376(6954): eabl4290. DOI: 10.1126/science.abl4290. Co-first authors, these authors contributed equally. ‡Co-second authors, these authors contributed equally. #Co-senior and co-corresponding author.

John M. Rouhana*, Jiali Wang*, Gokcen Eraslan, Shankara Anand, Andrew R. Hamel, Brian Cole, Aviv Regev, François Aguet, Kristin G. Ardlie, Ayellet V. Segrè. ECLIPSER: identifying causal cell types and genes for complex traits through single cell enrichment of e/sQTL-mapped genes in GWAS loci. BioRxiv, 2021. *Co-first authors. doi: https://doi.org/10.1101/2021.11.24.469720
