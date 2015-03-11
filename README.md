# Massive genomic variation and strong selection in Arabidopsis thaliana lines from Sweden

> Long Q, Rabanal FA, Meng D, Huber CD, Farlow A, Platzer A, Zhang Q, Vilhjálmsson BJ, Korte A, Nizhynska V, Voronin V, Korte P, Sedman L, Mandáková T, Lysak MA, Seren Ü, Hellmann I, Nordborg M (2013) Massive genomic variation and strong selection in Arabidopsis thaliana lines from Sweden. Nat Genet 45: 884–890.

 
## Genotype files

In total we identiﬁed around 4.5M SNPs, 576k short indels, 23k transposable elements, 7.7k CNVs, as well as 3.8k other structural variants (larger than 200 bp). The following ﬁles are available:

* [Raw data files](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=search_obj&m=&s=&term=SRP012869&go=Search) (externally hosted at NCBI trace SRA)
* SNPs
 * [Original](SNPs/original) (+ annotations)
 * [Imputed](SNPs/imputed) (+ annotations)
 * [ReadMe](SNPs/)
* Small indels
 * [Original](Indels/original) (+ annotations)
 * [Imputed](Indels/imputed)
 * [ReadMe](Indels/)
* [Large structural variants](SV/)


## Programs

### LAE-finder

* [Program](programs/LAE-finder/bin)
* [Source](programs/LAE-finder/src)
* [README](programs/LAE-finder/README.MD)

### Denovo-SV

* [Program and readme](programs/denovosv)

### MACH related scripts

These are the scripts that format genotypes files and submit MACH jobs to our cluster. Please note that these scripts are not intended to be easily usable elsewhere.

* [Scripts](programs/mach_scripts)


## Other files

* [Lines used in this study](https://www.google.com/fusiontables/DataSource?docid=1G500e9BsyalkbxWyyzrXlUa7odlwRGqUkf_k7h0#rows:id=1) (Google Docs)
* [Alignment and assembly statistics for all sequenced lines](files/sequencing_statistics.csv)
* Flow cytometry and repeat copy number estimates for the 180 lines and 36 worldwide lines
 * [Flow cyto estimates 180 lines](files/genome_size180.csv)
 * [Flow cyto estimates 36 worldwide lines](files/genome_size36.csv)
* [Multiple alignments for candidate genes from Fig 2](files/candidate_alignments.zip)
* [Summaries of putative sweep regions](files/sweeps.xlsx)
* [PCR primers used in this study](files/primers.xlsx)
* [Predicted genotype for 1306 lines with respect to the transposition on chromosome 1, and the large inversion on chromosome 4 (Supplementary Fig. 17).](files/sweep_knob.csv)
* [Genes in the swept transposition on chromosome 1](files/genes_in_transposition.xlsx)
* [North and south labels](files/northAndSouth.csv)

