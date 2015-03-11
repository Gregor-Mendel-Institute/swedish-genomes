there are three types of files in the subfolders:

info file columns:
1. chromosome;
2. position;
3. average Q-value
4. reference allele;
5. non-reference counts: the number of individuals that are called as alternative allele;
6. alternative allele;
7. long-range LD flag: The flag denoting whether there is a SNP that are in long range LD with this SNP, indicating potential mapping artifacts.


csv-files (chr1-5)
The first line is the header for all the sample IDs. In the remaining lines, the first two columns are chromosome and coordinate, followed by the actual SNP calls. “N” means that the call for this site is missing.


annotations file columns:
1. chromosome
2. coordinate
3. reference allele
4. ancestral allele (“N” denotes no alignment, “-” denotes deleted, “n” means that there is alignment but with low identity sequences surrounding the site)
5. alternative allele
6. non-reference counts: the number of individuals that are called as alternative allele
7. synonymous or non-synonymous changes ('*' denotes non-coding sequences.)
8. gene(s) located. If two genes separated by ";" it means this SNP is located in two overlapping genes
9. gene function: protein coding genes or RNA genes or TE genes, etc.
10. in exon or intron. If in multiple genes, they are separated by ";". If there are multiple isoforms of the same genes, the details will be listed with “ix :” denoting the x-th isoform. Different isoforms are separated by ",".
11. if exon, in coding sequence or UTR
12. the codon, if in coding sequence. The capitalized letter denotes the focal SNP. Codons that are in negative strand genes are reverse-complemented.


only in the original folder is an info file, because the imputed SNPs are a subset of the original (and the SNP info is in both cases the same).

