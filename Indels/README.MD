# Indels

There are three types of files in the subfolders, and these are based on the SNP-formats, extended for indels:


info file columns:
1. chromosome
2. coordinate
3. length of the event
4. alternative allele (with respect to the reference genome)
5. non-reference counts: the number of individuals that are called as alternative allele
6. ancestral allele (“N” denotes no alignment, “0” denotes ancestral reference allele, “1” denotes ancestral non-reference allele, “c” or “C” denotes cases where there is alignment, but the ancestral allele is neither the reference nor the alternative allele.)
7. read-pair-support: how many read pairs support this event. (This format was designed to be shared by both small and large structural variants. But for the small indels called from local realignment, there is no read-pair support. For some pipelines that do not leverage read pair (e.g., coverage counting), this info is not available.)
8. bp-range: the range of breakpoint when repetitive sequences make it difficult to find the exact position.
9. call-method: the method(s) that called this event.


## csv-files (chr1-5)
The format is the same for the SNP file, with "-1", "1" and "0" denoting the no-call, alternative allele and reference allele, respectively.


## annotations file columns:
The same format as for SNPs, but with two lines per polymorphism for the starting and ending point of that indel, respectively.


## Comments:
Be aware that the coordinates for indels are according to the vcf-format (see e.g. http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40 ). In brief that means the coordinate points at the base pair before the event: an insertion takes place between the bp of the coordinate and the next one, a deletion takes place from the next bp (so the bp of the locus is the last one still present).

