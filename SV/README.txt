formats in this folder:

SV-info columns (from the title line)
chr . . . . . . chromosome (1-5)
loc . . . . . . locus
len . . . . . . length of the event
event_type_ref  either a shortcut for the event or a sequence in case of insertion and if the calling method reported it
non_ref_counts  how many lines have the non-reference alternative called
anc_status .  . ancestral allele ("N" denotes no alignment, "0" denotes ancestral reference allele, "1" denotes ancestral non-reference allele, "c" or "C" denotes cases where there is alignment, but the ancestral allele is neither the reference nor the alternative allele.)
read_pair_support   how many read pairs support this event. (This format was designed to be shared by both small and large structural variants. But for the small indels called from local realignment, there is no read-pair support. For some pipelines that do not leverage read pair (e.g., coverage counting), this info is not available.)
bp_range1 and 2 . . if the method reports it, the coordinates of the 2 breakpoints; should be quite similar to locus and len.
four_gamete_left  ...
four_gamete_right  ...
call_method . . . . methods supporting this event


SV-csv/table
The format is the same for the SNP file, with "-1", "1" and "0" denoting the no-call, alternative allele and reference allele, respectively.


TEs_info columns (from the title line)
chr+loc . . . . locus
len . . . . . . The length of the corresponding reference event.
event_type_ref  The class of this event annotated (resp. the item/TE)
non_ref_counts  The amount of individuals sharing this event.
anc_status . .  unused
read_pair_support  The total number of all supporting read pairs of all individuals.
bp_range1 + 2 .. unused
four_gamete_left + _right  unused
call_method  For TE-Locate here is written 'PairEndTE', used if merged with other data in this format.
Orientation  'parallel', 'inverse' or 'uncertain': The orientation according to the reference sequence.
#pPairs . .  The number of read pairs supporting parallel orientation. Not used if the orientation is 'uncertain'.
#iPairs . .  The number of read pairs supporting inverse orientation. Not used if the orientation is 'uncertain'.
new/old . . 'new' or 'old'. 'old' if the item is called at the locus in the reference, 'new' otherwise. Note that at higher hierarchical levels all locations of this item are meant, e.g. any Copia called at a Copia locus in the reference is called 'old' as the item's name is the only distinction.


TEs_table
the same as the other csvs, except that there are the number of supporting reads in each cell.


comments:
- As for the indel data the coordinates for indels are according to the vcf-format (see e.g. http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40 ). In brief that means the coordinate points at the base pair before the event: an insertion takes place between the bp of the coordinate and the next one, a deletion takes place from the next bp (so the bp of the locus is the last one still present).
- the 'event_type_ref' column in the SV-info file contains the info from the best source for this event (sources are the different tools, these are as a list in the 'call_method' column). The -<sequence> in the 'event_type_ref' column comes only from Pindel (http://bioinformatics.oxfordjournals.org/content/25/21/2865.short).

