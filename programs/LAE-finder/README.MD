LAE-finder 1.0

is brought to you by:
        Alexander Platzer ( alexander.platzer@gmi.oeaw.ac.at )

LAE-finder is a tool to confirm or discard known large events for certain
individuals.

###############################################################################

1. Prerequisites
   1.1 Copy the package
   1.2 SAM files
   1.3 Known large events

2. Processing
   2.1 Checking of events
   2.2 Filter read pairs
   2.3 Find support or restriction of events
   2.4 Collect the result files
   
3.  Output format

4.  Usage of the LAE-finder executable


###############################################################################

1. Prerequisites

1.1 Copy the package

    Copy groupINVs.pl, groupTLs.pl and LAE-finder to one folder. The location
    is later referred as main folder.


1.2 SAM files

    Generate SAM files of the read-pairs of your accessions and move them to
    one folder, e.g. data/ .
    ! The SAM files must be sorted lexically !
    In Linux you can use the 'sort' command for this, e.g.:    
    
    sort --temporary-directory=. <sam file> > <sorted sam file>

    See for SAM file format : http://samtools.sourceforge.net
    Aligners producing this format are e.g. : 
           bwa ( http://bio-bwa.sourceforge.net/ )
           SMALT ( http://www.sanger.ac.uk/resources/software/smalt/ )
           segemehl ( www.bioinf.uni-leipzig.de/Software/segemehl/ )
           ...


1.3 Known large events

    A list of known large events in the format:
    
    Type	chr1	BP1.1	BP1.2	BP2.1	BP2.2	Chr2	BP3.1	BP3.2	
    INV	1	8590024	8593298	9475376	9481018
    INV	1	9294309	9295509	9317112	9317112
    TL	2	8343119	8333012	8753433	8755591	2	8121925	8123161
    TL	2	6947830	6985392	7155186	7155377	2	8753433	8755591	
    ...
    
    With '.1' and '.2' the start and the end of a range of a possible 
    breakpoint is meant. For translocations BP3 is the breakpoint where the
    sequence is inserted.
    
    
2. Processing

2.1 Checking of events

    It is a good idea to check the events first e.g. with
    
    ./LAE-finder nothing ../SVs_02.dat ../data/adal_1.9321.realigned.sort.bam.sam.lsort
    
    the second parameter means nothing is done, but the events are checked in either case.
    
    
2.2 Filter read pairs

    The later analysis/counting uses just the read pairs within the events. So
    the first steht is to filter for these, e.g. :
    
    ./LAE-finder filter ../SVs_02.dat ../data/adal_1.9321.realigned.sort.bam.sam.lsort

    The resulting sam file gets a '.sel' to its name, so in this case -> adal_1.9321.realigned.sort.bam.sam.lsort.sel


2.3 Find support or restriction of events

    In the beforehand filtered reads is searched for events:

    ./LAE-finder SVcall SVs_02.dat data/adal_1.9321.realigned.sort.bam.sam.lsort.sel
    
    if no error occures, there are two result files, one named with appended 
    '.INVs.dat' one with appended '.TLs.dat', so in this case:
    
    adal_1.9321.realigned.sort.bam.sam.lsort.INVs.dat
    adal_1.9321.realigned.sort.bam.sam.lsort.TLs.dat
    

2.4 Collect the result files

    If lots of accessions are analyzed, for each accession there are two result
    files. To have an better overview per event the two scripts might be helpful:

    perl groupINVs.pl
    perl groupTLs.pl

    There will be then one file per event named: INVxxx.dat resp. TLxxx.dat,
    where xxx stands for the numeration.


3.  Output format

    Beside the first column, which is either the (arbitrary, but consistent
    between individuals) numbered event (e.g. INV3, TL45) or the id of the
    individual (which is just the string before the first dot), the output
    format is,
    for inversions:

    ID : The accession
    <INV coordinates>
    #reads supporting
    #continuousCovsingle
    #continuousCovpair
    meanCov
    (the last three are respectively for side A and B)
    
    for translocations:
    ID : The accession
    <TL coordinates>
    #reads supp from / to / deletion
    #reads supp. forward / backward
    #stacks in BP
    #overlapping supp reads
    best stack direction ratio from / to
    #continuous cov single
    #continuous cov pair
    mean cov
    (the last 3 are for all three sides)
    bestBPC (start and end of BP1)


4.  Usage of the LAE-finder executable

    ./LAE-finder <subcommand> <SV-list file> <sam file>

    in detail:
    <subcommand>     can be 'filter', 'SVcall' or 'filterAndSVcall'; anything
                     else is also valid, in this case only the starting checks
                     are done
    <SV-list file>	a file with the known large events, format descriped in
                     section 1.3
    <sam file>       the sam file for this individual. See section 1.2

    the <subcommand> in detail:
    
    'filter'           section 2.2
    'SVcall'           section 2.3
    'filterAndSVcall'  'filter' and 'SVcall' at once
    <anything else>    just the feedback of the starting checks (or nothing if
                       everything looks correct at this stage), the parameter 
                       can be e.g. 'nothing' as in section 2.1


