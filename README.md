# vcf2MK

Overview: vcf2MK is a program for annotated polymorphism and divergence data within coding sequences in a vcf formatted DNA polymorphism file. This program is distributed with a number of scripts for assembling MK tables among other things.  

Downloading and Compilation: The vcf to MK package can be acquired via github at https://github.com/russcd/vcf2MK, using the command: “git clone https://github.com/russcd/vcf2MK.git” . The program can then be compiled via the command: “c++ -o vcf2MK –-std=c++0x vcf2MK.cpp” 

Usage: To use vcf2MK, the user must supply the program with three files. First, the user must supply a vcf file of polymorphism data from within a focal population (using the –v option). The user must also supply an aligned outgroup file, either in vcf of maf formats (using the –o option). This file must have the correct file extension (either .vcf or .maf). All vcf files supplied to vcf2MK must contain only the first eight columns of the vcf file. This can be accomplished using the command “cut –f1-8” as shown below. Both the ingroup data, and outgroup must be aligned to the same reference genome. Finally, the user must supply an annotation file, in gff, gff3 or gtf formats. Again, the annotation file extension must be correct.

An example usage is this:

mkfifo in.vcf out.vcf
cut –f1-8 ingroup.vcf > in.vcf &
cut –f1-8 outgroup.vcf > out.vcf &
./vcf2MK –v in.vcf –o out.vcf –c annotation.gff > transcript_data.txt

Little quality filtration is built into the program. Instead, the user is expected to provide vcf2MK with a heavily filtered vcf file for both ingroup and outgroup vcfs.  

Output: The program output is tab-delimited data for all variable sites within a transcript supplied in the annotation file. The format is chromosome, position, number of nonsynonymous divergences, number of synonymous divergences, number of nonsynonymous polymorphisms, number of synonymous polymorphisms, number of genotyped chromosomes in vcf ingroup, frequency of polymorphic site (if a polymorphic site is present), and a comma separated list of transcripts that share this position. The following are two example output lines. 
 
2L	70661	0	1	0	1	192	0.932292	FBtr0306536,FBtr0306538,FBtr0306539,FBtr0078100,FBtr0306537,

2L	70669	1	0	0	0	192	NA	FBtr0306536,FBtr0306538,FBtr0306539,FBtr0078100,FBtr0306537,

Columns 3 and 4 may contain any value between 0 and 1. The reason for this is that when more than a single divergent site is present within a codon, there may be multiple mutational trajectories between the two codons. vcf2MK will report the average nonsynonymous and synonymous counts of all possible mutational paths between the two codons. Because, in most cases, phase is not known for polymorphic sites, in instances where there are two or more polymorphisms within a single codon, each site is assumed to arise independently on a codon consisting of the highest frequency bases at the other two sites within the codon sequence. 

Eample Usage: Example files are provided in the example directory. To test the program, use the command “./vcf2MK –v Dmelanogaster.vcf –o Dsimulans.vcf –c Dmelanogaster.gff > Dmelanogaster.detailed”

Scripts: Two scripts are contained with the vcf2MK package. detailed2site_averages.pl will take the output from vcf2MK and aggregate information for all tanscripts that contain a given site, and output the average dn, ds, pn, and ps counts for that site across all transcripts. The package also contains a script, detailed2mk.pl, that will aggregate site counts across each transcript and output totals. By default this script will discard all sites with minor allele frequencies less than 0.1. The output of this script is tab-delimited with the transcript id, dn, ds, pn, and ps. 

