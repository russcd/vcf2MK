# vcf2MK

Overview: vcf2MK is a program for annotated polymorphism and divergence data within coding sequences in a vcf formatted DNA polymorphism file. This program is distributed with a number of scripts for assembling MK tables among other things.  

Downloading and Compilation: The vcf to MK package can be acquired via github at https://github.com/russcd/vcf2MK, using the command: “git clone https://github.com/russcd/vcf2MK.git” . The program can then be compiled via the command: “g++ -o vcf2MK -std=c++0x vcf2MK.cpp” 

Usage: To use vcf2MK, the user must supply the program with three files. First, the user must supply a vcf file of polymorphism data from within a focal population (using the –v option). The user must also supply an aligned outgroup file, either in vcf of maf formats (using the –o option). This file must have the correct file extension (either .vcf or .maf). All vcf files supplied to vcf2MK must contain only the first eight columns of the vcf file. This can be accomplished using the command “cut –f1-8” as shown below. Both the ingroup data, and outgroup must be aligned to the same reference genome. Finally, the user must supply an annotation file, in gff, gff3 or gtf formats. Again, the annotation file extension must be correct.

An example usage is this:

mkfifo in.vcf out.vcf

cut –f1-8 ingroup.vcf > in.vcf &

cut –f1-8 outgroup.vcf > out.vcf &

./vcf2MK –v in.vcf –o out.vcf –c annotation.gff > transcript_data.txt

Little quality filtration is built into the program. Instead, the user is expected to provide vcf2MK with a heavily filtered vcf file for both ingroup and outgroup vcfs.  

Output: The program output is tab-delimited data for all variable sites within a transcript supplied in the annotation file. The format is 

1. chromosome
2. position
3. outgroup allele
4. ancestral (or  highest frequency) allele
5. polymorphic allele
6. number of nonsynonymous divergences
7. number of synonymous divergences
8. number of nonsynonymous polymorphisms
9. number of synonymous polymorphisms
10. number of genotyped chromosomes in vcf ingroup
11. frequency of polymorphic site (if a polymorphic site is present)
12. comma separated list of transcripts that share this position. 
13. The fold-degeneracy of that position in the codon containing the ingroup ancestral allele.

The following are two example output lines. 
 
2L	70293	A	G	T	1	0	1	0	192	0.00520833	FBtr0306539,FBtr0078100,FBtr0306536,FBtr0306537,FBtr0306538,	0

2L	70316	C	C	A	0	0	0	1	192	0.00520833	FBtr0306539,FBtr0078100,FBtr0306536,FBtr0306537,FBtr0306538,	4

Eample Usage: Example files are provided in the example directory. To test the program, use the command “./vcf2MK –v Dmelanogaster.vcf –o Dsimulans.vcf –c Dmelanogaster.gff > Dmelanogaster.detailed”

Options: "-i" will output the full set of sites with coverage in both ingroup vcf and outgroup alignment regardless of whether there is a polymorphic or divergent site. 

Scripts: Two scripts are contained with the vcf2MK package. detailed2site_averages.pl will take the output from vcf2MK and aggregate information for all tanscripts that contain a given site, and output the average dn, ds, pn, and ps counts for that site across all transcripts. The package also contains a script, detailed2mk.pl, that will aggregate site counts across each transcript and output totals. By default this script will discard all sites with minor allele frequencies less than 0.1. The output of this script is tab-delimited with the transcript id, dn, ds, pn, ps, dos, and alpha. 

