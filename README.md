ReadMe.txt

Code and data relating to:

Sequence effects on patterns of variation observed from whole-genome sequenced UK Biobank participants David Curtis d.curtis@ucl.ac.uk

Data consists of counts of variants obtained from 500,000 UK Biobank participants

Contents of counts directory

Note that variants are specified as pentanucleotide background followed by base which replaces central base So ACGTTA represents AC[G>A]TT, CGTCGC represents CG[T>C]CG, etc.

counts.chr.x.txt - variant counts chr is chromosome or total x = 0 - SNP variants with allele count >=3 x = 1 - variants with allele count =2 x = 2 - singleton variants with allele count =1

backgroundCounts.chr.txt - pentanucleotide counts in the reference sequence

frequencies.chr.x.txt - variant frequencies, variant count divided by count of reference sequence background

Contents of src directory

*.?pp - C++ source makeExe.sh - bash script to compile C++ executables, requires DLIB from https://dlib.net/ makeVariantPlots.20260214.R - code to produce plots
