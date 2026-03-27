README.md

Code and data relating to:

Sequence effects on mutation frequencies investigated in whole-genome sequenced UK Biobank participants

David Curtis 

d.curtis@ucl.ac.uk

Data consists of counts of variants obtained from 500,000 UK Biobank participants

Contents of counts directory

backgroundCounts.total.txt 

Pentanucleotide counts in the reference sequence

Note that variants are specified as pentanucleotide background followed by base which replaces central base So ACGTTA represents AC[G>A]TT, CGTCGC represents CG[T>C]CG, etc.

acBinCounts.total.txt 

The first column is the variant. The next column is the count of singleton variants. The next doubletons, then AC=3, AC=4, etc. The final column is AC>200 pooled.

acBinFrequencies.overall.txt  

Same as above but frequencies calculated as counts divided by background counts.

Contents of src directory

*.?pp  

C++ source files

makeExe.sh 

Bash script to compile C++ executables, requires DLIB from https://dlib.net/ 

makeVariantPlots.20260323.R 

Code to produce plots
