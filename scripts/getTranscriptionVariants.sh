allChrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"

for chr in $allChrs
do
  cat minus/locs.CAAT_G.$chr.txt minus/locs.CTAT_G.$chr.txt plus/locs._ATTGC.2.txt plus/locs._ATAGC.2.txt | \
  awk ' BEGIN { OFS="\t" } {  if ( $6 == "1" ) { print $0 } }' | \
  sort -k 2 -n | \
  awk ' BEGIN { OFS="\t" } {   print $1,$2,$2,$1":"$2 } ' > \
  transcriptVars.template.range.$chr.lst
done

for chr in $allChrs
do
  plink \
  --bed /SAN/ugi/UGIbiobank/data/downloaded/ukb23155_c${chr}_b0_v1.bed \
  --bim /SAN/ugi/UGIbiobank/data/downloaded/UKBexomeOQFE_chr${chr}.bim \
  --fam /SAN/ugi/UGIbiobank/data/downloaded/ukb23155_c22_b0_v1_s200632.fam \
  --extract range transcriptVars.template.range.$chr.lst \
  --recode A \
  --out transcriptVars.template.$chr
done

exit

for chr in $allChrs
do
  cat plus/locs.CAAT_G.$chr.txt plus/locs.CTAT_G.$chr.txt minus/locs._ATTGC.2.txt minus/locs._ATAGC.2.txt | \
  awk ' BEGIN { OFS="\t" } {  if ( $6 == "1" ) { print $0 } }' | \
  sort -k 2 -n | \
  awk ' BEGIN { OFS="\t" } {   print $1,$2,$2,$1":"$2 } ' > \
  transcriptVars.range.$chr.lst
done

for chr in $allChrs
do
  plink \
  --bed /SAN/ugi/UGIbiobank/data/downloaded/ukb23155_c${chr}_b0_v1.bed \
  --bim /SAN/ugi/UGIbiobank/data/downloaded/UKBexomeOQFE_chr${chr}.bim \
  --fam /SAN/ugi/UGIbiobank/data/downloaded/ukb23155_c22_b0_v1_s200632.fam \
  --extract range transcriptVars.range.$chr.lst \
  --recode A \
  --out transcriptVars.$chr
done
