#$1 is R3 fire $2 F3 file $3,$4 lower and upper limit $5 temp pairing filr $6 reference sequence $7 cutoff for unique starting point  
~/readsmap/pairing $1 $2 $3 $4 V=$5 G=2 D=10 T=2 I=3 i=7 d=7 S=$6 P=1 z=10 >indel.res 
#~/readsmap/pairing R3.25.3 F3.25.3 800 3500 V=pairing.dat S=ref.fa G=2 D=10 T=2 I=2 i=7 d=7
~/readsmap/mpindel_parse indel.res 1 > indel.pas
sort -n indel.ps > indel.sort
~/readsmap/mpindel_summ indel.sort $7 1 > indel.sum$7
