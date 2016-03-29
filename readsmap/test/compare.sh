./fasta2match.pl -g genome1 -r testreads1 -t 20 -e 2 -d ./
diff testreads1.ma.20.2 testreads1.ma.20.2.gold
./fasta2match.pl -g genome2 -r testreads2 -t 20 -e 2 -d ./
diff testreads2.ma.20.2 testreads2.ma.20.2.gold
