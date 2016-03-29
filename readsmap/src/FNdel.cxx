#include "simu.h"

main(int argc, char *argv[])
{

    if (argc < 10) {
        fprintf(stderr, "simu binfile nreads totalnum(distribution count reads)precision(int) maxlen seqfile maxindellen minendlen inordel(0=del,1=ins)\n");
        exit(1);
    }

    char *binfile = argv[1];
    int nreads = atoi(argv[2]);
    int totalnum = atoi(argv[3]);
    int prec = atoi(argv[4]);
    int maxlen = atoi(argv[5]);
    char *seqfile = argv[6];
    int nn = totalnum/prec+1;
    int dellen = atoi(argv[7]); 
    int ends = atoi(argv[8]);
    int iord = atoi(argv[9]);

    init_code(binfile, nn, prec, maxlen);

    FastaFile fafile(SEQTYPE_NT);

    if (!fafile.Open(seqfile,"r"))
        return 0;
    FastaSeq faseq;
    if (!fafile.Read(faseq)) return 0;
    int seqlen = faseq.Length();
    seqlen = firstcl(faseq.Sequence(), seqlen);
    char forwardf[1000], reversedf[1000];
    sprintf(forwardf, "%s.for", seqfile);
    sprintf(reversedf, "%s.rev", seqfile);
    FILE *ff = fopen(forwardf, "w");
    FILE *rf = fopen(reversedf, "w");
    int i;
    for (i = 0; i < nreads; i++) {
        char f[100], r[100];
        int s1, s2;
	if (iord==0)
            produce_del_read(faseq.Sequence(), seqlen, nn, f,r, s1, s2, dellen, ends);
	else 
	    produce_ins_read(faseq.Sequence(), seqlen, nn, f,r, s1, s2, dellen, ends);
        fprintf(ff, ">1_1_%d_%d_forward\n%s\n", i, s1, f);
        fprintf(rf, ">1_1_%d_%d_reversed\n%s\n", i, s2, r);
    }
    fclose(ff);
    fclose(rf);
}

