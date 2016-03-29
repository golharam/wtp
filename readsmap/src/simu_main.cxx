#include "simu.h"
#include <sys/time.h>
#include <time.h>

static int secondoftime()
{

    time_t clock, c1;
    struct tm *ttm;
    clock = time(&c1);
    ttm = localtime(&clock);
    char str[100];
    strftime(str, 100, "%s", ttm);
    return atoi(str);
}

main(int argc, char *argv[])
{
    if (argc < 7) {
        fprintf(stderr, "simu binfile nreads totalnum(distribution count reads) precision(int) maxlen seqfile [erate]\n");
        exit(1);
    }

    char *binfile = argv[1];
    int nreads = atoi(argv[2]);
    int totalnum = atoi(argv[3]);
    int prec = atoi(argv[4]);
    int maxlen = atoi(argv[5]);
    char *seqfile = argv[6];
    int nn = totalnum/prec+1;
    int sec = secondoftime();
    srand48(sec);

    if (argc > 7) erate = (double) atof(argv[7]);
    init_code(binfile, nn, prec, maxlen);

    FastaFile fafile(SEQTYPE_NT);

    if (!fafile.Open(seqfile,"r"))
        return 0;
    FastaSeq faseq;
    if (!fafile.Read(faseq)) return 0;
    int seqlen = faseq.Length();
    seqlen = firstcl(faseq.Sequence(), seqlen);
    int i;
    char forwardf[1000], reversedf[1000];
    sprintf(forwardf, "%s.for", seqfile);
    sprintf(reversedf, "%s.rev", seqfile); 
    FILE *ff = fopen(forwardf, "w");
    FILE *rf = fopen(reversedf, "w");
    for (i = 0; i < nreads; i++) {
        char f[100], r[100];
        int s1, s2;
        produce_read(faseq.Sequence(), seqlen, nn, r,f, s1, s2);
        fprintf(ff, ">1_1_%d_forward\n%s\n", i, f);
	fprintf(rf, ">1_1_%d_reversed\n%s\n", i, r);
    }
    fclose(ff);
    fclose(rf);
}

