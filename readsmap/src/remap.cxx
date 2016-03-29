#include "fasta-io.h"
#include "rescue.h"

main(int argc, char *argv[])
{
    char *refseq = argv[1];
    char *readfile = argv[2]; 
    char *qfile = argv[3];
    int thres =200;
    int beg=0, end=-1, chrom = 0;
    int gapalign = 0;
    int ins = 0, del = 0;
    int hi = 0, hd = 0;
    char *rfile  = NULL;
    int match_count = 0;
    char *masking = NULL;
    int filter = 0;
    int iomode = 0;
    float SNPrate = 0.0;
    int amballow = 0;  
    int adjas1 = 0;
    if (argc <= 2) { 
	fprintf(stderr, "remap(0.1) refseq readfile qfile [M=masking] [P=I/Omode]\n");
	exit(1);
    }	
    int i;
    for (i=3; i<argc; ++i)
    {
        if (argv[i][1] == '=')         // X=value
        {
            if (argv[i][2] == 0)
                fprintf(stderr,"Missing value for %s\n",argv[i]);
	    else if (argv[i][0] == 'M') 
		masking = argv[i]+2;
	    else if (argv[i][0] == 'P')
		iomode = atoi(argv[i]+2);
	    else if (argv[i][0] == 'U') {
                SNPrate = atof(argv[i]+2);
            } else if (argv[i][0] == 'H') {
                amballow = atoi(argv[i]+2);
            } else if (argv[i][0] == 'A') {
		adjas1 = atoi(argv[i]+2); 
	    }
	}
    }
                         
        FastaFile fafile(SEQTYPE_NT);
        FastaSeq faseq;
   if (iomode == 0) {
        if (!fafile.Open(refseq,"r") || !fafile.Read(faseq)) {
	    fatal("Can not open fasta file\n"); 
        } 
    }
        remapping *r;
        r = new remapping();
	r->setmatching(adjas1, SNPrate, amballow, qfile, 0);
/*
        if (adjas1 == 0 && SNPrate <= 0.00001 && amballow == 0) {
                mapping_machine_color *mp = new mapping_machine_color();
                r->mp = mp;
        } else {
                if (adjas1<=1) {
                    mapping_machine_color_adj *mp = new mapping_machine_color_adj();
                    if (adjas1 == 1) {
                        mp->setSNPrate(2.0, amballow);
                    } else {
                        mp->setSNPrate(SNPrate, amballow);
                    }
                    r->mp = mp;
                } else {
                    mapping_machine_color_adja *mp = new mapping_machine_color_adja();
                    r->mp = mp;
                }
       } 
*/
    if (iomode == 0) {
        r->ProcessSeq(faseq.Label(),faseq.Sequence());
    } else {
	r->ProcessFile(refseq);
    }
/*
	char *a = NULL;
	r->mp->setQfile(qfile);
	r->mp->setMatrix(a);
*/
	r->Fmask(masking);
    FILE *fp = ckopen(readfile, "r");
    FILE *fpq = ckopen(qfile, "r");
    char line[1000000], line1[10000];
    char readseq[1000];
    while (fgets(line, sizeof line, fp)) {
	if (line[0] == '#') continue;
	while (fgets(line1, sizeof line1, fpq) && line1[0] == '#') 
		;
	if (fgets(readseq, sizeof readseq, fp) && fgets(line1, sizeof line1, fpq )) {
	    char *p = strchr(line, ',');
	    if (p) *p = 0; 
	    if (!p && (beg > end && rfile == NULL)) continue;
	    char rp[1000];
	    strcpy(rp, readseq);
	    r->setprobe("F", readseq, strlen(readseq)-2);
	    char *q = strchr(line, '\n');
	    if (q) *q = 0;
	    print_line(line);
	    if (r->mp->quality_value(line1)) {
		if (p) print_line(p+1);
		int outsize = printf("\n%s", rp);
		check_printf;
		continue;
	    }
	    while (p){
		p++;
		int b;
		if (iomode == 0) sscanf(p, "%d", &b);
		else {
			sscanf(p, "%d_%d", &chrom, &b);
			chrom--;
		}
	 	int s = r->re_align(chrom, b);
		int outsize; 
		if (iomode == 0) outsize = printf(",%d.%d", b, s/10);
		else outsize = printf(",%d_%d.%d", chrom+1, b, s/10);   
		check_printf;
		p = strchr(p, ',');
	    }
	    int outsize = printf("\n%s", rp);
	    check_printf;
	} else break; 	
    }  
} 
