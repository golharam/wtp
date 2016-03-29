#include "fasta-io.h"
#include "rescue.h"
#include "zcompress.h"
#define MAXREGION 5000000

main(int argc, char *argv[])
{
    char *refseq = argv[1];
    char *readfile = argv[2]; 
    int thres =200;
    int beg=0, end=-1, chrom = 0;
    int gapalign = 0;
    int ins = 0, del = 0;
    int hi = 0, hd = 0;
    char *rfile  = NULL;
    int match_count = 0;
    char *masking = NULL;
    float SNPrate = 0.0;
    int amballow = 0;
    int adjas1 = 0;
    char *qfile = NULL;
    int filter = 0;
    int iomode = 0;
    int outputSeq = 0;
    int numRe = 0;
    if (argc <= 2) { 
	fprintf(stderr, "rescue(1.6 5/29/2009) refseq readfile [T=threshold(200)] [G=gapalign(0)] [R=beg_end][I=ins][i=hi][D=del][d=hd][F=regionfile][f=0/1][M=mask][P=0/1][U=snprate][H=0/1][A=0/1][Q=qfile][S=0/1]\n");
	fprintf(stderr, "when G=2 allow 1 indel, use I, D to limit the maximum size of insertion and deletion resp, and d and i to limit minimum half diagnol lengths.\n");
	fprintf(stderr, "P=0 input/output position.mis; P=1 i/o chrom_pos.mis\n");
	fprintf(stderr, "A=0/1 whether to treat VA as 1\n");
	fprintf(stderr, "H=0/1  whether to allow IUB mapping\n");
	fprintf(stderr, "S=0/1 whether output read seq\n");
	exit(1);
    }	
    int i;
    for (i=3; i<argc; ++i)
    {
        if (argv[i][1] == '=')         // X=value
        {
            if (argv[i][2] == 0)
                fprintf(stderr,"Missing value for %s\n",argv[i]);
            else if (argv[i][0] == 'T')
                thres = atoi(argv[i]+2);
            else if (argv[i][0] == 'G')
		gapalign = atoi(argv[i]+2);
	    else if (argv[i][0] == 'R')
		sscanf(argv[i]+2, "%d_%d", &beg, &end);
	    else if (argv[i][0] == 'I') 
		ins = atoi(argv[i]+2);
	    else if (argv[i][0] == 'D')
                del = atoi(argv[i]+2);
	    else if (argv[i][0] == 'd')
		hd = atoi(argv[i]+2);
            else if (argv[i][0] == 'i')
                hi = atoi(argv[i]+2);
	    else if (argv[i][0] == 'F') 
		rfile = argv[i]+2;
            else if (argv[i][0] == 'f')
		filter = atoi(argv[i]+2);
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
            } else if (argv[i][0] == 'Q') {  
		qfile = argv[i]+2;
	    }  else if (argv[i][0] == 'S') {
		outputSeq = atoi(argv[i]+2);
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
        rescue *r;
        if (gapalign == 1) {
            r = new rescue();
 	} else if (gapalign == 2  && (del >0 || ins > 0)) {
	    rescue_oneindel *r1 = new rescue_oneindel();
	    r1->setMaxIndel(ins, del, SSIZE, thres, hi, hd); 
	    r = r1;
        } else {
            rescue_nogap *r1 = new rescue_nogap();
            r = r1;
        }
    r->setmatching(adjas1, SNPrate, amballow, qfile, 0);
    if (iomode == 0) {
        r->ProcessSeq(faseq.Label(),faseq.Sequence());
    } else {
	r->ProcessFile(refseq);
    }
        r->setMatrix(NULL);
	r->Fmask(masking);
	int scale = 10;
	if (qfile) scale = 100;
	r->setThreshold(thres, scale);
    thres *= 10;
    genFile gFp; gFp.setUsePipe();
    gFp.open(readfile,"r");
    FILE *fp = gFp.getFILE();  //ckopen(readfile, "r");
    char line[1000000];
    char readseq[1000];
    int *regionC, *regionB = NULL, *regionE; 
    while (fgets(line, sizeof line, fp)) {
	if (line[0] == '#') continue;

	char readseqNoCR[1002];
	{
	  readseqNoCR[0] = '[';
	  strcpy ( readseqNoCR+1, readseq );
	  size_t CRPos = strlen(readseq);
	  readseqNoCR[CRPos] = ']';
	} // Section added by EFT

	if (fgets(readseq, sizeof readseq, fp)) {
	    char *p = strchr(line, ',');
	    if (p) *p = 0; 
	    if (!p && (beg > end && rfile == NULL)) continue;
	    r->setprobe("F", readseq, strlen(readseq)-2);
	    char *q = strchr(line, '\n');
	    if (q) *q = 0;
	    printf("%s", line);
	    if (beg < end) {
		int le, ri;
		int s1 = r->align(chrom, beg, end, 10, le, ri);
		int s, m;
                if (r->check_score(s1, m, s)) {
		    printf(",%d.%d.%d", le, ri-le+1, s/10);
		    if (m > 0) printf("(%d:%d_%d)%s", m,r->get_range_beg(),
				      r->get_range_end(),
				      readseqNoCR);
		}
		printn();
		if (outputSeq) printf("%s\n", readseq);
		continue;
	    } 
	    if (rfile) {
		if (filter && p) continue;
		if (regionB == NULL) {
		    if (iomode == 0) regionC=new int[MAXREGION];
		    regionB = new int[MAXREGION];
		    regionE = new int[MAXREGION];
		    FILE *fp = fopen(rfile, "r");
		    int b, e, le,ri;
		    char line[10000];
		    numRe = 0; 
		    while (fgets(line, sizeof line, fp)) {
			if (numRe >= MAXREGION) fatal("Too many regions in region File\n");
		    	if (iomode == 0) {sscanf(line, "%d %d", regionB+numRe, regionE+numRe);}
		    	else sscanf(line, "%d %d %d", regionC+numRe, regionB+numRe, regionE+numRe);
			if (regionB[numRe]>regionB[numRe]) continue;
			numRe++;
		    }
		    fclose(fp);
		}
		int i;
		for (i = 0; i < numRe; i++) {
		    chrom = 1;
		    int b, e, le,ri;
		    if (iomode == 1) chrom = regionC[i];
		    b = regionB[i]; e = regionE[i];
		    int s1 =r->align(chrom-1, b, e, 10, le, ri);
		    int s, m;
                    if (r->check_score(s1, m, s)) {
			if (iomode == 0) printf(",%d.%d.%d", le, ri-le+1, s/10);
			else printf(",%d_%d.%d.%d", chrom, le, ri-le+1, s/10);
			if (m > 0) printf("(%d:%d_%d)%s", m,r->get_range_beg(),
				      r->get_range_end(),
				      readseqNoCR);
			match_count++;
		    }
		}
		printn();
                if (outputSeq) printf("%s\n", readseq);
		continue;
	    }

	    while (p){
		p++;
		int b, e, le, ri;
		if (iomode == 0) sscanf(p, "%d_%d", &b, &e);
		else {
			sscanf(p, "%d_%d_%d", &chrom, &b, &e);
			chrom--;
		}
		if (b <= e) {
	 	    int s1 = r->align(chrom, b, e, 10, le, ri);
		    int m, s;
		    if (r->check_score(s1, m, s)) {
			if (iomode == 0) printf(",%d.%d.%d", le, ri-le+1, s/10);
			else printf(",%d_%d.%d.%d", chrom+1, le, ri-le+1, s/10);
			if (m > 0) printf("(%d:%d_%d)%s", m,r->get_range_beg(),
					  r->get_range_end(),
					  readseqNoCR);
		    }
		}
		p = strchr(p, ',');
	    }
	    printn();
            if (outputSeq) printf("%s\n", readseq);
	} else break; 	
    }  
    if (rfile) fprintf(stderr, "match_count=%d\n", match_count);
} 
