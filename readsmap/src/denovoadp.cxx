#include "util.h"
#include "zutil.h"
#define ChrPerLine 70

class Color2Seq {
  public:
    Color2Seq(int alen);
    ~Color2Seq();
    int  inputReads(FILE *file, int cost);
    int inputBaseReads(FILE *file, int cost);
    void outputSeq();
  protected:
    int assemlen;
    int **Co, **Ba;
    char *seq, *cseq, **color, *notN;
    int **pa;
    char color2seq[128][128];
    char code[128];
    int first, last;
    int numR;
}; 

Color2Seq::Color2Seq(int alen)
{
    assemlen = alen;
    first = alen+1;
    last = 0;	
    Co = new int*[3*assemlen];
    Ba = Co+assemlen;
    pa = Ba+assemlen;
    seq = new char[assemlen*2+2];
    cseq = seq+assemlen+1;
    Co[0] = new int[4*3*assemlen];
    memset(Co[0], 0, sizeof(int)*4*3*assemlen);
    Ba[0] = Co[0]+4*assemlen;
    pa[0] = Ba[0]+4*assemlen; 
    color = new char*[assemlen];
    color[0] = new char[5*assemlen];
    notN = color[0]+4*assemlen;
    memset(color[0], 0, 5*assemlen);
    int i;
    for (i = 1; i < assemlen; i++) {
	Co[i] = Co[i-1]+4; 
  	Ba[i] = Ba[i-1]+4;
	pa[i] = pa[i-1]+4;
	color[i] = color[i-1]+4;
    }

    color2seq['A']['0'] = 'A';
    color2seq['A']['1'] = 'C';
    color2seq['A']['2'] = 'G';
    color2seq['A']['3'] = 'T';
    color2seq['C']['0'] = 'C';
    color2seq['C']['1'] = 'A';
    color2seq['C']['2'] = 'T';
    color2seq['C']['3'] = 'G';
    color2seq['G']['0'] = 'G';
    color2seq['G']['1'] = 'T';
    color2seq['G']['2'] = 'A';
    color2seq['G']['3'] = 'C';
    color2seq['T']['0'] = 'T';
    color2seq['T']['1'] = 'G';
    color2seq['T']['2'] = 'C';
    color2seq['T']['3'] = 'A';
    init_scode(code);
}

Color2Seq::~Color2Seq()
{
    delete [] Co[0];
    delete [] seq;
    delete [] Co;
    delete [] color[0];
    delete color;
}

void Color2Seq::outputSeq()
{
    int i;
    int CinLine=0;
    if (numR == 0) return;
    if (last >= assemlen) fatal("contig longer than maximum\n");
    pa[last+1][0] = pa[last+1][1] =pa[last+1][2]=pa[last+1][3]=0;
    char c[] = "ACGT";
    for (i = last; i >= first; i--) {
	int j;
	for (j = 0; j < 4; j++) {
	    char B = c[j];
	    char N;
	    int k, min=500000000, index = 0;
	    for (k = 0; k < 4; k++) {
		N = color2seq[B]['0'+k];
		int x = code[N];
		int sc = pa[i+1][x]+Co[i][k];
		if (k==0 || sc < min){ min = sc; index =x;} 
	    }  
	    pa[i][j] = min+Ba[i][j];
	    color[i][j] = index;
	}
	if (pa[i][0] > 10000000) {
	    pa[i][1] -= pa[i][0];
	    pa[i][2] -= pa[i][0];
            pa[i][3] -= pa[i][0];
    	    pa[i][0] = 0;
	}	    
    }
    int min = pa[first][0], w = 0;
    for (i = 1; i < 4; i++) {
	if (min > pa[first][i]) { min = pa[first][i]; w = i;}
    }
    printf("%c", c[w]); CinLine = 1;
    for (i = first+1; i <= last; i++) {
	w = color[i-1][w];
        if (CinLine >= ChrPerLine) {
            printf("\n");
	    CinLine=0;
	}
	printf("%c", (notN[i])? c[w]:'N');
	CinLine++;
    }
    printf("\n");
}

int Color2Seq::inputBaseReads(FILE *fp, int cost)
{
    if (fp == NULL) return 0;
    char line[10000], line1[1000];
    int y = 0;
    //numR = 0;
    while (fgets(line, sizeof line, fp)) {
        if (line[0] != '=') continue;
        //printf(">%s", line+1);
        y = 1;
        break;
    }
    if (y == 0) return 0;
    while (fgets(line, sizeof line, fp)) {
        if (line[0] == '#') continue;
        if (line[0]== '-') {
            break;
        }
        if (line[0] != '>') fatal("format of input wrong\n");
        if (!fgets(line1, sizeof line1, fp)) break;
        char *c = strchr(line, ',');
        if (!c) continue;
        int pos = atoi(c+1);
        int p = abs(pos);
        int len = strcspn(line1, "\n\r");
        if (pos >= 0) {
            if (pos < first) first = pos;
            if (pos+len-1>last) last = pos+len-1;
        } else {
            if (-pos > last) last = -pos;
            if (-pos-len+1 < first) first = -pos-len+1;
        }
	if (first < 0 || last >= assemlen) 
                fatal("contig exceed maximum length or reads match to negative position\n");

	int i;
	for (i = 0; i < len; i++) {
	    if (pos >= 0) { 
		p = pos+i;
	    } else {
		p = -pos-i;
	    }
	    Ba[p][0]+=cost; Ba[p][1]+=cost; Ba[p][2]+=cost; Ba[p][3]+=cost;
	    int x = code[line1[i]];
	    if (pos < 0) x = 3-x;
	    Ba[p][x]-=cost;
	}
	numR++;
    }
    return 1;
}


int Color2Seq::inputReads(FILE *fp, int cost)
{
    char line[10000], line1[1000];
    int y = 0;
    numR = 0;
    while (fgets(line, sizeof line, fp)) {
	if (line[0] != '=') continue;
	printf(">%s", line+1);
	y = 1;
	break;
    }
    if (y == 0) return 0;
    while (fgets(line, sizeof line, fp)) {
	if (line[0] == '#') continue;
	if (line[0]== '-') {
	    break;
	}
        if (line[0] != '>') fatal("format of input wrong\n");
        if (!fgets(line1, sizeof line1, fp)) break;
	char *c = strchr(line, ',');
	if (!c) continue;
	int pos = atoi(c+1);
	int p = abs(pos);
        int len = strcspn(line1, "\n\r");
        if (pos >= 0) {
	    if (pos < first) first = pos;
	    if (pos+len-3>last) last = pos+len-3; 
	} else {
	    if (-pos > last) last = -pos;
	    if (-pos-len+2 < first) first = -pos-len+2;
	}
	Ba[p][0]+=cost; Ba[p][1]+=cost; Ba[p][2]+=cost; Ba[p][3]+=cost;
	int x = code[color2seq[line1[0]][line1[1]]];
	if (pos < 0) x = 3-x;
	//printf("Base %d\n", x);
	Ba[p][x]-=cost; 
	notN[p] = 1;
	int i;
	for (i = 2; i < len; i++) {
	    if (pos >= 0) { 
		p = pos +i-2;
	    } else {
		p = -pos-i+1;
	    }
	    if (p < assemlen  && p >= 0) {
		Co[p][0]+=cost; Co[p][1]+=cost; Co[p][2]+=cost; Co[p][3]+=cost;
		Co[p][line1[i]-'0']-=cost;
		//printf("%d", line1[i]-'0');
	    } else {
		fatal("contig exceed maximum length or reads match to negative position\n");
	    }
	    notN[p] = 1;
	}
	//printf("\n");
	numR++;
    }
    return 1;
}


main(int argc, char *argv[])
{
    if (argc < 2) {
	fatal("denovoadp (V2.0 Nov 19,2008) readsfile max_assem_length [B=base_reads][C=##][M=##]\n Here C is the penalty for color, M is the penalty for base in a base read, defaults for both 1\n");
    }
    int assemlen = atoi(argv[2]);
    char *rfile = argv[1];
    FILE *fp = ckopen(rfile, "r"); 
    char *base_file = NULL;
    FILE *fpp = NULL;
    int i;
    int colorCost = 1, BaseCost = 1; 
    for (i = 3; i < argc; i++) 
	if (argv[i][1] == '=') {
	    if (argv[i][0] == 'B') {
		base_file = argv[i]+2;
	    } else if (argv[i][0] == 'C') {
		colorCost = atoi(argv[i]+2);
	    } else if (argv[i][0] == 'M') {
		BaseCost = atoi(argv[i]+2);
	    }
	}
    if (base_file ) fpp =ckopen(base_file, "r");
    do { 
    	Color2Seq *g = new Color2Seq(assemlen);
	int n;
    	if ((n=g->inputReads(fp, colorCost))==0) break; 
	g->inputBaseReads(fpp, BaseCost);
    	g->outputSeq();
	delete g;
    } while (1);
}
