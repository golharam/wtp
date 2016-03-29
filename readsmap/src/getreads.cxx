#include "util.h"
#include "zutil.h"

#define MAXSIZE 20000

class read {
  public:
    char *seq;
    char *def;
    read *next;
    int dir;
};

static void addread(char *seq, char *def, char *anchor, int *beg, int *end, read **reads, int n, int L, int H)
{
    char *p = anchor;
    int len = strcspn(def, "\n");
    char c = def[len];
    def[len]= 0;
    while (p = strchr(p, ',')) {
	p++;
	int pos = atoi(p);
	int i;
	for (i = 0; i < n; i++) {
	    if ((pos >= 0 && pos >= beg[i]+L && pos <= end[i]+H) ||
		(pos < 0 && pos >= -end[i]+L && pos <= -beg[i]+H)) {
		read *r = new read();
		r->seq= strsave(seq);
		r->def = strsave(def);
		r->next = reads[i];
		r->dir = (pos >= 0)? 1: -1; 
		reads[i] = r;
	    }
	}
    } 
    def[len] = c;
}

static void getreads(int a, int b, int c, char *readseq, char *def, FILE *fp)
{
    char line1[100000];
    do {
	if (fgets(line1, sizeof line1, fp)== NULL) fatal("Not find read\n");
	int x, y, z;
	sscanf(line1, ">%d_%d_%d", &x, &y, &z);
	if (a == x && b==y && c==z) {
	    strcpy(def, line1);
	    if (fgets(line1, sizeof line1, fp)== NULL) fatal("No read seq\n");
	    strcpy(readseq, line1);
	    return;
	}
    } while (1);
}

main(int argc, char *argv[])
{
    if (argc < 10) {
	fatal("getreads Rfile Ffile pairfile lowcovfile rlen Low High PATH colorfile mini_over\n");
    }
    int rlen = atoi(argv[5]);
    int L = atoi(argv[6]), H = atoi(argv[7]);
    char *Rfile = argv[1];
    char *Ffile = argv[2];
    char *pfile = argv[3];
    char *lowfile = argv[4];
    char *path = argv[8];
    char *colorfile=argv[9];
    int over = atoi(argv[10]);
    FILE *fp = ckopen(lowfile, "r");
    char line[100000];
    int beg[MAXSIZE], end[MAXSIZE];
    read *reads[MAXSIZE];
    int n = 0;
    while (fgets(line, sizeof line, fp)) {
	int f, r;
	sscanf(line, "%d %d", &f , &r);
	if (n >= MAXSIZE) fatal("Out of bound\n");
	beg[n] = f;
	end[n] = r;
	reads[n] = NULL;
	n++;
    }
    fclose(fp);
    fp = ckopen(pfile, "r");
    FILE *fpf = ckopen(Ffile, "r");
    FILE *fpr = ckopen(Rfile, "r");
    char readf[1000];
    char readr[1000];
    char deff[100000];
    char defr[100000];
    while (fgets(line, sizeof line, fp)) {
	if (strchr(line, ',')) continue;
	int a, b, c;
	sscanf(line, ">%d_%d_%d", &a, &b, &c);
	getreads(a, b, c, readf,  deff,  fpf);
	getreads(a, b, c, readr, defr, fpr);
	addread(readr, defr, deff, beg, end, reads, n, L, H);
	addread(readf, deff, defr, beg, end, reads,n, -H, -L);
    }
    fclose(fp);
    fclose(fpf);
    fclose(fpr);
    int i;
    for (i = 0; i < n; i++) {
	char command[10000];
	//printf("=%d %d\n", beg[i], end[i]);
	printf("%d %d\n", beg[i], end[i]);
	FILE *fp = ckopen("temp", "w");
	read *r;
	for (r = reads[i]; r; r = r->next) {
	    fprintf(fp, "%s", r->def);
	    if (r->dir == -1) fprintf(fp, "@");
	    fprintf(fp, "\n%s", r->seq);
	}
	fclose(fp);
	sprintf(command, "%s/miniassem %d %d temp %s %d %d", path, beg[i]-rlen*2, end[i]+rlen*2, colorfile, rlen-1, over);
	printf("%s\n", command);
	system(command);
    }
}
