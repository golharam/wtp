#include "simu.h"
int L=25;
double erate = 0.05;
char seq2color[128][128];
char Aadiffer = 'a'-'A';
int *Index, *length;

static void sets2c(char s[][128], char a, char b, char c)
{
    s[a][b] = s[a+Aadiffer][b] = s[a][b+Aadiffer]= s[a+Aadiffer][b+Aadiffer] = c
;
}

void mutate(char *seq, int len, double error_rate)
{
    int i;
    for (i = 0; i < len; i++) {
	if (drand48() > error_rate)  continue;
	double r = drand48()*4;
	int rr = (int) r;
	seq[i] = '0'+((seq[i]-'0'+rr)% 4);
    }
}

int firstcl(const char *seq, int length)
{
    char *s = strchr(seq, '.'); 
    if (s == NULL) return length;
    return s-seq;
}

void init_code(char *binfile, int &nn, int prec, int maxlen)
{
    Index  = new int[nn];
    length = new int[nn+1];
    int i;
    for (i = 0; i < 128; i++) {
        seq2color['.'][i] = '4';
    }
    sets2c(seq2color, 'A','A', '0');
    sets2c(seq2color, 'C','C', '0');
    sets2c(seq2color, 'T','T', '0');
    sets2c(seq2color, 'G','G', '0');
    sets2c(seq2color, 'A','C', '1');
    sets2c(seq2color, 'C','A', '1');
    sets2c(seq2color, 'G','T', '1');
    sets2c(seq2color, 'T','G', '1');
    sets2c(seq2color, 'A','G', '2');
    sets2c(seq2color, 'G','A', '2');
    sets2c(seq2color, 'C','T', '2');
    sets2c(seq2color, 'T','C', '2');
    sets2c(seq2color, 'A','T', '3');
    sets2c(seq2color, 'T','A', '3');
    sets2c(seq2color, 'C','G', '3');
    sets2c(seq2color, 'G','C', '3');

    FILE *fp = fopen(binfile, "r");
    if (fp == NULL) exit(1);

    char line[1000];
    int cur = 0,  num = 0, j = 0;
	i = 0;
    Index[i]= 0;
    cur += prec;
    i++;
    while (fgets(line, sizeof line, fp)) {
        int a, b;
        sscanf(line, "%d %d", &a, &b);
	if (a > maxlen) break;
        num += b;
        while (cur < num) {
            Index[i] = j;
            cur += prec;
            i++;
        }
        length[j] = a;
        //printf("%d %d\n", j, a);
        j++;
    }
    length[j] = maxlen;
    nn = i;
}

int get_distance(int nn)
{
    double r = drand48();
    int rd = r * nn;
    int j = Index[rd];
    r = drand48()*(length[j+1] -length[j]);
    int dis = length[j]+r;
    return dis;
}

void makeseq(char *f, int s, int len, const char *seq)
{
    f[0] = seq[s];
    int i;
    for (i = 0; i < len; i++) {
	f[i+1] = seq2color[seq[s+i]][seq[s+i+1]];
    }
    f[i+1] = 0;
}

void makeseqwdel(char *f, int s, int len, const char *seq, int sdel, int ldel)
{
    f[0] = seq[s];
    int i;
    for (i = 0; i < sdel; i++) {
        f[i+1] = seq2color[seq[s+i]][seq[s+i+1]];
    }
    f[i+1] = seq2color[seq[s+i]][seq[s+i+ldel+1]];
    i++;
    for (;i < len; i++) {
	f[i+1] = seq2color[seq[s+i+ldel]][seq[s+i+1+ldel]];
    }
    f[i+1] = 0;
} 
void makeseqwins(char *f, int s, int len, const char *seq, int sdel, int ldel)
{
    f[0] = seq[s];
    char acode[] = "0123";
    int i, j;
    for (i = 0; i < sdel; i++) {
        f[i+1] = seq2color[seq[s+i]][seq[s+i+1]];
    }
    for (j = 0; j < ldel; j++) {
	f[i+1] = acode[(int)(drand48()*4)];
	i++;
    }
    for (;i < len; i++) {
        f[i+1] = seq2color[seq[s+i-ldel]][seq[s+i+1-ldel]];
    }
    f[i+1] = 0;
}

void produce_read(const char *seq, int seqlen, int nn, char *f, char *r, int &s1, int &s2)
{
    int dis = get_distance(nn);
    s1 = drand48()*(seqlen-dis-L*2);
    s2 = s1+L+dis;
    makeseq(f, s1, L,  seq);
    makeseq(r, s2, L,seq);
    mutate(f+1, L, erate);
    mutate(r+1, L, erate);
}

void produce_del_read(const char *seq, int seqlen, int nn, char *f, char *r, int &s1, int &s2, int maxdel, int ends)
{
    int dis = get_distance(nn);
    s1 = drand48()*(seqlen-dis-L*2);
    s2 = s1+L+dis;
    double a = drand48()*2;
    int dellen = drand48()*maxdel+1;
    int middle = drand48()*(L-ends-ends)+ends; 
    if (a > 1) {
	makeseq(f, s1, L,  seq);
	makeseqwdel(r, s2, L, seq,  middle, dellen);
    } else {
        makeseq(r, s2, L,  seq);
        makeseqwdel(f, s1, L, seq,  middle, dellen);
    }
    mutate(f+1, L, erate);
    mutate(r+1, L, erate);
}
void produce_ins_read(const char *seq, int seqlen, int nn, char *f, char *r, int &s1, int &s2, int maxdel, int ends)
{
    int dis = get_distance(nn);
    s1 = drand48()*(seqlen-dis-L*2);
    s2 = s1+L+dis;
    double a = drand48()*2;
    int dellen = drand48()*maxdel+1;
    int middle = drand48()*(L-ends-ends)+ends;
    if (a > 1) {
        makeseq(f, s1, L,  seq);
        makeseqwins(r, s2, L, seq,  middle, dellen);
    } else {
        makeseq(r, s2, L,  seq);
        makeseqwins(f, s1, L, seq,  middle, dellen);
    }
    mutate(f+1, L, erate);
    mutate(r+1, L, erate);
}

