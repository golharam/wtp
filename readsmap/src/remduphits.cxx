#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "zutil.h"

#include "zcompress.h"

#define MAX_SEQ 5000000
#define CHUNK 1024
#define CSHIFT 10
#define MAXCHUNK 10000000

class chrom {
 public:
    int cnum;
    int *chunk;
}; 

static int num[8];
static long long hits_b[8][10000];
static long long  start = 0;
static int mul = 0, chrom_start = 0;
static long long *offset;
static chrom *carray = NULL;

int get_chrom(long long bb, long long  &br)
{
    long long b = llabs(bb)+ start;
    chrom *x = carray+(b >>CSHIFT);
    int c = x->cnum;
    if (c < 0) {
	c = x->chunk[b & ((long long) (CHUNK-1))];
    }
    b -= offset[c];
    if (bb >0) br = b;
    else br = -b;
    return c;
}

void ssset(chrom *s, int b, int e, int c)
{
    if (b == 0) {
	if (e == CHUNK) {s->cnum = c; return;}
	s->cnum = -1;
	s->chunk = new int[CHUNK];
    }
    int i;
    for (i = b; i < e; i++)
	s->chunk[i] = c;
} 

void setchrom(char *file)
{
    if (mul == 0) return;
    carray = new chrom[MAXCHUNK];
    offset = new long long[MAX_SEQ+1]; 
    FILE *fp =  fopen(file, "r");
    char a[10000];
    long long last = 0; 
    int i;
    int chunk_num = 0, index_num = 0, chrom_n = 1;
    while (fgets(a, 10000, fp)) {
	long long b = atoll(a);
	if (chrom_n > MAX_SEQ+1) fatal("Too many sequence\n");
	offset[chrom_n] = last; 
	if (b < last) {fprintf(stderr, "wrong index file\n"); exit(1);}
	int cn = b >> CSHIFT;
 	int num = b & ((long long)(CHUNK-1));
	if (cn == chunk_num) {
	    ssset(&carray[cn], index_num, num, chrom_n);
	} else {
	    ssset(carray+chunk_num, index_num, CHUNK, chrom_n);
	    for (i = chunk_num+1; i < cn; i++) 
		carray[i].cnum = chrom_n;
   	    ssset(carray+cn, 0, num, chrom_n);
	}
	chrom_n++;
	last = b;
	chunk_num = cn;
	index_num = num;
    } 
    fclose(fp);
}

int add(char *s)
{
    long long b;
    int f;
    char *tmp;
    sscanf(s, "%lld.%d", &b, &f);
    tmp = strchr(s, ':');
    int i;
    int h = (int) (llabs(b) & ((long long)7));
    for (i = 0; i < num[h]; i++) {
	if (b == hits_b[h][i]) return 0;
    }
    int outsize;
    if (mul == 0) {
	if (chrom_start == 0) f = f/10;
      if (b > 0) 
    	outsize = printf(",%lld.%d",  b + start, f);
      else outsize = printf(",%lld.%d", b - start, f);
    } else {
      long long	 bb;
      int ch = get_chrom(b, bb);
      outsize = printf(",%d_%lld.%d", ch+chrom_start,  bb, f/10);
    }
    check_printf;
    if (tmp) {
	char *p = strchr(s, ',');
	if (!p || p > tmp) { 
	    p = strchr(tmp, ',');
	    if (p) *p = 0;
	    outsize = printf("%s", tmp);
	    check_printf;
	    if (p) *p = ',';
	}
    }
    if (i >= 10000) {
	fprintf(stderr, "hash table full\n");
	exit(1);
    }
    hits_b[h][i] = b;
    num[h]++;
    return 1;
}

main(int argc, char *argv[])
{
    FILE *fp;
    if (strcmp(argv[3], "stdin")==0) fp= stdin;
    else fp = fopen(argv[3], "r");
    if (fp == NULL) exit(0); 
    genFile fp1;  //FILE *fp1 = NULL;
    if (argc > 7) fp1.open(argv[7], "r");
    int limit = atoi(argv[2]);
    int compact = 0;
    if (limit < 0) {limit *= -1; compact =1;} 
    char *indexfile = argv[6]; 
    start = atoll(argv[1]);
    mul = atoi(argv[4]);
    chrom_start = atoi(argv[5]);
    setchrom(indexfile);
    if (fp == NULL) return 1;
    char line[1000000], line1[10000], lastid[10000];
    int notprint = 0;
    while (fgets(line, sizeof line, fp)) {
	char *q = line;
	q = strchr(line ,'\n');
	if (q) *q = 0;
	q = strchr(line, ',');
	notprint = 0;
	if (compact && !q) notprint =1;
	int find = 1;
	int numhit = 0;
	if (!fp1.getIsNull() && fp1.gets(line1, sizeof line1)) {
	    while (line1[0] == '#' || line1[0] != '>') {
		if (fp1.gets(line1, sizeof line1) == NULL) {
		    find = 0;
		    break;
		}
	    }
	} else find = 0;
	if (notprint){
	    if (find) {
		strcpy(lastid, line1);
		fp1.gets(line1, sizeof line1);
	    } else {
		strcpy(lastid, line);
	    }
	    continue;
	}
	int outsize;
	if (find) {
	    char line2[1000];
	    sscanf(line1, ">%s", line2);
	    char *p = strchr(line2, ',');
	    if (p) *p = 0;
	    outsize = printf(">%s", line2);
	    check_printf;
	    fp1.gets(line1, sizeof line1);
	} else {
	    if (q) *q = 0;
	    char *qq = strchr(line, '\n');
	    if (qq) *qq = 0;
	    print_line(line);
	}
	if (q) {
	    memset(num, 0, 8*sizeof(int));
	    do {
		q++;
		numhit += add(q);
		if (numhit >= limit) break;
		q = strchr(q, ',');
	    } while (q);
	}
	printn();
	if (find && !compact) {
	    print_line(line1);
	}
    }
    if (notprint) print_line(lastid);
    if (fp != stdin) fclose(fp);
    if (!fp1.getIsNull()) fp1.close();
}
    
