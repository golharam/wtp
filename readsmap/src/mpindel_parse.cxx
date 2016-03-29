
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "zcompress.h"

#define MAXSIZE 100000
#define MSIZE 20

int readin(int &n, genFile &gF_fp, char *id, char *readseq, int order)
{
    FILE *fp = gF_fp.getFILE();
    char line[1000000];
    do {
	if (fgets(line, sizeof line, fp) == NULL) return 0;
	if (line[0] != '>') continue;
	break;
    } while (1);
  
    if (fgets(readseq, 1000, fp) == NULL) return 0;
    return 1;
}

int lenofread(char *a)
{
    int i = strcspn(a, "\r\n")-1;
    return i;
}

int findReadLength(char *fasta_fn) {
  genFile gF_fasta;
  gF_fasta.setUsePipe();
  gF_fasta.open(fasta_fn,"r");
  char id[1000], read[1000];
  int n1;
  readin(n1, gF_fasta, id, read, 1);
  int read_length = lenofread(read); 
  return read_length;
}

int main(int argc, char *argv[])
{
    if (argc == 1) {fprintf(stderr, "mpindel_parse file LenR LenF [state]\nCan use r3:r3_csfasta, f3:f3_csfasta instead of LenR LenF\n"); exit(1); }

    genFile gF_fp;
    gF_fp.setUsePipe();
    gF_fp.open(argv[1],"r");
    FILE *fp = gF_fp.getFILE(); //fopen(argv[1], "r");

    int length_r,length_f;

    if(argv[2][0]=='r') {
      length_r = findReadLength(argv[2]+3);
      fprintf(stderr, "R3 read length detected as %d\n", length_r);
    } else {
      length_r = atoi(argv[2]);
    }
    if(argv[3][0]=='f') {
      length_f = findReadLength(argv[3]+3);
      fprintf(stderr, "F3 read length detected as %d\n", length_f);
    } else {
      length_f = atoi(argv[3]);
    }
    int state = 0; 
    if (argc >4) state = atoi(argv[4]); 

    char line[100000];
    while (fgets(line, sizeof line, fp)) {
	char *p = strchr(line, '|');
	if (!p) continue;
	char *q = strchr(line, '(');
	if (!q) continue;
	char *d = strchr(line,',');
	if (strchr(p+1, '|')) continue;
	char *x; 
	int length;
	if (q > p) { length  = length_f; x = p+1;}
	else {x = d+1; length = length_r;}
	int i, j, m, s, c = 0;
        if (state) {
	    sscanf(x, "%d_%d.%d.%d(%d", &c, &i, &j, &m, &s);
	} else {
	    sscanf(x, "%d.%d.%d(%d", &i, &j, &m, &s);
	}
	i += s;
	if (i < 0) {
	    i = -i;
	    i--;
	    if (j-length+1 >0) 
	    	i -= (j-length+1);
	}
	long long xx = (((long long) c) << 32)+i;
	j = length-1-j;
	printf("%lld\t%d\t%d\t%s",xx, j, m, line);
    }
    return 0;
}


