#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "zutil.h"
#include <math.h>

#define MAXLEN 10
static int rlen = 0;
static double miss[MAXLEN];
static double A = 1.0;

static void ana(char *line)
{
    double x = 0;
    char *p = strchr(line, ',');
    if (!p) {
	p = strchr(line, '\n');
	*p = 0;
    } else {
	*p = 0;
	do {
	    p++;
	    int n;
	    sscanf(p, "%*d.%d", &n);
	    //n = n/10;
	    x += 1.0/miss[n];
	} while (p=strchr(p, ','));
    }
    printf("%s,%f\n", line, (float) (x-1.0/miss[0]));
}

static void calmiss()
{
    int i;
    int j = 1;
    int x = rlen;
    double p = 3.0;
    if (x > MAXLEN) x = MAXLEN;
    miss[0] = A;
    for (i = 1; i < x; i++) {
	j *= rlen-i+1;
	j /= i;
	miss[i] = p*j*A;
	p *= 3;
	//printf("%f\n", miss[i]);
    }
}


int main (int argc, char **argv)
{
    if (argc < 3) 
	fatal("seqstat file readlen\n");
    char *rfile = argv[1];
    rlen = atoi(argv[2]);
	A=1.0;
    FILE *fp = ckopen(rfile, "r");
    if (rlen <=0) fatal("length nonpositive\n");
    calmiss();
    char line[500000];
    while (fgets(line, sizeof line, fp)) {
	ana(line);
    }
}
