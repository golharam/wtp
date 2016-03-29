#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "zutil.h"
#include <time.h>

void matrixunit(FILE *f, int i) 
{
    fprintf(f, "Matrix %d\nDefault 0 10\nEnd\n", i);
}

void matrixzero(FILE *f, int i) 
{
    fprintf(f, "Matrix %d\nDefault 0 0\nEnd\n", i);
}

void convert(char *line, char *line1, int len)
{
    int i;
    for (i = 0; i < len; i++) line1[i] = '1';
    line1[i] = 0;
    char *p = line;
    int a;
    while (p) {
	a = atoi(p);
	line1[a-1] = '0';
	p = strchr(p, ' ');
	if (p) p++;
    }
    printf("%s\n", line1);
}

main(int argc, char *argv[])
{
    if (argc < 6) fatal("makepattern pat templatefile numMis outMatrix outTemplate\n");
    char *pat = argv[1];
    char *file = argv[2];
    int len = strlen(pat), zero = 0, i, j;
    int numMis = atoi(argv[3]);
    char *outMa = argv[4];
    char *outT = argv[5];

    char *p;
    for (p = pat+1; *p; p++) {
	if (*p == '0') zero++;
    }
    FILE *fp = ckopen(file, "r");
    char line[1000];
    int found = 0;
    while (fgets(line, sizeof line, fp)) {
	if (line[0] = 't') {
	    int a, b, c;
	    sscanf(line, "template %d %d %d", &a, &b, &c);
	    //printf("%d %d %d %d %d\n", a, b, len, zero, numMis);
	    if (a == len-zero-1 && b == numMis) {
		found = 1;
		break;
	    }
	}
    }
    FILE *fpp =  ckopen(outMa, "w");
    fprintf(fpp, "Header\nNumMatrices %d\nEnd\n", len);
    for (i = 0; i < len; i++) {
	if (pat[i] == '0') {
	    matrixzero(fpp, i+1);
	} else {
	    matrixunit(fpp, i+1);
	}
    }
    fclose(fpp);
    if (!found) {
        fatal("template DB does not have the right templates\n");
    }
    fpp = ckopen(outT, "w");
    char line1[1000];
    while (fgets(line, sizeof line, fp)) {
	if (line[0] == 't') break;
	convert(line, line1, len-zero-1);
	for (i = 1, j = 0; i < len; i++) {
	    if (pat[i]=='1') {
		fprintf(fpp, "%c", line1[j++]);
	    } else {
		fprintf(fpp, "0");
	    }
	}
	fprintf(fpp, "\n");
    }
    fclose(fp);
    fclose(fpp);
}


