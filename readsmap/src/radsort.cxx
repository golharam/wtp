#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "zutil.h"

#define MS 2000
#define RAD 32
#define RADBIT 5
#define MAXSIZE 1000000
#define PSIZE 30000000 

static char header[1000];

void openwrite(FILE **fp)
{
    int i;
    char file[1000];
    for (i = 0; i < RAD; i++) {
	sprintf(file, "out%d.%s", i, header);
	fp[i] = ckopen(file, "w");
    }
}

void closefiles(FILE **fp)
{
    int i;
    for (i = 0; i < RAD; i++) {
	fclose(fp[i]);
    }
}

int Read(FILE *fp, char *label, char *seq)
{
    do {
	if (fgets(label, MAXSIZE, fp) == NULL) return 0;
	if (label[0] == '>') break;
    } while (1);
    if (fgets(seq, MAXSIZE, fp) == NULL) return 0;
    return 1;
}

void movefile()
{
    int i;
    char file[1000];
    for (i = 0; i < RAD; i++) {
	sprintf(file, "mv -f out%d.%s in%d.%s", i ,header, i, header);
	system(file);
    }
}

int num(char *label, int &d)
{
    int a, b, c;
    sscanf(label, ">%d_%d_%d", &a, &b, &c);
    d = c;
    return (a<<11)+b-1;
}
    

void output(int i, FILE **fp, char *label, char *seq)
{
    int d;
    int id = num(label, d);
    //printf("id%d\n", id);
    int j = (id >> (i*5)) & 31;
    fprintf(fp[j], "%s%s", label, seq);
}

void readfile(int j, char *a)
{
    sprintf(a, "in%d.%s", j ,header);
}

int main (int argc, char **argv)
{
    if (argc < 3) 
	fatal("radsort readsfile numPass\n");
    char *seqfile = argv[1];
    int n = atoi(argv[2]);
    int i, j;
    if (n <= 0) return 1;
    if (n > 32/RADBIT+1) n = 32/RADBIT+1;

    char tmpfile[100]; 
    sprintf(tmpfile, "out0.XXXXXX");
    mktemp(tmpfile);
    strcpy(header, tmpfile+5);

    FILE *fp[RAD], *fr[RAD];
    openwrite(fp);
    FILE *faseq = ckopen(seqfile, "r");
    char Label[MAXSIZE], Seq[MAXSIZE];
    while (Read(faseq, Label, Seq))
    {
	output(0, fp, Label, Seq);
    }
    fclose(faseq);
    closefiles(fp);
    
    for (i = 1; i < n; i++) {
	movefile();
	openwrite(fp);
	for (j = 0; j < RAD; j++) {
	    char a[100];
	    readfile(j, a);
	    faseq = ckopen(a, "r");
	    while (Read(faseq, Label, Seq))
	    {
		output(i, fp, Label, Seq);
	    }
	    fclose(faseq);
	}
	closefiles(fp);
    }
    movefile();
    //sort the last number
    int lastn;
    char *la[MS];
    char se[MS][200];
    char *pool = new char[PSIZE];
    if (pool == NULL) fatal("Not enough memory\n");
    int used;
    int nn[MS];
    for (i = 0; i < RAD; i++) {
	int j = 0;
	lastn = -1;
	char a[100];
	readfile(i, a);
	FILE *fp = ckopen(a, "r");
	int r;
	used = 0;
	do {
	    r = Read(fp, Label, Seq);
	    int l;
	    int x;
	    if (r && (x=num(Label, l)) == lastn) {
		if (j >= MS) fatal("radsort error: MS size too small");
	    } else {
		// slow sort;
		int k, g, m, w;
		for (k = 0; k < j; k++) {
		    m = 100000;
		    for (g = 0; g < j; g++) {
			if (m > nn[g]) {
			    m = nn[g];
			    w = g;
			}
		    }
		    printf("%s%s", la[w], se[w]);
		    nn[w] = 100001;
		}
		j = 0;
		if (!r) break;
		lastn = x;
		used = 0;
	    }
	    la[j] = pool+used;
            used += strlen(Label)+1;
	    if (used > PSIZE) fatal("radsort error:pool memory over limit\n");
	    strcpy(la[j], Label);
	    strcpy(se[j], Seq);
	    nn[j] = l;
	    j++;
	} while (r);
	fclose(fp);
    }
    for (i = 0; i < RAD; i++) {
	char a[100];
	sprintf(a, "rm -f in%d.%s", i,  header);
	system(a);
    }
    return 0;
}
