#include "fasta-io.h"
#include "util.h"
#include "zutil.h"

int readin(FILE *fp, char *id, char *readseq, int order)
{
    char line[1000000];
    do {
	if (fgets(line, sizeof line, fp) == NULL) return 0;
	if (line[0] != '#') break;
	if (order == 2) printf("%s", line);
    } while (1);
    char *p = line;
    if (!order) {
      int a = compare_id(id, line);
      if (a > 0) {
	fprintf(stderr, "At least one of the read fine is out of order\n%s%s\n", id, line);
	exit(1);
      }
    }
    strcpy(id, line);
    if (fgets(readseq, 1000, fp) == NULL) return 0;
    return 1;
}

main(int argc, char **argv)
{
    char *file1 = NULL;
    char *file2 = NULL;
    char *pfile = NULL;
    char *misfile1 = NULL;
    char *misfile2 = NULL;

    int i;
    char cmt[10000];
    strcpy(cmt, "#");
    for (i = 0; i <argc; i++) {
	strcat(cmt, argv[i]);
	strcat(cmt, " ");
    }
    strcat(cmt,"\n");
    for (i=1; i<argc; ++i)
    {
	if (argv[i][0] == '-')         // X=value
	{
	    if (argv[i][1] == 0)
		fprintf(stderr,"Missing flag after -\n");
	    else if (argv[i][1] == 'R') {
		i++;
		if (i < argc) file2 = argv[i]; //discont_input = 0;
	    } else if (argv[i][1] == 'F')  { 
		i++;
		if (i < argc) file1 = argv[i];
	    } else if (argv[i][1] == 'P') {
		i++;
		if (i < argc) pfile = argv[i]; 
	    } else if (argv[i][1] == 'r') {
		i++;
		if (i < argc) misfile2 = argv[i]; 
	    } else if (argv[i][1] == 'f') {
		i++;
		if (i < argc) misfile1 = argv[i]; 
	    }
	}
    }
    if ( !file1 or !file2 or !pfile ) {
	fprintf(stderr, "\nprepairing (V0.2) -R R3file -F F3file -P Pairfile [-r R3unpaired] [-f F3unpaired]\n\n");
	exit(1);
    } 
    char id1[100000], id2[100000], read1[1000], read2[1000];
    FILE *fp1 = ckopen(file1, "r");
    FILE *fp2 = ckopen(file2, "r");
    FILE *pp  = NULL;
    pp= ckopen(pfile, "w");
    fprintf(pp, "%s", cmt);
    FILE *mp1 = NULL;
    FILE *mp2 = NULL;
    if (misfile1) mp1 = ckopen(misfile1, "w");
    if (misfile2) mp2 = ckopen(misfile2, "w");

    readin(fp1, id1, read1, 1); 
    readin(fp2, id2, read2, 2);
    int order = 0;
    unsigned long int count = 0;
    do {
	int a = compare_id(id1, id2);
	if (a < 0) {
	    if ( misfile1 ) {
	        fprintf(mp1, "%s%s", id1, read1);
	    } else {
	        fprintf(pp, "%s%s", id1, read1);
	        count++;
	    }
	    if (!readin(fp1, id1, read1, order)) break;
	} else if (a > 0) {
	    if ( misfile2 ) {
	        fprintf(mp2, "%s%s", id2, read2);
	    } else {
	        fprintf(pp, "%s%s", id2, read2);
	        count++;
	    }
	    if (!readin(fp2, id2, read2, order)) break;
	} else {
	    if ( pfile ) {
	        fprintf(pp, "%s%s%s%s", id2, read2, id1, read1);
	        count++;
	    } else {
                printf("%s%s%s%s", id2, read2, id1, read1);
            }
            if (!readin(fp1, id1, read1, order)) break;
            if (!readin(fp2, id2, read2, order)) break;
	}
    } while(1);
    fclose(fp1);
    fclose(fp2);
    if (pfile) fclose(pp);
    if (misfile1) fclose(mp1);
    if (misfile2) fclose(mp2);
    printf("Total output %d reads in %s\n", count, pfile);
    exit(0);
}
