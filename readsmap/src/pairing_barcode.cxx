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
    if (argc<2) {
	fprintf(stderr, "pairing_barcode(V0.1) -r R3file -f F3file\n ");
	exit(1);
    } 
    char *file1 = NULL;
    char *file2 = NULL;

    int i, ins = 0, del = 0, hi = 0, hd = 0;
    printf("#");
    for (i = 0; i <argc; i++) {
	printf("%s ", argv[i]);
    }
    printf("\n");
    for (i=1; i<argc; ++i)
    {
	if (argv[i][0] == '-')         // X=value
	{
	    if (argv[i][1] == 0)
		fprintf(stderr,"Missing flag after -\n");
	    else if (argv[i][1] == 'f') {
		i++;
		if (i < argc) file2 = argv[i]; //discont_input = 0;
	    } else if (argv[i][1] == 'r')  { 
		i++;
		if (i < argc) file1 = argv[i];
	    }
	}
    }
    char id[100000], id1[100000], read1[1000], read2[1000];
    int n1, n2;
    FILE *fp1 = ckopen(file1, "r");
    FILE *fp2 = ckopen(file2, "r");

    readin(fp1, id, read1, 1); 
    readin(fp2, id1, read2, 2);
    int order = 0;
    do {
	int a = compare_id(id, id1);
	if (a < 0) {
	    if (!readin(fp1, id, read1, order)) break;
	} else if (a > 0) {
	    if (!readin(fp2, id1, read2, order)) break;
	} else {
	    printf("%s%s%s%s", id1, read2, id, read1); 
            if (!readin(fp1, id, read1, order)) break;
            if (!readin(fp2, id1, read2, order)) break;
	}
    } while(1);
    fclose(fp1);
    fclose(fp2);
}

