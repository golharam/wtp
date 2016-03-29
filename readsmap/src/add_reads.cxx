#include "util.h"
#include "zutil.h"

main(int argc, char *argv[])
{
    if (argc < 3) {
	fatal("add_reads readsfile match_header_file state\nstate=0 include only reads in matchfile, state=1 include everything in readsfile\nThe two files must both be in correct order\n");
    }
    FILE *fpread = ckopen(argv[1], "r");
    FILE *fm = ckopen(argv[2], "r");
    int state = atoi(argv[3]);
    char line[1000000], line2[1000000];
    while (fgets(line, sizeof line, fm)) {
	if (line[0] != '>') continue;
	int a = -1;
	while (fgets(line2, sizeof line2, fpread)) {
	    if (line2[0] != '>') continue;
	    a = compare_id(line, line2);
	    if (a > 0) {
		if (state == 1) {
		    printf("%s", line2);
		    if (fgets(line2, sizeof line2, fpread))
                        printf("%s", line2);
		}
		continue;
	    }
	    break;
	}
	if (a < 0) {
	    fprintf(stderr, "%s%s\n", line, line2);
	    fatal("reads file does not contain one of the match headerline or oreder not correct\n");
	}
 	if (fgets(line2, sizeof line2, fpread)) {
	    printf("%s%s", line, line2);	
	} else {
	    fatal("read file contain a header without read sequence\n");
	}
    }
    fclose(fm);
    if (state) {
	while (fgets(line2, sizeof line2, fpread)) {
	    char *q = strchr(line2, ',');
	    if (q) {
		*q = '\n';
		*(q+1) = 0;
	    }
	    printf("%s", line2);
	}
    }
    fclose(fpread);

}
