#include "util.h"
#include "zutil.h"

main(int argc, char *argv[])
{
    if (argc==1) fatal("schema2adj schemaFile\n");
    FILE *fp = ckopen(argv[1], "r");
    char line[10000];
    while (fgets(line, sizeof line, fp)) {
	if (line[0] == '#') {printf("%s", line); continue;}
	char *q = line;
	int state = 1;
	while (*q) {
	    if (state == 0) {
		if (*q == '1') {
		    *q = '0';
		    state = 1;
		}
	    } else {
		if (*q == '0') state = 0;
	    }
	    q++;
	}
	printf("%s", line);
    }
}	  
