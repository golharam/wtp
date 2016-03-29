#include "util.h"
#include "zutil.h"

main(int argc, char *argv[])
{
    printf("#Schema data base\n");
    char *list = argv[1];

    FILE *fp = ckopen(list, "r");
    char line[10000];
    while (fgets(line, sizeof line, fp)) {
	int a = -1, b = -1;
	char s[100];
	s[0] = 0; 
	sscanf(line, "schema_%d_%d_%s", &a, &b, s);
	if (a < 0 || b < 0)  fatal("file name not in right format");
	printf("$%d %d ", a, b);
	if (s[0] && strncmp(s, "adj", 3) == 0) {
	    printf("1\n");
	} else {
	    printf("0\n");
	}
	fflush(stdout);
	char com[10000];
	sprintf(com, "cat %s", line); 
	int i = system(com);
	if (i != 0) {
	    fprintf(stderr, "command %s failed\n", com);
	    exit(1);
	}
    }
} 
