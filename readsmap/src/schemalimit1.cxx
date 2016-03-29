#include "util.h"
#include "zutil.h"

main(int argc, char *argv[])
{
    if (argc==1) fatal("schemalimit1 schemaFile NumOnes\n");
    FILE *fp = ckopen(argv[1], "r");
    int numOnes = atoi(argv[2]);
    char line[10000];
    printf("#schemalimit1 %s %d\n", argv[1], numOnes);
    while (fgets(line, sizeof line, fp)) {
	if (line[0] == '#') {printf("%s", line); continue;}
	char *q = line+strlen(line);
	int count = 0;
	while (q >= line) {
		if (*q == '1') count++;
		if (count > numOnes) *q = '0';
		q--;
	}
	printf("%s", line);
    }
}
