#include "util.h"
#include "zutil.h"

static void addzero(int x)
{
    while (x > 0) {
	printf("0");
	x--;
    }
}

main(int argc, char *argv[])
{
    char *f1 = argv[1];
    char *f2 = argv[2];
    int len = atoi(argv[3])-1;
    FILE *fp1 = ckopen(f1, "r");
    FILE *fp2 = ckopen(f2, "r");
    char line[10000];

    while (fgets(line, sizeof line, fp1)) {
	if (line[0] =='#') continue;
	char *c = strchr(line, '\n');
	*c = 0;
	int len1 = strlen(line);
	printf("%s", line);
	addzero(len-len1);
	printf("\n");
    }
    while (fgets(line, sizeof line, fp2)) {
	if (line[0] =='#') continue;
        char *c = strchr(line, '\n');
        *c = 0;
        int len1 = strlen(line);
        addzero(len-len1);
        printf("%s", line);
        printf("\n");
    }
}

