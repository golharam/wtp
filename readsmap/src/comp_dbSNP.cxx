#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#define margin 10
int getnew_yoru(char *line, int s, FILE *fp, int &c, int &p)
{
    if (!fgets(line, s, fp)) return 0;
    int cc, pp; 
    sscanf(line, "%d:%d", &cc, &pp);
    c = cc; p =pp; 
    return 1;
}

int getnew_list(char *line, int s, FILE *fp, int &c, int &p)
{
    if (!fgets(line, s, fp)) return 0;
    int  pp;
    char ch[100];
    sscanf(line, "%*s %*s %*s %s %d", ch, &pp);

    int cc = -1;
    if (ch[0] == 'X') cc = 23;
    else if (ch[0] == 'Y') cc = 24;
    else if(ch[0] == 'M') cc= 25;
    else cc = atoi(ch); 
    c = cc;
    p =pp-1; //could make dbSNP pos 0 based here by subtracting 1
    return 1;
}


main(int argc, char *argv[])
{
    if (argc == 1) {
	fprintf(stderr, "comp result_file dbSNPfile\nRequires both files to be sorted by Chromosome then by position.\n");
	exit(1);
    }
    FILE *fp1 = fopen(argv[1], "r");
    FILE *fp2 = fopen(argv[2], "r");
    char line1[10000], line2[10000];
    int c1, c2, p1, p2;
    if (!getnew_yoru(line1,sizeof line1, fp1, c1, p1)) exit(0); 
    if (!getnew_list(line2,sizeof line2, fp2, c2, p2)) exit(0);
    do {
	if (c1 == c2 && p1 > p2-margin && p2 > p1-margin) {
	    printf("+%s-%s", line1, line2);
	    if (!getnew_yoru(line1,sizeof line1, fp1, c1, p1)) exit(0);
	    continue;
	}
	if (c1 < c2 || c1==c2 && p1 < p2+margin) {
            if (!getnew_yoru(line1,sizeof line1, fp1, c1, p1)) exit(0);
            continue;
	}
	if (!getnew_list(line2,sizeof line2, fp2, c2, p2)) exit(0);
    } while (1);
}
