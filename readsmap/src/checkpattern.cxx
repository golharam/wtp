#include "util.h"
#include "zutil.h"
#define maxsize 10000

int check(int *num, int n, int len, char **pat, int m, int len1, int adj)
{
    int i, j;
    for (i = 0; i < m; i++) {
	for (j = 0; j < n; j++) {
	    if (num[j] == 0) continue;
	    if (num[j] > len1) continue;
	    if (pat[i][num[j]-1] == '1') break;
	    if (adj == 0 || num[j] > len1-1) continue;
	    if (pat[i][num[j]] == '1') break;
	}
	if (j >= n) return 1;
    }
    //for (i = 0; i < n;i++) printf("%d\n", num[i]);
    //exit(1);
    return 0;
}

int numberofones(char *l)
{
    char *lp = l;
    int co = 0;
    while (*lp) {
	if (*lp == '1') co++;
	lp++;
   }
    return co;
} 

main(int argc, char *argv[])
{
    if (argc == 1) {
	fatal("checkpattern +/-nummis length schema_file [unit]\n");
    }
    int n = atoi(argv[1]);
    int adj = 0;
    int len = atoi(argv[2]);
    FILE *fp = ckopen(argv[3],"r");
    if (n < 0) { n = -n; adj = 1;}
    char *unit = NULL, *zero;
    char *unito[len], *zeroo[len];
    if (argc > 4) {
	unit = argv[4];
	int i;
        zero = strsave(unit);
        for (i = 0; i < strlen(unit); i++)
            zero[i] = '0';

	if (strchr(unit,',')) {
	    char *q = unit;
	    for (i = 0; i < len; i++) {
		unito[i] = q;
		zeroo[i] = zero+(q-unit);
		q = strchr(q,',');
		if (q) {*q = 0; *(zero+(q-unit)) = 0; q++; } else break;
	    }
	    if (i < len-1) fatal("unit not long enough\n");
	} else {
	    for (i = 0; i < len; i++) {
		unito[i] = unit;
		zeroo[i] = zero;
	    }
	}
    }
    char *pat[maxsize];
    char line[10000];
    int len1; 
    int m = 0,max1 = 0, min1=10000 ;
    while (fgets(line, sizeof line, fp)) {
	if (m >= maxsize) fatal("too many lines in schema file\n");
	if (line[0] == '#') continue;
	pat[m] = strsave(line);
	len1 = strlen(line)-1;
	if (len1 != len-1) fatal("line of different size\n");
	int x = numberofones(line);
	if (x > max1) max1= x;
	if (x < min1) min1 = x; 
	m++;
    } 
    int num[n];
    int i, j ;
    long long found =0, notfound=0;
    for (i = 0; i < n; i++) num[i] = i;
    do {
	if (check(num, n, len, pat, m, len1, adj)) found++; 
	else notfound++;
	if (unit) {
	    int k;
	    for (k = 0, j = 0; j < len; j++) {
		if (k < n && j == num[k]) {
		    printf("%s", unito[j]);
		    k++;
		} else {
		    printf("%s", zeroo[j]);
		}
	    }
	    printf("\n");
	}
	
	for (j = n-1; j>=0; j--) {
	    if (num[j] < len-1- (n-1-j)) {
		num[j]++;
		int k;
	    	for (k = j+1; k < n; k++) 
		    num[k] = num[k-1]+1;
		break;
	    }
	}
    } while (j >= 0);
    printf("found =%d \nnotfound=%d\n%d%% found\n", found, notfound, (int) (100*found/(found+notfound)));
    printf("total lines=%d\nmax one=%d\tmin ones=%d\n", m,max1, min1);
}
