
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

main(int argc, char *argv[])
{
    FILE *fp = fopen(argv[1], "r");
    int cut = atoi(argv[2]);
    int limit_item = atoi(argv[3]);
    int state = 0, item;
    if (argc > 4) state = atoi(argv[4]);
    if (argc < 3) {fprintf(stderr, "%s indel_result num_uniq max_output [state]\n", argv[0]); exit(1);} 
    char line[100000], line1[10000];
    char outline[1000000];
    long long F[1000], R[1000];
    long long p = ((long long) -100);
    int  n = 0, m = 0, large=0, small=10000;
    outline[0] = 0;
    while (fgets(line, sizeof line, fp)) {
	int h, f; 
	long long i;
	sscanf(line, "%lld %d", &i, &h);
	if (i > p+2) {
	    if (n >= cut /*&& large-small <5*/) {
		if (state) {
		    int ch = (int) (p >> 32);
		    long long o = (((long long) 1) << 32)-1;
		    p = p & o;
		    printf("%d:%lld-%lld\t%d(%d)\t%d-%d%s\n", ch, p, p+2, n, m, small, large, outline);
		} else printf("%lld-%lld\t%d(%d)\t%d-%d%s\n", p, p+2, n,m, small,large,outline);
	    }
	    p = i;
	    large = small = h;
	    n = m= 0;
	    outline[0] = 0; 
	    item = 0;
	}
	if (item < limit_item) { 
	    char *p = strchr(line, '\n');
	    if (p) *p = 0;
	    p = strchr(line,'>');
            strcat(outline, "\t");
	    strcat(outline, p);
	    item++;
	}
	    int j;
	    if (large < h) large = h;
	    if (small > h) small = h;
	    long long s1, s2;
	    char *p = strchr(line, ',');
	    if (state) {
		int a, b;
		sscanf(p+1, "%d_%d", &a, &b);
		//printf("%d %d\n", a, b);
		s1 = (((long long) a)<<32)+b;
	    } else s1 = atoll(p+1);
	    p = strchr(p, '|');
            if (state) {
                int a, b;
                sscanf(p+1, "%d_%d", &a, &b);
                s2 = (((long long) a)<<32)+b;
            } else s2 = atoll(p+1);

	    for (j = 0; j < n; j++) {
		if (F[j] == s1 && R[j] == s2) break;
	    }
	    m++;
	    if ( j == n) {
		F[n] = s1; R[n] = s2;
		n++;
	    }
    }
}
