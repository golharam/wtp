#include "util.h"
#include "zutil.h"
#include <math.h>
#define low_ratio 350.0 //single*0.1 
#define single 1000
int max_mis = 3;
int err = 9; // equivalent to 0.1 error rat.

static void low_coverage(int *cover, int seqlen, int relen, int num_read)
{
    int cur_low_start = 0; 
    int cur_hi_start = 0;
    int cutoff_cov = (int) (low_ratio * (double) num_read* (double)relen/ (double) seqlen);
    fprintf(stderr, "%d\n", cutoff_cov);
    int j, k, i;
    for (i = 0; i < seqlen; i++) {
	if (cover[i] <= cutoff_cov) {
	    if (cur_hi_start < i) { //end of a high run
		if (i - cur_hi_start > 3) {
		    int len = cur_hi_start-cur_low_start;
		    if (len <= relen /*&& len > relen/6*/) {
			printf("%d %d\n", cur_low_start, cur_hi_start);
		    }
		    cur_low_start = i;
		}
	    }
	    cur_hi_start = i+1;
	} else { //hi coverage
	    if (cur_hi_start >= i) { // end of a low run
		cur_hi_start = i;
	    }
	}
    }
    if (i - cur_hi_start <= 3) {
	cur_hi_start = i;
    }
    int len = cur_hi_start-cur_low_start;
    if (len <= relen && len > relen/4) {
        printf("%d %d\n", cur_low_start, cur_hi_start-1);
    }
}

static void coverage(int *cover, int seqlen, int relen, char *mapfile, int & num_read)
{
    int i;
    int power[max_mis+1];
    int c = 1;
    power[max_mis] = 1;
    for (i = 1; i<=max_mis; i++) {
	c *= err;
	power[max_mis-i] = c;
    } 
    FILE *fp = ckopen (mapfile, "r");
    char line[200000];
    while (fgets(line, sizeof line, fp)) {
	char *p = line;
	if (strchr(line, ',')) num_read++;
	int n = 0;
	while (p = strchr(p, ',')) {
	    int mis;
	    p++;
	    sscanf(p, "%*d,%d", &mis);
	    n += power[mis];  
	}
	p = line;
	while (p = strchr(p, ',')) {
	    int pos, mis;
	    p++;
	    sscanf(p, "%d,%d", &pos, &mis);
	    int w = single *  power[mis]/ n; 
	    if (pos < 0) pos = -pos-relen;
	    int i;
	    for (i = 0; i < relen; i++) {
		cover[i+pos] += w;
	    }
	}
	if (fgets(line, sizeof line, fp)== NULL) break;
    }
    fclose(fp);
}

main(int argc, char *argv[])
{
    if (argc < 2) {
	fatal("lowcov readlen seqlen mapfile1 [mapfile2]\n");
    }
    char *mapfile1 = argv[3];
    char *mapfile2 = NULL;
    if (argc >4) mapfile2 = argv[4];
    int relen = atoi(argv[1]);
    int seqlen = atoi(argv[2]);
    int *cover = new int[seqlen];
    CHECK_MEM(cover); 
    memset(cover, 0, seqlen*sizeof(int));
    int num_read = 0;
    coverage(cover, seqlen, relen, mapfile1, num_read);
    if (mapfile2) coverage(cover, seqlen, relen, mapfile2, num_read);
    low_coverage(cover, seqlen, relen, num_read);
}
