#include "util.h"
#include "zutil.h"
#define MAXLEN 3000
#define MAXNLINE 1000000
static char **schemalns;

class SS {
  public:
    SS() {start = NULL;}  
    double score() {
    	int i; 
	double j = 1.0;
	for (i = 27; i > maxOne+minOne; i-- ) {
	     j *= 1.3;
	}
	return j*num;
    }
    void output() {
	int i;
	for (i = 0; i < num; i++) {
	    printf("%s\n", start[i]);
	}
    }
    void output(int m) {
        int i;
        for (i = 0; i < num; i++) {
            printf("%s,%d\n", start[i], m);
        }
    }
    char **start;
    int num;
    int maxOne;
    int minOne;
};

static SS *schema[MAXLEN];

static int num_of_ones(const char *a)
{
    int x = 0;
    while (*a) {
	if (*a == '1') x++;
	a++;
    }
    return x;
}


static int ReadSchemas(FILE *tmpfile, int len, int mis, int use_adj)
{
    char line[10000];
    SS *s = NULL;
    int LN = 0;
    while (fgets(line, sizeof line, tmpfile)) {
	if (line[0] == '#') continue;
        if (line[0] == '$') {
            int a, b, c;
            a = b = c = -1 ;
            sscanf(line+1, "%d %d %d", &a, &b, &c);
	    s = NULL;
            if (c < 0) continue;
            if (a <= len && b <= mis && (b<=1 || c == use_adj)) {
		//find a useful schema
		s = schema[a]+b;
		s->start = schemalns+LN;
		s->num = 0;
		s->maxOne = 0;
		s->minOne = a; 
	    }
        } else if (s) {
	    char *t = strchr(line, '\n');
  	    if (t) *t = 0;		
	    schemalns[LN] = strsave(line);
	    LN++;
	    s->num += 1; 
	    int num_ones = num_of_ones(line);
	    if (num_ones > s->maxOne) s->maxOne = num_ones;
	    if (num_ones < s->minOne) s->minOne = num_ones;   
	}
    }
}

static void fill(char *s, int n, char f)
{ 
    int i;
    for (i = 0; i < n; i++) s[i] = f;
}

static void fillzero(char *s, int n) {
	fill(s, n, '0');
}

static void fillone(char *s, int n) {
        fill(s, n, '1');
}


static void extend(int len, int mis, int slen)
{
    SS *s = schema[slen]+mis;
    char **sc = new char*[s->num];
    int i; 
    for (i = 0; i < s->num; i++) { 
	sc[i]= new char[len+1];
	fillzero(sc[i], len-slen);
	strcpy(sc[i]+len-slen, s->start[i]);
    }
    SS *t = schema[len]+mis;
    t->num = s->num;
    t->minOne = s->minOne;
    t->maxOne = s->maxOne;
    t->start = sc;
}

static void merge(int len, int mis, int flen, int fmis, int rlen, int rmis)
{
    SS *f = schema[flen]+fmis;
    SS *r = schema[rlen]+rmis;
    SS *t = schema[len]+mis;
    int i;
    t->num = f->num + r->num;
    char **sc = new char*[t->num];
    t->start = sc;
    for (i = 0; i < f->num; i++) {
	sc[i] = new char[len];
	sc[i][len-1] = 0;
	strcpy(sc[i], f->start[i]);
	fillzero(sc[i]+flen-1, len-flen);
    }
    for (i = 0; i < r->num; i++) {
	char *line = sc[i+f->num] = new char[len]; 
	line[len-1] = 0;
	fillzero(line, len-rlen);
	strcpy(line+len-rlen, r->start[i]);
    }
    t->minOne = (f->minOne < r->minOne)? f->minOne : r->minOne;
    t->maxOne = (f->maxOne > r->maxOne)? f->maxOne : r->maxOne;
}


static int build_new_simple(int len, int mis, int adj) 
{
    if (len < 14) return 0;
    if (schema[len][mis].start) return 1;
    int i, j;
    double bestscore = 1e+50;
    int flen = -1, fmis = -1, rlen = -1, rmis = -1;
    for (i = len-1; i >= 15; i--) {
	SS *s = schema[i]+mis;
	double ss;
	if (s->start) {
	    ss= s->score();
	    if (ss > bestscore) {
		bestscore = ss;
		flen = i;
		fmis = mis;
	    }
	}
    }
    for (i = len-14; i >= len/2-1; i--) {
	for (j = mis-1; j >= 0; j--) {
	    SS *s =schema[i]+j;
	    if (s->start == NULL) continue;
	    int a = len-i+(1-adj);
	    int b = mis-j-1;
	    SS *t = schema[a]+b;
	    if (t->start == NULL) continue;
	    double ss = s->score()+t->score();
	    if (ss < bestscore) {
		bestscore = ss;
		flen = i; fmis = j;
		rlen = a; rmis = b;
	    }
	}
    }
    if (flen < 0) return 0;
    if (rlen <0) {
	extend(len, mis, flen);
    } else {
	merge(len, mis, flen, fmis, rlen, rmis);
    }
    return 1;
}

static int build_new_rec(int len, int mis, int adj) 
{
    if (len < 14) return 0;
    if (schema[len][mis].start) return 1;
    if (mis <= 2) {
	return build_new_simple(len, mis, adj);
    }
    int m1 = (mis-1)/2, m2 = mis-1-m1;
    int flen = (len+1-adj)*m1/(mis-1), rlen = len-flen+1-adj;
    if (build_new_rec(flen, m1, adj) && build_new_rec(rlen, m2, adj)) {
	merge(len, mis, flen, m1, rlen, m2);
	return 1;
    }
    return 0;
}
 
	  

main(int argc, char *argv[])
{
    if (argc < 5) {
	fatal("newSchema (+/-)length Mis 0/1[Adj or not] DB\n To get k-best schemas use negative number to specify length\n");
    } 	
    int kbest = 0;
    int length = atoi(argv[1]);
    if (length < 0) { kbest = 1; length = -length;}
    int Mis = atoi(argv[2]);
    int adj = atoi(argv[3]);
    if (adj != 0) adj = 1;
    if (length > MAXLEN) fatal("Read length too long\n");
    int x = (int) (length*0.2);
    if (x <= Mis) fatal("Miss level too high\n"); 
    int i;
    schema[0] = new SS[(length+1)*(Mis+1)]; 
    for (i = 1; i <= length; i++) {
	schema[i] = schema[i-1]+Mis+1;
    }
    schemalns = new char*[MAXNLINE];
    FILE *fp = ckopen(argv[4], "r");
    char line[10000];
    if (fgets(line, sizeof(line), fp)== NULL) exit(1);
    if (strncmp(line, "#Schema data base", 17) != 0) fatal("Input file format not right\n"); 
    ReadSchemas(fp, length, Mis, adj);
    if (kbest) {
	int m;
	for (m = 0; m <= Mis; m++) {
	    if (length < 100 &&  build_new_simple(length, m, adj)) {
        	schema[length][m].output(m);
    	    } else if (build_new_rec(length, m, adj)) {
        	schema[length][m].output(m);
	    } else {
		fprintf(stderr, "Not complete in building k-best at %d_%d\n", length, m);
	    }
	}
	return 0;
    }
    if (length < 100 &&  build_new_simple(length, Mis, adj)) {
	schema[length][Mis].output();
    } else if (build_new_rec(length, Mis, adj)) {
	schema[length][Mis].output();
    }
}
