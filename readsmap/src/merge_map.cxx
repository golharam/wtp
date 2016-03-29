#include "util.h"
#include "zutil.h"
#include "zcompress.h"

#define MAX_FILE 5000
#define maxline 500000 

long long  numMas = -1;
long long begin =0, end = -1; 
int addchrom = 0;
int incremental = 0;

/* 
Notes:
8/18/08
    incremental merging:
	0. default, regular merging, every read in input read file will be
	 output. If defline have hit list, it will be removed and replaced by
	 new hits found.
	1. Incremental mapping, if a read has a hit in the defline, it will not
	  be mapped this run, and will not be output.
	2. Incremental 2, if a read has a hit in defline, it will not be mapped 
	  this run, but the hits in the defline line will be output along with
	  the read.   
*/


class triple {
  public:
    void set(int x, int y, int z) {
	a = x; b= y; c= z;
    }
    void copy(triple *x) { set(x->a, x->b, x->c);}
    int a, b, c;
};

int compare_id(triple *t1, triple *t2)
{
    int a = t1->a, b= t1->b, c= t1->c;
    int x = t2->a, y = t2->b, z=t2->c;
    if (a < x) return -1;
    if (a > x) return 1;
    if (b < y) return -1;
    if (b > y) return 1;
    if (c > z) return 1;
    if (c < z) return -1;
    return 0;
}


class node {
  public:
    node(int sn);
    void readin(char *line);
    triple *id;
    char match[maxline];
    int seqno;
    node *next;
};

node::node(int x)
{
    seqno=x;
    id = new triple();
    next = NULL;
}

void node::readin(char *line)
{
    int ll = strcspn(line, "\r\n");
    line[ll] = 0;
    if (strlen(line) > maxline) fatal("line too long\n");  
    strcpy(match, line);
    int a, b, c;
    sscanf(line, ">%d_%d_%d", &a, &b, &c);
    id->set(a, b, c);
}

int compare_node(node *a, node *b) {
    return compare_id(a->id, b->id);
}


class heap {
  public:
    heap() {num = 0; head = NULL;}
    void insert(node *);
    node *delete_min();
  protected:
    int num;
    node *head;
};

void heap::insert(node *a)
{
    if (head == NULL) {
	a->next = NULL; num=1; head =a; return;
    }
    node *c = head, *prev = NULL;
    while (c && compare_node(c, a) < 0) {
	prev = c; c = c->next; 
    }
    if (prev) { prev->next = a;} else {head = a;}
    a->next = c;
    num++;
}

node *heap::delete_min()
{
    node *t = head;
    if (head) head = head->next;
    return t;
}

int getline(char *a, int s, FILE *fp)
{
    while (fgets(a, s, fp)) {
	if (a[0] == '>') return 1;
    }
    return 0;
}

static void getrest(FILE *fp)
{
    if (fp == NULL) return;
    char line[10000];
    while (fgets(line, sizeof line, fp)) {
	if (end < 0)  printf("%s", line);
	else {
	    if (line[0] == '>') numMas++;
	    if (numMas > end) return;
	    printf("%s", line);
	}
    }
}

static void getother(FILE *fp, triple *id)
{
    if (fp == NULL) return;
    triple aa;
    triple *idnew = &aa;
    char line[10000];
    int do_output = 1;
    while (fgets(line, sizeof line, fp)) {
	if (line[0] != '>') {
	    if (numMas < begin) continue;
	    if (line[0] != '#' && do_output) printf("%s", line);
	    continue;
	}
	numMas++;
	if (numMas < begin) continue;
	if (end >= 0 && numMas > end){
	    fprintf(stderr, "Bead id %d_%d_%d\n", id->a, id->b, id->c); 
	    fatal("Did not find ID from the master read file\n");
	}
 	int a, b, c;
	sscanf(line, ">%d_%d_%d", &a, &b, &c);
	idnew->set(a, b, c);
	int x = compare_id(idnew, id);
	if (x < 0) {
	    char *q = strchr(line,',');
	    do_output = 1;
	    if (q == NULL) {
		 printf("%s", line);
	    } else {
		if (incremental == 2) printf("%s", line);
		else if (incremental == 0) {
		    *q = 0;
		    printf("%s\n", line); 
		} else do_output = 0; 
	    }
            continue;
	}
	if (x == 0) return;
	fprintf(stderr, "Read A: Bead id %d_%d_%d\n", id->a, id->b, id->c);
    	fprintf(stderr, "Current Read B: Bead id %d_%d_%d\n", idnew->a, idnew->b, idnew->c);
	fatal("Can not find read A. Either it is not in read file, or it is after B, which means the read file is not sorted in lexicographical order\n");

    }
}

void output(char *p, node *m)
{
    char *q;
    if (addchrom == 0 || (q = strchr(p, ','))== NULL) { printf("%s", p); return;}
    *q = 0; 
    if (q != p) printf("%s", p);
    do {
	p = q+1;
	printf(",%d_", m->seqno);
	q = strchr(p,',');
 	if (q) *q = 0;
	printf("%s", p);	
    } while (q);
}

int main(int argc, char *argv[])
{
    //Individual chromosome files
    FILE *fp[MAX_FILE];
    memset(fp, 0, sizeof(FILE *)*MAX_FILE);

    //cmap and match file
    FILE *flist = NULL, *masterFile = NULL; 
    genFile gF_masterFile; gF_masterFile.setUsePipe();
    
    char line[maxline+10];
    heap *hp = new heap();
    int i = 0, chromList = 0, maxlist = 0;
    char *inputdir = ".";
    char *genericMatchfile = "match"; 

    for (i=1; i<argc; ++i)
    {
        if (argv[i][1] == '=')         // X=value
        {
            if (argv[i][2] == 0)
                fprintf(stderr,"Missing value for %s\n",argv[i]);
            else if (argv[i][0] == 'S') chromList = atoi(argv[i]+2);
	    else if (argv[i][0] == 'I') inputdir = argv[i]+2;
	    else if (argv[i][0] == 'G') genericMatchfile = argv[i]+2;
	    else if (argv[i][0] =='B') begin = atoll(argv[i]+2);
	    else if (argv[i][0] =='E') end  = atoll(argv[i]+2);
	    else if (argv[i][0] =='A') addchrom = atoi(argv[i]+2);
	    else if (argv[i][0] == 'i') incremental = atoi(argv[i]+2);
		else { fprintf(stderr,"unknown option %c\n", argv[i][0]);}
	} else {
	    if (flist == NULL) flist = ckopen(argv[i], "r");
	    else if (masterFile == NULL) {
	      gF_masterFile.open(argv[i], "r");
	      masterFile = gF_masterFile.getFILE(); //ckopen(argv[i], "r");
	    }
	}
    }
    if (flist == NULL) {
	fatal("merge_map(1.0 April 24,2009) flist [masterFile] [S=state] [I=inputdir(./)] [G=genericMatchfile (match)][B=begin][E=end][A=addchrom][i=incremental]\n S=0 filelist; S=1 chromosome list; i=0 no inc, 1 inc, 2 inc merge\n");
    }
    if (chromList== 0) {
      i = 1;
      while (fgets(line, sizeof line, flist)) {
	char fname[10000];
	if (i > MAX_FILE) fatal("too many files to merge\n");
	sscanf(line, "%s", fname);
	fp[i] = ckopen(fname, "r");
	if (getline(line, sizeof line, fp[i])) {
	    node *n = new node(i);  
	    n->readin(line);
	    hp->insert(n);
	}
	i++;
      }
      maxlist = i-1;	
    } else {
	i = 0;
	while (fgets(line, sizeof line, flist)) {
	    if (line[0] == '#') continue;
	    char a[1000], fname[10000];
	    sscanf(line, "%d %s", &i, a);
	    if (i > MAX_FILE) fatal("too many files to merge\n");
	    if (i > maxlist) maxlist = i;
	    if (fp[i] != NULL) fatal("ERROR: in cmap file, $cmap_file, the chromosome ID $chr_id is written more than once.  Fix the cmap file.\n");
	    sprintf(fname, "%s/chr%s/%s", inputdir, a, genericMatchfile);
	    fp[i] = ckopen(fname,"r");
	    //printf("file %d %s\n", i, fname);
	    if (getline(line, sizeof line, fp[i])) {
            	node *n = new node(i);
            	n->readin(line);
            	hp->insert(n);
            }
	}
    }
    triple *cur = new triple();
    cur->set(0, 0, 0); 
    node *m;
    m = hp->delete_min();
    if (!m) exit(0);
    getother(masterFile, m->id);
    output(m->match, m);
    cur->copy(m->id);
    if (getline(line, sizeof line, fp[m->seqno])) {
        m->readin(line);
        hp->insert(m);
    }
    char colorseq[10000];
    while (m = hp->delete_min()) {
	if (compare_id(cur, m->id) == 0) {
	    char *p = strchr(m->match, ',');
	    if (p) output(p, m);
	} else {
	    printf("\n");
	    getother(masterFile, m->id);
	    output(m->match, m);
	    cur->copy(m->id);
	}
	if (getline(line, sizeof line, fp[m->seqno])) {
	    m->readin(line);
	    hp->insert(m);
	}
    }
    printf("\n");
    getrest(masterFile);
    for (i = maxlist; i >=0; i--) {
	if (fp[i]) fclose(fp[i]);
    } 
}
