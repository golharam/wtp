#include "util.h"
#include "zutil.h"
#include "mytemplate.h"

#define THREE_PRIME 0
#define FIVE_PRIME 1

// the class T must has a field index, and a field weight, both public
template <class T>
class _heap 
{
   public:
    _heap(int m) {
	max_size = m;
	list = new T*[m];
	num = 0;
    }
    ~_heap() {
	delete [] list;
    }
    void renew() { num = 0;}
    void  move_up(int n) {
	if (n == 0) { return;} 
	T *x = list[n];
	int p; 
	do {
	    p = (n-1)/2;
	    if (list[p]->weight <= x->weight) {
		set(x, n);
		return;
	    }
	    set(list[p], n);
	    n = p;
	    if (n == 0) { set(x, n); return;}
	} while (1);
    }
    T *delete_min() {
 	if (num == 0) return NULL;
	T *n = list[0];
	num--;
	if (num == 0) return n;
	int i = 0;
	T *s = list[num];
	int w = s->weight;
	while (i*2+1 < num) {
	    int left = i*2+1;
	    int right = left+1;
	    int j;
	    j = (list[left]->weight < list[right]->weight)? left : right; 
  	    if (list[j]->weight >= w) break;
  	    set(list[j], i);
	    i = j;
	}
	set(s, i);
	return n;
    }
    void insert(T *n) {
	if (num < max_size) {
	    set(n,num);
	    move_up(num);
	    num++;
	} else {
	    fatal("Exceeding size for heap\n");
	} 
    }
    void decrease_key(T *n, int w) {
	int i = n->index;
 	if (list[i] != n) fatal("Index wrong\n");
	if (w > n->weight) fatal("Not decrease key\n"); 
	n->weight = w;
	move_up(i);
    }
  protected:
    T **list;
    int num, max_size;
    void set(T *n, int i) {
	n->index = i;
	list[i] = n;
    }
};

class seqread;

class edge {
  public:
    edge(int ov, seqread *rs, edge *n) {
	over_len = ov; opp = rs; next = n;
    }
    int over_len;
    seqread *opp;
    edge *next;
};

class seqread 
{
  public:
    seqread() {
	number = 0;
	colorseq = id = NULL; 
	visit = 0;
	overlap = overlap_b = NULL;
	prev = NULL;
	index = weight = 0;
    }
    void set(int i, char *s, char *name, int d) {
	number = i; id = strsave(name); dir = d;
	if (d == -1) {
	    int i, len = strcspn(s, "\n\r");
	    colorseq = new char[len+1];
	    for (i = 0; i < len; i++) {
		colorseq[i] = s[len-i-1];
	    }
	} else {
	    colorseq = strsave(s);
	}
    } 
    int visited(int ithrun) { return (visit >= ithrun);}
    void setvisit(int ithrun) { visit = ithrun;}
    int number, dir;
    char *colorseq;
    char *id; 
    edge *overlap, *overlap_b;
    seqread *prev;
    int index, weight;
  protected:
    int visit;
};

typedef _heap<seqread> heap;

class graph
{
  public: 
    graph() {max_num_reads = num_reads = 0; rlist = NULL; p_queue = NULL; ithrun = 0;}
    ~graph();
    void init(int maxNum, int rl);
    void inputReads(char *file); 
    int readslist(seqread* &list) { list = rlist; return num_reads;}
    void add_overlap(int oversize, seqread *a,  seqread *b);
    int shortest_path(seqread *start, seqread *end, int pos1, int pos2, int overstart, int oveend, char *res, int maxlen);
    int allshortest_path(seqread *start, int pos1, int overstart, char *res, int maxlen, int for_or_back);
    // res shall have enough space to store the result string, the size of 
    // res is maxlen. return 0 is no path, -1 is overflow res. 
  protected:
    seqread *rlist;
    int num_reads;
    int max_num_reads;
    char dir;
    heap *p_queue;
    int read_length;
    int ithrun;
};

class hashnode {
  public:
    hashnode(int h, seqread *sr) {rd = sr, hash = h;}
    void set(int h, seqread *sr) {rd = sr; hash = h;}
    int same(hashnode *p){ return ( seqr()== p->seqr());}
    hashnode *next;
    int hashvalue(int x) {return hash;}
    seqread *seqr() {return rd;}
  protected:
    seqread *rd;
    int hash;
};

typedef _hash_table<hashnode> hashtable;

class graph_assem : public graph
{
  public:
    graph_assem(){;}
    void check_overlap(int mo);
    void assembly(const char *refseq, int start, int end, int mini);
  protected:
    hashtable *htb;
};
