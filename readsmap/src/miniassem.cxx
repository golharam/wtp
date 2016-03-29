#include "miniassem.h"

void graph::init(int mn, int rl)
{
    if (rlist) delete [] rlist;
    max_num_reads = mn+mn;
    rlist = new seqread[max_num_reads];
    num_reads = 0;
    if (p_queue) delete p_queue;
    p_queue = new heap(mn); 
    read_length = rl;
}

graph::~graph()
{
    if (rlist) delete [] rlist;
    if (p_queue) delete p_queue;
}

void graph::inputReads(char *file)
{
    FILE *fp = ckopen(file, "r");
    char line[10000], line1[1000];
    int i = 0;
    num_reads = 0;
    while (fgets(line, sizeof line, fp)) {
	if (line[0] == '#') continue;
	if (line[0] != '>') fatal("format of input wrong\n");
 	if (!fgets(line1, sizeof line1, fp)) break; 
	if (num_reads > max_num_reads) fatal("too many reads\n");
	int dir = 1; 
 	if (strchr(line, '@')) { dir = -1;}
	rlist[num_reads].set(num_reads, line1+2, line+1, dir);
	num_reads++;
    }
}
// edge from's 3' end overlap to's 5' end, 
void graph::add_overlap(int overlap, seqread *from, seqread *to)
{
    edge *e = new edge(overlap, from, to->overlap);
    to->overlap = e;
    edge *f = new edge(overlap, to, from->overlap_b);
    from->overlap_b = f; 
}
//find a shortest path to all the reads from a starting read, then find the longest path. 
int graph::allshortest_path(seqread *start, int pos1,  int overlapstart,  char *res, int maxlen, int forb)
//forb = 0 means going forward
{
    ithrun+=2;
    seqread *n, *last = NULL;
    p_queue->renew();
    edge *ep;

    for (n = start; n; ) {
	last = n;
	ep = (forb == 0)? n->overlap_b : n->overlap;
        for (; ep; ep = ep->next) {
            seqread *t = ep->opp;
            //fprintf(stderr, "%d ", t->number);
            if (!t->visited(ithrun-1)) {
                t->setvisit(ithrun-1);
                t->weight = n->weight - ep->over_len + read_length;
                t->prev = n;
                //fprintf(stderr, "inserted\n");
                p_queue->insert(t);
            } else if (t->visited(ithrun)) {
                //fprintf(stderr, "already visited\n");
                continue;
            } else {
                //fprintf(stderr,"visited not processed\n");
                int w = n->weight - ep->over_len+read_length;
                if (t->weight > w) {
                    p_queue->decrease_key(t, w);
                    t->prev = n;
                }
            }
        }
        n = p_queue->delete_min();
        if (n) n->setvisit(ithrun);
    }
    // from last to get the longest path.
	char *s = res;
	n = last;
	char *tempstring = n->colorseq;
	char tmp[1000];
	if (forb == 0){ reverse(n->colorseq, read_length, tmp); tempstring = tmp;}
            strncpy(s, tempstring, read_length);
            s += read_length;
            while (n != start) {
                int a = n->weight-n->prev->weight;
                n = n->prev;
                if (s-res+a >maxlen) return -1;
		if (forb == 1) {
		    tempstring = n->colorseq;
		} else {
		    reverse(n->colorseq, read_length, tmp); tempstring = tmp
;
		}
                strncpy(s,tempstring+read_length-a, a);
                s += a;
            }
            s -= read_length-1-overlapstart-1;
            *s = 0;
	    if (forb==0) reverse_seq(res, s-res);
            return s - res;
}

int graph::shortest_path(seqread *start, seqread *end, int pos1, int pos2, int overlapstart, int overlapend, char *res, int maxlen)
{
    ithrun+=2; 
    seqread *n;
    p_queue->renew();
    edge *ep;
    for (n = end; n; ) {
	for (ep = n->overlap; ep; ep = ep->next) {
	    seqread *t = ep->opp;
	    //fprintf(stderr, "%d ", t->number);
	    if (!t->visited(ithrun-1)) {
		t->setvisit(ithrun-1);
		t->weight = n->weight - ep->over_len + read_length;  
		t->prev = n;
		//fprintf(stderr, "inserted\n");
		p_queue->insert(t);
	    } else if (t->visited(ithrun)) {
		//fprintf(stderr, "already visited\n");
		continue;
	    } else {
		//fprintf(stderr,"visited not processed\n"); 
		int w = n->weight - ep->over_len+read_length;
		if (t->weight > w) {
		    p_queue->decrease_key(t, w);
		    t->prev = n;
		}
	    }
	}
	n = p_queue->delete_min();
	if (n) n->setvisit(ithrun);
	//if (n) fprintf(stderr, "delete min %d\n", n->number);
	if (n == start) {
	    // report chain
	    char *s = res;
	    strncpy(s, n->colorseq+overlapstart, read_length-overlapstart);
	    s += read_length-overlapstart;
	    while (n != end) {
		int a = n->weight-n->prev->weight;
		n = n->prev;
		if (s-res+a >maxlen) return -1;
		strncpy(s, n->colorseq+read_length-a, a);
		s += a;
	    }
	    s -= read_length-1-overlapend+pos2-pos1-1;
	    *s = 0;
	    return s - res;
	}
    }
    return 0; 
}

static int hashof(const char *seq, int b, int len)
{
    int i;
    int hash = 0;
    for (i = b; i < b+len; i++) {
	hash = (hash << 2) |( seq[i]-'0');
    }
    return hash;
}

void graph_assem::check_overlap(int mini_over)
{
    int real_over = mini_over;
    if (mini_over > 13) mini_over = 13;  
    int hsize = (1 << (2*mini_over));
    int mask = hsize-1; 
    if (htb) delete htb;
    htb = new hashtable(hsize); 
    int i;
    for (i = 0; i < num_reads; i++) {
	seqread *n = rlist+i;
	int h = hashof(n->colorseq, 0, mini_over);	
	//fprintf(stderr, "%d %d %s\n", i, h, n->colorseq);
	hashnode *np = new hashnode(h, n);
	CHECK_MEM(np);
	htb->insert(np, h);
    }
    for (i = 0; i < num_reads; i++) {
        seqread *n = rlist+i;
	int j;
	int h = hashof(n->colorseq, 0, mini_over-1);
	hashnode *hn;
	for (j = mini_over-1; j < read_length; j++) {
	    h = (h << 2) | (n->colorseq[j]-'0');
	    h = h & mask;
	    for (hn = htb->find(h); hn; hn = hn->next) {
		int k;
		seqread *m = hn->seqr();
		for (k = j+1; k < read_length; k++) {
		    if (n->colorseq[k] != m->colorseq[k-j-1+mini_over]) break;
		}
		if (k >= read_length && read_length-j+mini_over-1 >= real_over) {
		    //fprintf(stderr, "add overlap\n");
		    add_overlap(read_length-j+mini_over-1, n, m);	    
		}
	    }
	}
    }
}

static void addto(int *Pos, int *Opos, int *Read, int &num, int maxnum, int pos, int rp, int r)
{
    int a;
    //printf("%d %d %d\n", pos, rp, r);
    if (num >= maxnum) fatal("too many hits \n");
    for (a = num-1; a >=0; a--) {
	if (Pos[a]> pos) {
	    Pos[a+1]=Pos[a]; 
	    Read[a+1] = Read[a];
	    Opos[a+1] =Opos[a];
	} else {
	    break;
	}
    }
    Pos[a+1] = pos;
    Read[a+1] = r;
    Opos[a+1] = rp;
    num++;
}

// this is run after check_overlap so htb is already built for all the reads
void graph_assem::assembly(const char *refseq, int start, int end, int mini_over)
{
    int real_over = mini_over;
    if (mini_over > 13) mini_over = 13;
    int hsize = (1 << (2*mini_over));
    int mask = hsize-1;
    hashtable *backhalf = new hashtable(hsize);
    int i;
    for (i = 0; i < num_reads; i++) {
        seqread *n = rlist+i;
        int h = hashof(n->colorseq, read_length-mini_over, mini_over);
        hashnode *np = new hashnode(h, n);
        backhalf->insert(np, h);
    }
    int numH = (end-start)*50;
    int forN, forP[numH], forR[numH], revN, revP[numH], revR[numH], forO[numH], revO[numH];
    int h = hashof(refseq, start, mini_over-1);
    for (i = start+mini_over-1; i <= end; i++) {
	h = (h << 2) | (refseq[i]-'0');
	h = (h & mask);
	hashnode *hn;
	for (hn = htb->find(h); hn; hn = hn->next) {
	    seqread *n = hn->seqr();
	    int pos, repos;
	    for (pos = i+1, repos = mini_over; repos < read_length; repos++, pos++) {
		if (refseq[pos] != n->colorseq[repos]) break;
	    }
	    //printf("first ");
	    if (repos >= real_over)
	    	addto(forP, forO, forR, forN, numH, pos, repos, n->number); 
	}
	for (hn = backhalf->find(h);  hn; hn = hn->next) {
            seqread *n = hn->seqr();
            int pos, repos;
            for (pos = i-mini_over, repos = read_length-mini_over-1; repos >= 0 && i >=start; repos--, pos--) {
                if (refseq[pos] != n->colorseq[repos]) break;
            }
	    //printf("back ");
	    if (read_length-repos-1 >= real_over)
                addto(revP, revO, revR, revN, numH, pos, repos, n->number);
        }
    }
    // merge two list of overlapping to find start and end read to do a shortestpath; 
    //fprintf(stderr, "%d %d\n", forN, revN);
    int j;
    i = j = 0;
    int pos =  forP[i], revp = revP[j]; 
    char res[10000];
    for (i = 0; i < forN; i++)  {
	for (j = 0; j < revN; j++) {
	    int diff = revP[j]-forP[i];
	    if (diff > -5 && diff < 5) {
		//fprintf(stderr, "%d %d %d %d\n", i, j, forP[i], revP[j]);
		int x = shortest_path(rlist+forR[i], rlist+revR[j], forP[i]-1,  revP[j]+1, forO[i], revO[j], res, 3000);
		if (x > 0) {
		    printf("len =%d start after %d \n%s\n", x,forP[i], res);
		    return;
		}
	    } 
	}
    }
}


