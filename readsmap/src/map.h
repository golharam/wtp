#ifndef __map_h__
#define __map_h__ 1

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "zutil.h"
#include <math.h>
#include <pthread.h>

class make_hash {
 public:
    make_hash(){N=0;hashv=0; bitmap = 0;}
    virtual void addone(int a) {;}
    virtual void renew(){;}
    int hash() {if (N==0) return hashv; return -1;}
    void setpat(char *p, int w) { 
	pat = strsave(p); wd = w;
	int i;
	bitmask2 =(long long) 0;
	bitmask = (((long long )1) << wd)-1;
	//bitmask = (1 << wd)-1;
	//fprintf(stderr, "%lld %d\n", bitmask, wd);
	for (i = 0; i < wd; i++) {
	   bitmask2 <<= 1;
	   if (pat[i] == '1')  bitmask2++; 
	}
    }
    void addbase(char a) {
	bitmap <<= 1; 
	if (a == AMBIG)  {
	    bitmap++;
	}
	bitmap &= bitmask;
	if (bitmap & bitmask2) N=1;  else N=0;
    }
    void sethash(const char *seq, char *scode) {
	int i;
	hashv = 0; N=0;
	for (i = 0; i < wd; i++) {
	    if (pat[i] == '0') continue;
	    char a = scode[seq[i]];
	    if (a == AMBIG) {N=1; return;}
	    hashv <<=2;
	    hashv |= a;
	}
    }
 protected:
    int N, hashv;
    int wd;
    char *pat;
    long long bitmap, bitmask, bitmask2;
};

class makehash_3one_discont : public make_hash {
 public:
    makehash_3one_discont();
    virtual void renew() {
	N=0; bitmap = 0; 
	hashv = 0;
	h1 = h2 = h3=0;
    }
    void set_param(int zero1, int one1, int zero2, int one2, int zero3, int one3);
    virtual void addone(int a);
 protected:
    unsigned long long h1, h2, h3;
    int wd, shift2, shift1, shift3, shift2n;
    unsigned long long mask1, mask2, mask2n, mask3n, mask3;
};

class makehash_simple_discont : public make_hash {
 public:
    makehash_simple_discont();
    virtual void renew() {
	N=0; 
	hashv = 0;
	h1 = h2 = 0;
	bitmap = 0;
    }
    void set_param(int zero1, int one1, int zero2, int one2);
    virtual void addone(int a);
 protected:
    unsigned long long h1, h2;
    int wd; 
    int shift2, shift1;
    unsigned long long  mask1, mask2, mask;
};
class makehash_sdiscont : public make_hash {
 public:
    makehash_sdiscont();
    virtual void renew() {
        N=0;
        hashv = 0;
        h1 = h2 =h3= 0;
    }
    void set_param(int zero1, int one1, int zero2, int one2);
    virtual void addone(int a);
 protected:
    unsigned long long h1, h2, h3;
    int wd, shift2, shift1, shift3;
    unsigned long long  mask1, mask2, mask3;
};

class makehash_gen_discont : public make_hash {
 public:
    makehash_gen_discont();
    virtual void renew() {
	N=0; 
	hashv = 0;
	bitmap = 0;
	memset(h, 0, sizeof(long long)*n_ones);
    }
    void set_param(char *pat);
    virtual void addone(int a);
 protected:
    unsigned long long  h[100];
    int wd, n_ones, shift[100], shiftn[100];
    unsigned long long mask[100], maskn[100];
};

class makehash_cont : public make_hash {
 public:
    makehash_cont() {
	maskn = 0;
	wd = 0;
	N=0;
	hashv = 0;
    }
    void setwordsize (int w) {
	wd = w;
	maskn = (1 << (2*wd))-1;
    }
    virtual void renew() {
	N=0;
	hashv = 0;
    }
    virtual void addone(int a) {
	hashv <<= 2;
	hashv &= maskn;
	if (a == AMBIG) N= wd;
	else {
	    if (N >0) N--;
	    hashv |= a;
	}
    }
 protected:
    int wd;
    unsigned int maskn;
};

class seqc {
 public:
    seqc(const char *n, const char *seq, const char *seq1, int l, int d) {
	seqname = strsave(n);
	sequence = strsave(seq);
	if (seq1) seq2 = strsave(seq1); else seq2 = NULL;
	len = l;
	dir = d;
    }
    char *seqname;
    char *sequence;
    char *seq2;
    int len;
    int dir;
};

class  _hashtable {
 public:
    _hashtable(int t) {
	index = new unsigned int[t];
	maxcount = 0;
	memset(index, 0, sizeof(unsigned int)*t);
	plist = NULL;
	derived = 0;
	mt = NULL; num_m =  0;
    }
    _hashtable(_hashtable *copy) {
	index = copy->getIndex();
	plist = copy->getplist();
	derived = 1;
	maxcount = 0; 
        mt = NULL; num_m =  0;
    }
    void setMutex(pthread_mutex_t *m, int mm) {mt= m; num_m =mm-1;} 
    void setSeqlen(long long  l) {
	plist = new unsigned int[l];
	//memset(plist, 0, sizeof(unsigned int)*l);
    }
    void insert(unsigned int pos, int hash) {
	if (hash < 0) return;
	int y;
	if (mt) {
	    y = hash & num_m;
	    pthread_mutex_lock(&mt[y]);
	}
	if (plist == NULL) {
	    long long  x = (index[hash]++);
	    if (x > maxcount) maxcount = x; 
	    return;
	}
	unsigned int i = index[hash];
	plist[pos] = i;
	index[hash] = pos+1;
	if (mt) pthread_mutex_unlock(&mt[y]);
    }
    long long findfirst(int hash) {
	if (hash < 0) return ((long long) -1);
	return (current = ((long long ) index[hash])-1);
    }
    long long findnext() {
	return (current = ((long long) plist[current])-1);
    }
    long long max_count() { return maxcount;}
    unsigned int *getIndex() { return index;}
    unsigned int *getplist() { return plist;}
 protected:
    unsigned int *index;
    unsigned int *plist;
    long long  current;
    long long  maxcount;
    int derived;
    pthread_mutex_t *mt;
    int num_m; 
};

class HHH {
 public:
    HHH() {;}
    virtual HHH *share(){return NULL;}
    make_hash *setonepattern(const char *pat, _hashtable* &htb, int &m_wsize, int &tzero);
    virtual void setpattern(char *pat, int &ms, int &tailingzero ) {;}
    virtual long long  findfirst(){return 0;}
    virtual long long findnext(){return 0;}
    virtual int addone(char a){return 0;}
    virtual void insert(unsigned pos1){;}
    virtual void setseq(long long length) {;}
    virtual long long max_count() {return 0;}
    virtual void setMutex(pthread_mutex_t *m, int p) {;}
};

class HHH1 : public HHH {
 public:
    HHH1() {htb = NULL; hashitem = NULL;}
    HHH1(_hashtable *h) { htb = new _hashtable(h); hashitem = NULL;}
    virtual HHH *share() {
	HHH1 *jj = new HHH1(htb);
	return jj;
    } 
    virtual void setpattern(char *pat, int &m_wsize, int &tzero) {
	hashitem = setonepattern(pat, htb, m_wsize, tzero);
    }
    virtual void setseq(long long length) {
	htb->setSeqlen(length);
    }
    virtual long long findfirst();
    virtual long long findnext();
    virtual int addone(char a){
	//if (a == AMBIG)return -1; 
	hashitem->addone(a);
	return 1;
    }
    virtual void insert(unsigned int pos1) {
	int h = hashitem->hash();
	//hashnode *n = new hashnode(h, seq, pos1);
	//fprintf(stderr, "sequence %d %d\n", h, pos1);
	htb->insert(pos1, h);	
    }
    virtual long long max_count() { return htb->max_count();}
    virtual void setMutex(pthread_mutex_t *m, int p) {
	htb->setMutex(m, p);
    }
 protected:
    _hashtable *htb;
    make_hash *hashitem;
};

class HHHm : public HHH {
 public:
    HHHm() {currentN = N = 0;}
    HHHm(_hashtable **h, int n) {
	N= n;
	int i;
	for (i= 0; i < N; i++) {
	    htb[i] = new _hashtable(h[i]);
	}
	currentN = 0;
    }
    virtual HHH *share() {
	HHHm *hm = new HHHm(htb, N); 
	return hm;
    }
    virtual void setpattern(char *pat, int &m_wsize, int &tzero);
    virtual long long findfirst();
    virtual long long findnext();
    virtual int addone(char a) {
	//if (a == AMBIG)return -1; 
	int i;
	for (i = 0; i < N; i++) hashitem[i]->addone(a);
	return 1;
    }
    virtual void setseq(long long length) {
	int i;
	for (i = 0; i < N; i++) 
	    htb[i]->setSeqlen(length);
    }
    virtual void setMutex(pthread_mutex_t *n, int p) {
	int i; for (i = 0; i < N; i++) htb[i]->setMutex(n, p);
    }
    virtual void insert(unsigned int pos1) {
	int i;
	for (i = 0; i < N; i++) {
	    int h = hashitem[i]->hash();
	    //hashnode *n = new hashnode(h, seq, pos1);
	    //printf("sequence %d %d\n", hh, pos1);
	    htb[i]->insert(pos1, h); 
	}
    }
 protected:
    _hashtable *htb[100];
    make_hash *hashitem[100];
    int  currentN, N;
};

class matrix {
 public:
    matrix(int w, int aow, int scale) {
	setAOW(aow);
	setScale(scale);	
	matrix_main(w);
    }
    matrix(char *matrix_file, int w, int aow, int scale) {
	setAOW(aow);
        setScale(scale);
	FILE *fp = ckopen(matrix_file, "r");
	matrix_read(fp, w);
    }
    matrix (FILE *fp, int w, int aow) {
	setAOW(aow);
	matrix_read(fp,w);
    }
    void matrix_read(FILE *fp, int w);
    void matrix_main(int w);
    ~matrix();
    int *ma(int, int, int);
    int *mar(int, int, int);
    void setwindow(int a) {wordsize = a;};
    void revising();
    void setScale(int s) { Scale =s;}
    int gap(int i);
    int Scale;
 protected:
    void setAOW(int a) { AOW= a;}
    int num, num2;
    int wordsize;
    int *m;
    int *mr;
    int *m2, *mr2;
    int sym;
    int AOW;
};

class unique_hits {
  public:
    unique_hits() {offset = ((long long) 0);}
    void setOffset(long long s) {offset = s;}
    virtual void addhits(long long pos){;}  
    virtual int addline(char *line, int &x, int len) {return 0;}
    virtual int is_new(long long pos, long long pos_e, int reversed, int dir, long long slen) {return 1;}
    virtual void renew() {;}
  protected:
    long long offset;
};

#define MMLIST 1000

class unique_hits_check : public unique_hits {
  public:
    virtual void addhits(long long pos);
    virtual int is_new(long long pos, long long pos_e,int reversed, int dir, long long slen);
    virtual int addline(char *line, int &x, int len);
    virtual void renew();
  protected:
    int num[16];
    long long list[16][MMLIST];
}; 

class mapping_machine {
 public:
    mapping_machine();
    ~mapping_machine();
    int ReadprobeFile (const char *fname);
    int ProcessSeq (const char *seq_label, const char *seq_data);
    long long init_hash(const char *seq_label, const char *seq_data);
    int build_hashing(long long b, long long e);
    void SetThreshold (int t) ;
    void SetWordSize (char *pat);
    void share(mapping_machine *copy) {
        hash = copy->getHash()->share();
	copy->getSeq(seqseq, seqcolor);
	totallen = strlen(seqseq);
    }
    void getSeq(char *&seq1, char *&seq_c) {
	seq1 = seqseq; seq_c = seqcolor;
    }
    virtual int left_local_align(const char *seq1, const char *seq2, int match, int mismatch, int &len, int dir, int lstate);
    int left_local_align(const char *seq1, const char *seq2, int match, int mismatch, int &len, int dir, int s, int t);
    int left_local_align(int i, const char *seq2, int match, int mismatch, int &len, int dir, int s, int t)  {
	return left_local_align(probe_seq+i, seq2, match, mismatch, len, dir, s, t);
    }
    int num_ext_allowed;
    void setMatrix();
    void setPrefixLen(int len) {prefix_len = len;}
    void getQ();
    void setMatrix(char *);
    void setMatrix(FILE *);
    void setgapsize(int a) {gapsize = a;}
    void setfirstlast(long long f, long long r) {firstr  = f; lastr= r;}
    void probe_init(int len);
    void setdisplay(int a) {display_alig = a;}
    void setTempout(char *t) {if (t) tempout = strsave(t);}
    void setFilter(int f) { filter = f;}
    void setColorspace(int c) {if(c==0) colorspace = 0;else colorspace=1;}
    void setBoth(int b) {both = b;}
    void setDir(int d) {dir = d;}
    void setShort(int s) {shortout = s;}
    void setHitlimit(int h) { hlimit = h;}
    void setNoperfect(int n) {noperfect = n;} // not report perfect hits if 1
    void setReversed(int r) {reversed = r;}
    void setReference(int r) {reference = r;}
    void setcolorcode(char *s); // other di-color code can be set by this
    void setXdrop(int x) {xdrop = x;} // not report hits xdrop more mismatches than the best hit 
    void setMismatch(int  m) {if (m <= 0) mismatch = m;}
    void setQfile(char *f) {  qfile =fopen(f, "r");}
    void setOutputfile(char *f) {
	if (f) {
	    outputfile = strsave(f);
	    outp = ckopen(outputfile, "w");
	}
    }
    void print_line(const char *s) {
	int outsize = fprintf(outp, "%s", s);
	check_printf;
    }
    void printn() {
	int outsize = fprintf(outp, "\n");
	check_printf;
    }
    void moveOuttoTemp() {
	char line[1000];
	sprintf(line, "mv -f %s %s", outputfile, tempout);
	//fprintf(stderr, "%s\n", line);
	system(line);
    } 
    void setUniqueCheck() {delete uc; uc = (unique_hits *) new unique_hits_check();}
    int buildprobe(const char *seq1);
    int header(unique_hits *uc, FILE *tpout, int shortout, int pos1, int filter);
    char _compl[128];
    char _scode[128];
    virtual void probe(char *, int);
    int quality_value(char *qstring);
    void setPsi(const char *seq);
    void setAverageQV(double qv) {averageQV = qv+5.0;}
    void setOffset(long long s) { offset = s;}
    int re_align(char *ref, const char *seq, int threshold) {
	int sum, ll, rr;
	seqseq = ref;
	sum = probe_first[ref[0]]; 
	if (ungappedmatch(seq, 0, threshold, sum, ll, rr, 0) == 1) return sum; 
	else return 10000;
    }
    int re_align(char *ref, const char *seq) {
	return re_align(ref, seq, 10000);
    }
    void getcolorseq(const char *ref, int len, char *color) {
	int i;
	for (i = 0; i< len-1; i++) color[i] = seq2color[ref[i]][ref[i+1]];  
    }
    HHH *getHash() { return hash;}
    void run_classification(const char *seq_data, const char *fname, const char *f1, const char *f2, long long b, long long e);
 protected:
    int  probe_main(char *probe_seq, matrix *pm, int nh);
    int real_length; // real length of the read present, real length may be longer than length, where the additional bases are not used in mapping. 
    int length;//read length in mapping
    int m_wsize;
    int hlimit;
    int noperfect;
    int do_classification;
    int xdrop; 
    matrix *mat;
    int threshold;
    long long offset;
    int amb_use_worst; 
    int gapsize;
    char probe_seq[1000];
    void display_align_core(const char *seq_data, int pos, int pos_e, int left, int right);
    void display_align(const char *seq_data,int pos, int pos_e);
    virtual int match(long long seed, int index, int wscore, int threshold, int &sum1, long long slen);
    virtual int ungappedmatch(const char *seq, long long seed, int threshold, int &sum1, int &ll, int &rr, char islocal);
    int gapmatch(const char *seq, long long & leftp, int threshold, int gs, int align) { return 0;}
    char probeid[100000];
    int display_alig;
    long long ReportHit(const char *seq_label, const char *seq_data, long long pos, long long pos_e, int sum1, long long slen);
    void modify(int, int, int); 
    int **psi, **psi_qvalue;
    char *tempout;
    int filter; //0 no 1 zero filter 2 anything
    int colorspace;
    int reference;
    int dir;
    int reversed;
    FILE *qfile;
    int tailingzero;
    int both;
    int shortout;
    char color2seq[128][4];
    char seq2color[128][128];
    HHH *hash;
    char *sname;
    char *seqseq;
    char *seqcolor;
    long long  totallen;
    int selfscore;
    char first_nuc;
    int *probe_first;
    int mismatch;
    long long firstr, lastr;
    int inc[100], end[100], benum;
    int score_adj;
    unique_hits *uc;
    double averageQV;
    char *outputfile;
    FILE *outp;
    char probe_prefix[1000];
    int prefix_len;
};

class mapping_machine_color : public mapping_machine
{
 public:
    mapping_machine_color();
    virtual void probe(char *, int);
    void setPsiMask(const char *seq, int *allzero, const char *mask);
    virtual int match(long long seed, int index, int wscore, int threshold, int &sum1, long long  slen);
};

class mapping_machine_color_adja : public mapping_machine_color
{
  public:
    mapping_machine_color_adja();
    int match(long long seed, int index, int wscore, int threshold, int &sum1, long long slen);
    virtual int ungappedmatch(const char *seq, long long seed, int threshold, int &sum1, int &ll, int &rr, char islocal);
    char amb_this[128];
};


class mapping_machine_color_adj : public mapping_machine_color_adja
{
  public:
    mapping_machine_color_adj();
    virtual int ungappedmatch(const char *seq, long long seed, int threshold, int &sum1, int &ll, int &rr, char islocal);
    virtual int left_local_align(const char *seq1, const char *seq2, int match, int mismatch, int &len, int dir, int lstate);
    virtual void probe(char *, int);
    char probe_class[1000]; 
    char probe_second; 
    void setSNPrate(double r, int amb);
    int consistent(char a, char b) {
	return (a == amb_other[b]);
    }
 protected:
    char amb_other[128];
    int  score_snp;
};

int localscore(int len, int mis, int mat, int numM);

#endif /* __map_h__ */
