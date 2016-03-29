#include "map.h"
#include "getComp.h"
#define SSIZE 500000
#define MAXCHROM 1000000

#include <vector>
#include <string>
using std::vector;
using std::string;

class node {
  public:
    long long pos;
    int numMis;
    int chrom;
    int len;
    int score;
    int r_start;
    char ext[150];
    void output(int iomode);
    void output(char *s, int iomode);
} ;

class rescue_basic
{
 public:
    rescue_basic() {threshold = 10000; VA = 0; mode = 0;}
    void setThreshold(int t, int sc) {threshold = t*sc+sc/2; scale = sc;}
    void setThreshold(int t)  { setThreshold(t, 10);}
    void setThreshold(int t, int tr, int tf, int sc) {
	scale = sc;
	threshold = t*sc;
	threshR = tr *sc;
	threshF = tf*sc;
    }
    void setThreshold(int t, int tr, int tf) {
	setThreshold(t, tr, tf, 10);
    }
    void setVA(int a) { VA = a;}
    void setMode(int m) {mode = m;}
    virtual void pair_rescue(node *list1, int n1, node *list2, int n2, int lower, int upper, char *read1, char *read2, int outmode, int len1, int len2i, int mr) {;}
    virtual void pair_rescue(node *list1, int n1, node *list2, int n2, int lower, int upper, char *read1, char *read2, int outmode, int len1, int len2i, int mr,vector<int> &F3_pd, vector<int> &R3_pd) {;}
    int check_score(int score, int &m, int &s) {
	if (score < 0) return 0;
	s = (score & ((1 <<15)-1));
	m = (score >> 15);
	if (s > threshold) return 0;
	return 1;
    }
    void set_unique_only (int u_in) { unique_only=u_in; };
  protected:
    int threshF, threshR, threshold;
    int VA;
    int scale;
    int mode;
	int unique_only;  //1 for unique outputs only in rescue
};

typedef rescue_basic rescue_no;

class rescue : public rescue_basic
{
 public:
    rescue();
    int ProcessSeq(const char *seq_label, const char *seq_data);
    void ProcessFile(char *seqfile);
    int get_refcolor(char *&a, char *&b, int &beg, int &end, int nc, int slen, int &dir);
    virtual void setprobe(const char *id, char *probe, int len);
    virtual int align(int chrom, int beg, int end, int gap, int &l, int &r, int &lscore);
    int align(int chrom, int beg, int end, int gap, int &l, int &r) {
	int ll;
	return align(chrom, beg, end, gap, l, r, ll);
    }
    void setmatching(int adj, float SNPrate, int amb_allow, char *qfile);
    void setmatching(int adj, float SNPrate, int amb_allow, char *qfile, int mis);
    virtual void setMatrix(char *);
    void setNum(int s) {numRes = s;}
    void Fmask(char *s) {mF = s;}
    void Rmask(char *s) {mR = s;}
    void setLocal(int mis_p, int len) {mis = mis_p; min_length = len;}
    virtual void pair_rescue(node *list1, int n1, node *list2, int n2, int lower, int upper, char *read1, char *read2, int outmode, int len, int len2i, int mr);
    virtual void pair_rescue(node *list1, int n1, node *list2, int n2, int lower, int upper, char *read1, char *read2, int outmode, int len, int len2i, int mr,vector<int> &F3_pd, vector<int> &R3_pd);
    mapping_machine_color *mp;
    int get_range_beg () const { return range_beg; };
    int get_range_end () const { return range_end; };
 protected:
    char *seq[MAXCHROM], *seqrev[MAXCHROM];
    char *seqseq[MAXCHROM], *seqseqrev[MAXCHROM];
    matrix *mat;
    int numC;
    char probeid[100000];
    char first;
    int probe_len, real_len;
    char *mask;
    int seqlen[MAXCHROM];
    char _compl[128], _scode[128];
    char color2seq[128][4];
    char seq2color[128][128];
    double S[SSIZE];
    int *psi[1000];
    char probe_seq[1000];
    int numRes;
    int allzero[128];
    char *mF, *mR;
    //void process(node *list1, int n1,int lower, int upper, char *read1, const char *id, int mode, int len, int len2);
    int process(node *list1, int n1,int lower, int upper, char *read1, char *read2, const char *id, int mode, int len, int len2, vector<int> &pairDists);
    char temp_out[10000];  //Put output in here first
    int range_beg, range_end;
    int mis, min_length;
};

class rescue_nogap : public rescue
{
  public:
    virtual int align(int chrom, int beg, int end, int gap, int &l, int &r, int &lscore);
    virtual void setprobe(const char *id, char *probe, int len);
   protected:
    void run_align(int, int &, int &, int &, int &, char *, char *);
    int hash[256];
    int hashitem[256];
    int numhit[SSIZE];
};

class halfalign
{
  public:
    int diag, len;
};
class operation
{
  public:
    operation() { ;}
    virtual void init() {;}
    virtual void add_diag(int d, int col, int mis, int r) {;}
    virtual int find_align(int d, int col, int mis, int &dt, int &colt, int r) {return -1;}
	int halflen;
};


class indel : public operation
{
  public:
    indel() {hf = NULL; header=end = NULL;}
    ~indel();
    virtual void init();
    virtual void add_diag(int d, int col, int mis, int reallen);
    virtual int find_align(int d, int col, int mis, int &dt, int &colt, int reallen);
    void setMaxIndel(int m) {maxIndel=m;}
    void setAll(int rlen, int mMis, int hlen);
  protected:
    int maxIndel;
    int maxMis;
    halfalign **hf;
    int *header, *end;
};

class rescue_oneindel : public rescue
{
  public:
    rescue_oneindel() {consist_check = new getComp();}
    virtual int align(int chrom, int beg, int end, int gap, int &l, int &r, int &lscore);
    void setMaxIndel(int maxIns, int maxDel, int rlen, int mMis, int hi, int hd);
  protected:
    operation *ins;
    operation *del;
    getComp *consist_check;
    void betteralign(int nm, int s, int len, int fdlen, int beg, int end, int &snm, int &start, int &alen, int &fdl);
};

class remapping : public rescue_nogap
{
  public:
    remapping();
    //int quality_value(char *qstring);
    int re_align(int chrom, int beg);
    //void modify(int, int, int);
    //int *psi_qual[1000];
    virtual void setMatrix(char *);
};

int lenofread(char *);
long long opp(long long, int);

int localscore(node *n, int mis, int mat);

int localscore(node *m, int mis);
