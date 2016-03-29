#include "util.h"
#include "fasta-io.h"
extern double erate;
void mutate(char *seq, int len, double error_rate);
int firstcl(const char *seq, int length);
void init_code(char *binfile, int &nn, int prec, int maxlen);
int get_distance(int nn);
void makeseq(char *f, int s, int len, const char *seq);
void produce_read(const char *seq, int seqlen, int nn, char *f, char *r, int &s1, int &s2);
void makeseqwdel(char *f, int s, int len, const char *seq, int sdel, int ldel);
void produce_del_read(const char *seq, int seqlen, int nn, char *f, char *r, int &s1, int &s2, int maxdel, int ends);
void produce_ins_read(const char *seq, int seqlen, int nn, char *f, char *r, int &s1, int &s2, int maxdel, int ends);

