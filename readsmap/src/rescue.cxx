#include "rescue.h"
#include "fasta-io.h"
#include "zcompress.h"

#define Aadiffer 'a'-'A'

long long opp(long long x, int len) { return -x-len;}  //must satisfy: opp( opp(x, len), len) = x

static void sets2c(char s[][128], char a, char b, char c)
{
    s[a][b] = s[a+Aadiffer][b] = s[a][b+Aadiffer]= s[a+Aadiffer][b+Aadiffer] = c;
}

rescue::rescue()
{
	unique_only = 0;  //0 (false) <- will report if there is multiple rescues, 1 <- will not
    init_compl(_compl);
    init_scode(_scode);
    mat = NULL;
    mask = mF = mR = NULL;
    numC = 0;
    int i, j;
    for (i = 0; i < 128; i++) {
	allzero[i] = 0;
	for (j = 0; j < 128; j++) seq2color[i][j] = 'N';
        seq2color[i]['.'] = 'X';
    }
    for (i = 0; i < 128; i++) seq2color['.'][i] = 'X';
    allzero['.'] = 10000;
    allzero['X'] = allzero['x'] = 10000;
    sets2c(seq2color, 'A','A', 'A');
    sets2c(seq2color, 'C','C', 'A');
    sets2c(seq2color, 'T','T', 'A');
    sets2c(seq2color, 'G','G', 'A');
    sets2c(seq2color, 'A','C', 'C');
    sets2c(seq2color, 'C','A', 'C');
    sets2c(seq2color, 'G','T', 'C');
    sets2c(seq2color, 'T','G', 'C');
    sets2c(seq2color, 'A','G', 'G');
    sets2c(seq2color, 'G','A', 'G');
    sets2c(seq2color, 'C','T', 'G');
    sets2c(seq2color, 'T','C', 'G');
    sets2c(seq2color, 'A','T', 'T');
    sets2c(seq2color, 'T','A', 'T');
    sets2c(seq2color, 'C','G', 'T');
    sets2c(seq2color, 'G','C', 'T');


    color2seq['A']['0'] = 'A';
    color2seq['A']['1'] = 'C';
    color2seq['A']['2'] = 'G';
    color2seq['A']['3'] = 'T';
    color2seq['C']['0'] = 'C';
    color2seq['C']['1'] = 'A';
    color2seq['C']['2'] = 'T';
    color2seq['C']['3'] = 'G';
    color2seq['G']['0'] = 'G';
    color2seq['G']['1'] = 'T';
    color2seq['G']['2'] = 'A';
    color2seq['G']['3'] = 'C';
    color2seq['T']['0'] = 'T';
    color2seq['T']['1'] = 'G';
    color2seq['T']['2'] = 'C';
    color2seq['T']['3'] = 'A';
}

void rescue::ProcessFile(char *refseq)
{
    FastaFile fafile(SEQTYPE_NT);
    FastaSeq faseq;

    if (!fafile.Open(refseq,"r")) return;
    while (fafile.Read(faseq)) {
	ProcessSeq(faseq.Label(), faseq.Sequence());
    }
}

int rescue::ProcessSeq(const char *seq_label, const char *seq_data)
{
    if (numC >= MAXCHROM) fatal ("too many chromosome\n");
    int len = seqlen[numC] = strlen(seq_data);
    /*
    seq[numC] = new char[len];
    if (seq[numC] == NULL) fatal("out of memory\n");
    seqrev[numC] = new char[len];
    if (seqseq[numC] == NULL)fatal("out of memory\n");
    */
    seqseq[numC] = strsave(seq_data);
    /*
    seqseqrev[numC] = new char [len+1];
    if (seqseqrev[numC]) fatal("out of memory\n");
    reverse(seqseq[numC], len, seqseqrev[numC], _compl);
    int i;
    char *s = seqseq[numC], *t = seqseqrev[numC];
    for (i = 0; i < len-1; i++) {
	seq[numC][i] = seq2color[*s][*(s+1)];
	seqrev[numC][i] = seq2color[*t][*(t+1)];
	s++; t++;
    }
    */
    numC++;
    return 1;
}

void rescue::setMatrix(char *mfile)
{
    if (mat) delete mat;
    if (mfile) {
	mat = new matrix(mfile, 10, 0, 1);
    } else mat = new matrix(0, 9, 1);
}

void remapping::setMatrix(char *mfile)
{
    if (mat) delete mat;
    if (mfile) {
        mat = new matrix(mfile, 10, 0, 10);
    } else mat = new matrix(0, 9, 10);
}

void rescue_nogap::setprobe(const char *id, char *probe, int len)
{
    //strncpy(probeseq, probe+2, len-1);
    char code[] = "ACGTN";
    mask = (id[0] == 'F')? mF : mR;
    strcpy(probeid, id);
    strcpy(probe_seq, probe);
    //probe = probe_seq;
    if (mis < 0) {
	real_len = len;
	probe_len = min_length;
	if (probe_len > real_len) probe_len = real_len;
    } else {
        real_len = len; probe_len = len;
    }
    first = color2seq[toupper(probe[0])][probe[1]];
    mp->probe_init(probe_len-1);
    if (probe[0] == '.' || probe[1] == '.') first = 'N';
    mp->probe(probe, -1);
    if (mask == NULL) mp->setPsi(probe);
    else mp->setPsiMask(probe, allzero, mask);
    mp->getQ();
}

void rescue::setprobe(const char *id, char *probe, int len)
{
    //strncpy(probeseq, probe+2, len-1);
    char code[] = "ACGTN";
    mask = (id[0] == 'F')? mF : mR;
    strcpy(probeid, id);
    strcpy(probe_seq, probe);
    //probe = probe_seq;
    probe_len = len;
    first = color2seq[toupper(probe[0])][probe[1]];
    if (probe[0] == '.' || probe[1] == '.') first = 'N';
    int i;
    if (mask == NULL) {
	for (i = 1; i < len; i++) {
	    if (probe[i+1] == '.') probe[i+1] = '4';
	    psi[i] = mat->ma(i+1, code[probe[i+1]-'0'], len);
    	}
    	psi[0] = mat->ma(1, first, len);
    } else {
        for (i = 1; i < len; i++) {
	    if (mask[i] == '0') psi[i] = allzero;
            else psi[i] = mat->ma(i+1, code[probe[i+1]-'0'], len);
        }
	if (mask[0] == '0') psi[0] = allzero;
	else psi[0] = mat->ma(1,first, len);
    }
}

static char temp1[SSIZE+100], temp2[SSIZE+100];
int rescue::get_refcolor(char * &ref, char *&color, int &beg, int &end,int nc,  int slen, int &dir)
{
    if (beg < 0) {
        ref = temp1;
        dir = -1;
        beg += slen-1;
        end += slen-1;
    } else {
        ref = seqseq[nc];
        dir = 1;
	ref+= beg;
    }
    if (beg < 0) beg = 0;
    if (end >= slen) end = slen-1;
    int len = end-beg+1;
    if (len >= SSIZE) fatal("region size exceed array\n");
    color = temp2;
    if (dir == -1) {
	reverse(seqseq[nc]-end+slen-1, len, ref, _compl);
    }
    int i;
    for (i = 0; i < len; i++) {
	color[i] = seq2color[ref[i]][ref[i+1]];
    }
    return len;
}

void rescue::setmatching(int adjas1, float SNPrate, int amballow, char *qfile) {
	setmatching(adjas1, SNPrate, amballow, qfile,0);
}

void rescue::setmatching(int adjas1, float SNPrate, int amballow, char *qfile, int mis_penalty)  {
	setVA(adjas1);
	mis = mis_penalty;
        if (adjas1 == 0 && SNPrate <= 0.00001 && amballow == 0) {
                mapping_machine_color *mp1 = new mapping_machine_color();
                mp = mp1;
        } else {
                if (adjas1<=1) {
                    mapping_machine_color_adj *mp1 = new mapping_machine_color_adj();
                    if (adjas1 == 1) {
                        mp1->setSNPrate(2.0, amballow);
                    } else {
                        mp1->setSNPrate(SNPrate, amballow);
                    }
                    mp = mp1;
                } else {
                    mapping_machine_color_adja *mp1 = new mapping_machine_color_adja();
                    mp = mp1;
                }
       }
	mp->setQfile(qfile);
	mp->setMismatch(mis);
	char *a = NULL;
	mp->setMatrix(a);
}

void rescue_nogap::run_align(int i, int &min, int &max, int &w, int &e, char *ref, char *color)
{
    char code[] = "ACGT";
    int sum = mp->re_align(ref+i, color+i, threshold);
    if (mis < 0) {
      if (sum <= threshold) {
	int l, score;
	score = mp->left_local_align(probe_len-1, color+i+probe_len-1, 10, mis, l, 1, code[probe_seq[probe_len]-'0'], color[i+probe_len-2]);
        score += (probe_len-1)*10+sum/10*(mis-10);
	//if (score > 0) fprintf(stderr, "sum %d score %d l=%d \n", sum, score, l);  // debug
	if (score > max) {
            max = score;
            min = sum;
            w = i;
            e = l;
        }
      }
    } else if (sum < min ) {
        w = i;
    	min = sum;
    }
}

int rescue_nogap::align(int nc, int beg, int end, int gap, int &left, int &right, int &lscore)
{
    if (nc >= numC) return -1;
    if (real_len < probe_len) return -1;

    char *color;
    char *ref;
    int dir;
    int min = 10000, w, max = 0, e;
    int mMis = threshold/scale;
    int slen = seqlen[nc];
    int len = get_refcolor(ref, color, beg, end, nc, slen, dir);
    mp->getcolorseq(ref, len+1, color);
    //fprintf(stderr, "beg%d end%d len %d probe_len %d min_len %d Mis %d\n",beg, end,len,  probe_len, min_length, mMis);
    //fprintf(stderr, "probe_seq %s\n ref %s\n", probe_seq, ref);
    if (len < probe_len) return -1;
    e = 0;
    int hs;
    hs = 4;
    int hsv = hs+VA;
    int limit = (probe_len-hs) - hsv*mMis;
    if (limit <=0 && probe_len <=32) {
	hs = 3;
	hsv =hs+VA;
 	limit = (probe_len-hs) - hsv*mMis;
    }
    if (mask == NULL && limit > 0) {
	int hsize = (1 << (2*hs));
	int mask = hsize-1;
	memset(numhit, 0, sizeof(int)*(len+1));
	memset(hash, 0, sizeof(int)*hsize);
	memset(hashitem, 0, sizeof(int)*probe_len);
	int i;
	int h = 0, N=hs-1;
	for (i = 1; i < probe_len; i++) {
	    int x = probe_seq[i+1]-'0';
	    if (x == 4) N = hs-1;
	    else {
		h = (((h << 2) & mask) | x);
		if (N > 0) { N--; }
		else {
		    hashitem[i] = hash[h];
		    hash[h] = i;
		}
	    }
	}
	N = hs-1; h = 0;
	int ii;
	for (ii = 0; ii < len; ii++) {
            int x = _scode[color[ii]];
            if (x == AMBIG) N = hs-1;
            else {
                h = (((h << 2) & mask) | x);
                if (N > 0) { N--;}
                else {
		    int p;
		    for (p = hash[h]; p; p = hashitem[p])  {
			int start = ii-p+1;
			if (start < 0 || start >= len-real_len) continue;
			numhit[start]++;
			if (numhit[start] != limit) continue;
			//extension
        		int k = 0;
			i = start;
			/*
        		if ((sum = psi[0][ref[i]]) > 0) {
            		    if (mMis == 0) continue;
            		    k++;
        		}
			*/
        		if (ref[i] == '.') break;
			/*
			int j;
        		for (j = 1; j < probe_len;j++) {
		            int a;
		            if ((a =psi[j][color[i+j-1]]) > 0) {
		                if (color[i+j-1] == 'X' ) break;
		                k++;
		                //if (k == 3 && j < 6) break;
		                if (k > mMis) break;
		                sum += a;
		            }
		        }
		        if (j >= probe_len && sum < min ) {
		            w = i;
		            min = sum;
		        }
			*/
			run_align(i, min, max, w, e, ref, color);
		    }
		}
	    }
	}

    	left = beg+w;
    	right = left+probe_len+e-2;
	lscore = max;
    	if (dir < 0) {
            left -= slen-1;
            right -= slen-1;
        }
	return min;
    }
    int i, j;
    for (i = 0; i < len-real_len; i++) {
        int sum = 0;
        int k = 0;
        if (ref[i] == '.') break;
	run_align(i, min, max, w, e, ref, color);
    }
/*
    int i, j;
    for (i = 0; i < len; i++) {
        S[i] = psi[0][ref[i]];
    }
    for (j = 1; j < probe_len; j++) {
        double  f=S[j-1];
	S[j-1] = 100000;
	int *p = psi[j];
	char *c = color+j-1;
        for (i = j-1; i < len; i++) {
            double x = f + p[*c++];
            f = S[i+1];
            S[i+1] = x;
        }
    }
    double min = 101000;
    int w = 0;
    for (i = 1; i <= len; i++) {
        if (min > S[i]) {
            min = S[i];
            w = i-1;
        }
    }
*/
    left = beg+w;
    right = left+probe_len+e-2;
    lscore = max;
    if (dir < 0) {
        left -= slen-1;
        right -= slen-1;
    }
    return (int) min;
}

remapping::remapping()
{
    /*
    int i;
    psi_qual[0] = new int[128*1000];
    for (i = 1; i < 1000; i++)  {
	psi_qual[i] = psi_qual[i-1]+128;
    }
    */
}

int remapping::re_align(int nc, int beg)
{
    if (nc >= numC) return -1;
    char *color;
    char *ref;
    int dir;
    int slen = seqlen[nc];
    int end = beg+probe_len-1;
    int len = get_refcolor(ref, color, beg, end, nc, slen, dir);
    mp->getcolorseq(ref, len+1, color);
    if (len < probe_len) return -1;
    return mp->re_align(ref, color);
    	/*
        int sum = psi[0][ref[0]];
	int j;
	//printf("(%d ", sum);
        for (j = 1; j < probe_len;j++) {
	    //printf("%d ", psi[j][color[j-1]]);
	    sum+=psi[j][color[j-1]];
        }
	//printf(")");
	return sum;
	*/
}
/*
static double q2p[] = {
    6.0,
    6.0,
4.0,
3.0,
2.2,
1.6,
1.24,
1.0,
0.7,
0.55,
0.4500000000000000,
0.36,
0.26,
0.2,
0.2,
0.1,
0.1,
0.1,
0.0,
0.0,
0.000000000000000,
0.00
};

void remapping::modify(int x, int i, int index)
{
    char a[]="ACGT";
    char c = a[probe_seq[i+1]-'0'];
    if (i == 0) c=first;
    int j, *old = psi[i], sum;
    double r, p;
    int *nw  = psi[i] = psi_qual[index];
    sum = old['A']+old['C']+old['G']+old['T'];
    if (x == 0) x= 10;
    p = q2p[x]/20.0;//  prob of wrong call
    r = ((double)(x+5))/((double) 20);
    //printf("%d %f %f\n", x, p, r);
    for (i = 0; i < 4; i++) {
        nw[a[i]] = (int) (r*old[a[i]]);
    }
    nw[c] = (int) (p*sum/3); // match value
    //printf("%d %d\n", nw[c], old[c]);
    nw['N'] = old['N'];
    nw['X'] = old['X'];
}


int remapping::quality_value(char *qline)
{
    int i, index_tempmem = 0;
    char *q = qline-1;
        for (i = 0; i < probe_len; i++) {
                if (q == NULL) break;
                q++;
                int x = atoi(q);
                q = strchr(q, ' ');
                if (x > 20) {x=20;}
                modify(x, i, index_tempmem);
                index_tempmem++;
	}
    if (index_tempmem>0) return 0;
    return 1;
}
*/
int rescue::align(int nc, int beg, int end, int gap, int &left, int &right, int &lscore)
{
    if (nc >= numC) return -1;
    char *color;
    char *ref;
    int dir;
    int slen = seqlen[nc];
    int len = get_refcolor(ref, color, beg, end, nc, slen, dir);
    /*
    if (beg < 0) {
        color = seqrev[nc];
        ref = seqseqrev[nc];
        dir = -1;
        beg += slen-1;
        end += slen-1;
    } else {
        color = seq[nc];
        ref = seqseq[nc];
        dir = 1;
    }
    if (beg < 0) beg = 0;
    if (end >= slen) end = slen-1;
    int len = end-beg+1;
    if (len >= SSIZE) fatal("region size exceed array\n");
    */
    if (len < probe_len) return -1;
    int i, j;
    S[0] = gap;
    for (i = 0; i < len; i++) {
	S[i] = psi[0][ref[i]]+((double) i / (double)len);
    }
    for (j = 1; j < probe_len; j++) {
	double e, f=S[0];
	S[0] += gap;
	e = S[0];
	for (i = 0; i < len; i++) {
	    double x = f + psi[j][color[i]];
	    f = S[i+1];
	    if (e > f) e = f;
	    e += gap;
	    if (e > x) e = x;
	    S[i+1] = e;
	}
    }
    double min = 101000;
    int w = 0;
    for (i = 1; i <= len; i++) {
	if (min > S[i]) {
	    min = S[i];
	    w = i-1;
	}
    }
    left = beg+ (int) ((min -(int)min)*len+0.5);
    right = w+beg;
    if (dir < 0) {
	left -= slen-1;
	right -= slen-1;
    }
    return (int) min;
}

void rescue::pair_rescue(node *list1, int n1, node *list2, int n2, int lower, int upper, char *read1, char *read2, int outputmode, int len1, int len2, int mr) {
  vector<int> F3_pd, R3_pd;
  F3_pd.reserve(n1);  //For performance--
  R3_pd.reserve(n2);  //Make array the right size
  int count = 0;

  if (n1 < mr) count += process(list1, n1, lower, upper, read2, read1, "F",  //read2 is F3 tag
		       outputmode, len2, len1, F3_pd);
  if (n2 < mr) count += process(list2, n2, lower, upper, read1, read2, "R",  //read1 is R3 tag
		     outputmode, len1, len2, R3_pd);

  if (count==1 && unique_only)
  	print_line(temp_out);
}


void rescue::pair_rescue(node *list1, int n1, node *list2, int n2, int lower, int upper, char *read1, char *read2, int outputmode, int len1, int len2, int mr,
			 vector<int> &F3_pd, vector<int> &R3_pd)
{
  int count = 0;

  if (n1 < mr) count += process(list1, n1, lower, upper, read2, read1, "F",  //read2 is F3 tag
		       outputmode, len2, len1, F3_pd);
  if (n2 < mr) count += process(list2, n2, lower, upper, read1, read2, "R",  //read1 is R3 tag
		     outputmode, len1, len2, R3_pd);

  if (count==1 && unique_only)
  	print_line(temp_out);
  else {
	  F3_pd.clear(); //no output, so no pairing distances
	  R3_pd.clear();
  }
}


int rescue::process(node *list, int n,int lower, int upper, char *readseq, char *readSeq2, const char *id, int outputmode, int len, int leno, vector<int> &pairDists )
{
    if (n <=0) return 0;
    int save_threshold = threshold;
    int dir = 1;
    int count= 0;
    int thresh = threshF;
    if (id[0] == 'R') { thresh = threshR; dir = -1;}
    setprobe(id, readseq, len);
    int i;
    for (i = 0; i < n; i++) {
	if (abs(list[i].pos) > seqlen[list[i].chrom]) {
	    fprintf(stderr, "chrom %d len%d slen%d Hit position out of bound, reference sequence may be wrong\n", list[i].chrom, list[i].pos, seqlen[list[i].chrom]);
	    exit(1);
	}
	int beg, end;
	long long pos = list[i].pos;
	if (mode == 1) pos = opp(pos, leno);
	if (mode != 1 && dir == 1) {
            beg = pos+lower-len;
            end = pos+upper;
            if (pos < 0 && beg >0) continue;
	} else {
	    beg = pos+leno-upper;
	    end = pos+leno-lower+len;
	    if (pos >= 0) {
		if (end < 0)  continue;
		if (beg < 0) beg = 0;
	    }
	}
	int score, l, r, m, ss;
 	threshold = save_threshold - list[i].numMis*scale;
	if (threshold < 0) continue;
	if (threshold > thresh) threshold = thresh;

	int score1 = align(list[i].chrom, beg, end, 10, l, r, ss);
	if (check_score(score1, m, score) == 0)  continue;
	count++;

	//Next small section to calculate paired distance
	int f3len = (id[0]=='F'?len:leno);
	int pairDist = pos - l; //n2->pos - n->pos;
	if(pairDist<0) { // On neg strand?
	  pairDist = f3len - pairDist;
	} else {
	  pairDist += f3len;
	}
	pairDists.push_back(pairDist);
	//End section

	char readseqNoCR[1002], readSeq2NoCR[1002];
	if (mis < 0) {
	    ss = (ss-(r-l+1)*10)/(mis-10);
	}
	if (m > 0 && mis >= 0) {
	  readseqNoCR[0] = '[';
	  strcpy ( readseqNoCR+1, readseq );
	  size_t CRPos = strlen(readseq);
	  readseqNoCR[CRPos] = ']';

	  readSeq2NoCR[0] = '[';
	  strcpy ( readSeq2NoCR+1, readSeq2 );
	  size_t otherCRPos = strlen(readSeq2);
	  readSeq2NoCR[otherCRPos] = ']';
	}
	//int outsize;
	if (outputmode == 0) {
	if (dir == 1)
	    if (mis < 0) {
		sprintf(temp_out, ",%d.%d%s|%d.%d:(%d.%d.%d)", list[i].pos, list[i].numMis, list[i].ext, l, score/10, r-l+1, ss, 0); // format to be decided.
	    } else if (m > 0)
	    	sprintf(temp_out,",%d.%d%s|%d.%d.%d(%d:%d_%d)%s",
				 list[i].pos, list[i].numMis, readSeq2NoCR,
				 l, r-l+1, score/10, m,
				 range_beg, range_end, readseqNoCR);
	    else
	    	sprintf(temp_out,",%d.%d|%d.%d.%d", list[i].pos, list[i].numMis, l, r-l+1, score/10);
	else // dir == 1
	    if (mis < 0) 
		sprintf(temp_out, ",%d.%d:(%d.%d.%d)|%d.%d%s",l, score/10, r-l+1, ss, 0,list[i].pos, list[i].numMis, list[i].ext);  
	    else if (m > 0)
	    	sprintf(temp_out,",%d.%d.%d(%d:%d_%d)%s|%d.%d%s",
				 l, r-l+1, score/10, m, range_beg, range_end, readseqNoCR,
				 list[i].pos, list[i].numMis, readSeq2NoCR);
	    else
	    	sprintf(temp_out,",%d.%d.%d|%d.%d", l, r-l+1, score/10,list[i].pos, list[i].numMis);
	} else { // outputmode 
        if (dir == 1)
	    if (mis < 0)
                sprintf(temp_out, ",%d_%d.%d%s|%d_%d.%d:(%d.%d.%d)",list[i].chrom+1,list[i].pos, list[i].numMis, list[i].ext, list[i].chrom+1,l, score/10, r-l+1, ss, 0); 
	    else if (m > 0)
	    	sprintf(temp_out,",%d_%d.%d%s|%d_%d.%d.%d(%d:%d_%d)%s", list[i].chrom+1,list[i].pos, list[i].numMis, readSeq2NoCR, list[i].chrom+1, l, r-l+1, score/10, m, range_beg, range_end, readseqNoCR);
	    else
	    	sprintf(temp_out,",%d_%d.%d|%d_%d.%d.%d", list[i].chrom+1,list[i].pos, list[i].numMis,list[i].chrom+1, l, r-l+1, score/10);
        else
            if (mis < 0)
                sprintf(temp_out, ",%d_%d.%d:(%d.%d.%d)|%d_%d.%d%s", list[i].chrom+1,l, score/10, r-l+1, ss, 0,list[i].chrom+1,list[i].pos, list[i].numMis, list[i].ext); 
	    else if (m > 0)
	    	sprintf(temp_out,",%d_%d.%d.%d(%d:%d_%d)%s|%d_%d.%d%s", list[i].chrom+1,l, r-l+1, score/10, m, range_beg, range_end, readseqNoCR, list[i].chrom+1, list[i].pos, list[i].numMis, readSeq2NoCR);
	    else
	    	sprintf(temp_out,",%d_%d.%d.%d|%d_%d.%d", list[i].chrom+1,l, r-l+1, score/10,list[i].chrom+1, list[i].pos, list[i].numMis);
	}
	//check_printf;
	if (unique_only == 0) printf("%s", temp_out);
	if (count >= numRes) break;
    } // end of for loop

    threshold = save_threshold;
	return count;
}

void rescue_oneindel::betteralign(int nm, int s, int len, int fdlen, int beg, int end,int &snm, int &start, int &alen, int &fdl){
    //printf("better %d %d %d\n", nm, s, len);
    // check consistency of indel
    if (nm >= snm) return;
    snm = nm;
    start = s;
    alen = len;
    fdl = fdlen;
    range_beg = beg;
    range_end = end;
}

int rescue_oneindel::align(int nc, int beg, int end, int gap, int &left, int &right, int &lscore)
{
    if (nc >= numC) return -1;
    char *color;
    char *ref;
    int dir;
    int fdlen = 0;
    int start, alen = -1, nn = 1000;
    int mMis = threshold/10;
    int slen = seqlen[nc];
    ins->init();
    del->init();
    nn = mMis+1;
    int len = get_refcolor(ref, color, beg, end, nc, slen, dir);
/*
    if (beg < 0) {
        color = seqrev[nc];
        ref = seqseqrev[nc];
        dir = -1;
        beg += slen-1;
        end += slen-1;
    } else {
        color = seq[nc];
        ref = seqseq[nc];
        dir = 1;
    }
    if (beg < 0) beg = 0;
    if (end >= slen) end = slen-1;
    int len = end-beg+1;
    if (len >= SSIZE) fatal("region size exceed array\n");
    //if (len < probe_len) return -1;
    color += beg;
    ref += beg;
*/
    int i, j;
    for (i = 0; i < len; i++) {
	int sum = 0;
	int k = 0;
	if ((sum = psi[0][ref[i]]) > 0) {
	    if (mMis == 0) continue;
	    k++;
	}
	if (ref[i] == '.') break;
	int endp = probe_len;
	if (endp+i> len) endp -= endp+i-len-1;
	for (j = 1; j < endp;j++) {
	    int a;
	    if ((a =psi[j][color[i+j-1]]) > 0) {
		if (color[i+j-1] == 'X' ) break;
		k++;
		sum += a;
		if (k*3 < (j-1)*2) {//long enough to enter data structure
		    int dd, dlen, nm;
                    del->add_diag(i, j, k-1, j);
		    int x;
                    if ((nm = ins->find_align(i, len+probe_len-j-1-i, k-1, dd, dlen, j))>=0) {
			int kk, t = dlen-len-probe_len+j+1+i, insert = i-dd;
			if (t > j-ins->halflen) t = j-ins->halflen;
			kk = j-t;
			for (x=kk; kk <= j; kk++) {
			    if (consist_check->is_compatible(probe_seq+kk+1,1+insert, color+i+kk-1, 1)) {
				x = kk;
				break;
			    }
			}
			if (kk > j) nm++;
		    	betteralign(nm+k-1, i, probe_len-1-i+dd, x-1, j-t-1, j-1,  nn, start, alen, fdlen);
		    }
		}
		if (k > mMis) break;
	    }
	}
	if (j == probe_len) {
	    // process alignment without indel
	    //printf("%d %d %d\n", k, i, probe_len-1);
	    betteralign(k-1, i, probe_len-1, 0, 0, 0, nn, start, alen, fdlen);
	}
	if (endp < probe_len) continue;
	sum = k = 0;
	for (j = probe_len-1; j >= 1; j--) {
	    int a;
            if ((a =psi[j][color[i+j-1]]) > 0) {
                k++;
		if (color[i+j-1] == 'X') break;
                sum += a;
		int reallen = probe_len-j;
                if (k*3 < (reallen-1)*2) {//long enough to enter data structure
		    ins->add_diag(i, len+probe_len-1-j-i, k-1, reallen);
		    int dd, dlen, nm;
                    if ((nm=del->find_align(i, j, k-1, dd, dlen,  reallen))>=0) {
			if (dd != i) {
			    int t = dlen-j, dele = i-dd, kk;
			    int x = j;
			    if (x < del->halflen) x = del->halflen;
			    if (dlen > 	probe_len-del->halflen) t = probe_len-del->halflen-j;
			    for (kk = x; kk <= j+t; kk++) {
				if (consist_check->is_compatible(probe_seq+kk+1, 1, color+i+kk-1-dele, 1+dele)) {
				    x = kk;
				    break;
				}
			    }
			    if (kk > j+t) nm++;
			    betteralign(nm+k-1, dd, probe_len-1+i-dd, x-1, j-1,j+t-1,nn, start, alen, fdlen);
			}
		    }
                }
                if (k > mMis) break;
            }
	}
    }
exitlabel:
    if (alen >= 0) {
	left = start+beg; right = left+alen-1;
        if (dir < 0) {
            left -= slen-1;
            right -= slen-1;
    	}
	return nn*scale+(fdlen<<15);
    }
    return -1;
}

void rescue_oneindel::setMaxIndel(int maxIns, int maxDel, int rlen, int mMis, int hi, int hd)
{
    if (maxIns == 0) {
	ins = new operation();
    } else {
	indel *i = new indel();
	i->setMaxIndel(maxIns);
	ins = i;
	i->setAll(rlen, mMis, hi);
    }
    if (maxDel == 0) {
	del = new operation();
    } else {
	indel *d = new indel();
	d->setMaxIndel(maxDel);
	del = d;
	d->setAll(rlen, mMis, hd);
    }
}

void indel::setAll(int rlen, int mMis, int hlen)
{
    halflen = hlen;
    header = new int[mMis+1];
    end = new int[mMis+1];
    maxMis = mMis;
    hf = new halfalign*[mMis+1];
    hf[0] = new halfalign[rlen*(mMis+1)];
    int i;
    for (i= 1; i<= mMis;i++) {
	hf[i] = hf[i-1]+rlen;
    }
}

void indel::init()
{
    memset(header, 0, sizeof(int)*(maxMis+1));
    int i;
    for (i = 0; i <=  maxMis; i++) end[i] = -1;
}

indel::~indel()
{
    if (header) delete [] header;
    if (end) delete [] end;
    if (hf) {
	if (hf[0]) delete[] hf[0];
	delete [] hf;
    }
}

int indel::find_align(int d, int len, int mis, int &dt, int &hlen, int reallen)
{
    if (reallen < halflen) return -1;
    int m = maxMis-mis;
    int i;
    for (i = 0; i <= m; i++) {
	halfalign *pn = hf[i];
	int j;
	for (j = header[i]; j <= end[i]; j++) {
	    if (pn[j].diag < d-maxIndel) header[i]++;
	    else break;
	}
	if (j <= end[i] && pn[j].len >= len) {
	    dt = pn[j].diag;
	    hlen = pn[j].len;
	    return i;
	}
    }
    return -1;
}

void indel::add_diag(int d, int len, int mis, int reallen)
{
    if (reallen < halflen) return;
    halfalign *pn = hf[mis];
    int i;
    for (i = end[mis]; i >= header[mis]; i--) {
	if (pn[i].len > len) break;
    }
    end[mis] = i+1;
    pn[i+1].diag = d;
    pn[i+1].len = len;
}

int lenofread(char *a)
{
    int i = strcspn(a, "\r\n")-1;
    return i;
}

int localscore(node *n, int mis, int mat)
{
    return localscore(n->len, mis, mat, n->score);
}

int localscore(node *m, int mis)
{
    return(m, mis, 10);
}

void node::output(int iomode)
{
    char temp[1000];
    output(temp, iomode);
    printf("%s", temp);
}

void node::output(char *s, int iomode)
{
    if (iomode == 0) {
	sprintf(s, "%lld.%d%s", pos, numMis, ext);
    } else {
	sprintf(s, "%d_%lld.%d%s", chrom, pos, numMis, ext);
    } 
}
