#include "fasta-io.h"
#include "rescue.h"
#include "zcompress.h"
#include <queue>

#include "pairTags.h"

#define MAXSIZE 100000
#define MSIZE 20

static int flag_unique = 0;
static int flag_rescue_unique = 0;
static int flag_rescue_was_set = 0;
static int discont_input = 0;
static char temp_out[10000];
static int LL[MSIZE], UU[MSIZE], mnum = 0;
static int max_out =100000;
static int mode = 0;
static int bestscore = 0, secondbest = 0;
static int scoredrop = 100000000;
static int min_length = 10000;
static  int mis_penalty = 0;
static int indelcut = 0;


static int startPanelNum = 0, endPanelNum = 1000000000;
static std::queue<int> panelsToDo;
bool usePanelList = false;

PairingDistanceStat pairingDistances;

int compare_id(char *id, char *id1, int order);

int same_sign(int x, int y)
{
    if (x < 0 && y > 0) return  0;
    return 1;
}

static void getLR(char *a, int &lower, int &upper)
{
    char *p = a-1;
    LL[0] = lower; UU[0] = upper; mnum++;
    do {
	int x, y;
	p++;
	sscanf(p, "(%d,%d", &x, &y);
	int outsize = printf("%d %d\n", x, y);
	check_printf;
	if (x >= y) continue;
	if (mnum >= MSIZE) fatal("Out of memory for LL UU array\n");
	LL[mnum] = x;
	UU[mnum] = y;
	if (x < lower) lower = x;
	if (y > upper) upper = y;
	mnum++;
    } while(p = strchr(p, ';'));
}

int compare(const void *p, const void *q)
{
    node *a = (node *) p;
    node *b = (node *) q;
    if (a->chrom < b->chrom) return -1;
    if (a->chrom > b->chrom) return 1;
    if (a->pos < b->pos) return -1;
    if (a->pos > b->pos) return 1;
    return 0;
}

int readin(node *list, int &n, genFile *gF_fp, char *id, char *readseq, int order)
{
    FILE *fp = gF_fp->getFILE();
    char line[1000000];
    do {
	if (fgets(line, sizeof line, fp) == NULL) return 0;
	if (line[0] != '>') continue;
	int p = panel_number(line);
	//if (p < startPanelNum) continue;  done with seekToPanelId

	if (usePanelList) { //if a list is given, which is sorted by the program,
	  if (p > panelsToDo.front()) {  //use those only in list are done.
	    panelsToDo.pop();
	    if (panelsToDo.empty()) return 0; //nothing more to do
	    gF_fp->seekToPanelId(panelsToDo.front());
	  }
	}
	if (p > endPanelNum) return 0;
	break;
    } while (1);
    n = 0;
    char *p = line;
    int a = compare_id(id, line, order);
    if (a > 0) {
	fprintf(stderr, "At least one of the read file is out of order\n%s%s\n", id, line);
	exit(1);
    }
    while (p = strchr(p, ',')) {
	if (n >= MAXSIZE) fatal("exceed memory\n");
	p++;
	long long a;
	int c = 1, b;
	int start = 0, len = 0, score = -1;
	if (discont_input == 0) sscanf(p, "%lld.%d:(%d.%d.%d", &a, &b, &len, &score, &start);
	else sscanf(p, "%d_%lld.%d:(%d.%d.%d", &c, &a, &b,&len, &score, &start);
	list[n].chrom = c-1;
	list[n].pos = a;
	list[n].numMis = b;
	list[n].r_start = start;
	list[n].score = score;
	list[n].len = len;
	if (score >=0) {
	    sprintf(list[n].ext, ":(%d.%d.%d)", len, score, start);
	    list[n].pos -= start;
	} else list[n].ext[0] = 0;
	n++;
    }
    if (n == 0) {
	sscanf(line, "%s", id);
    } else {
	p = strchr(line, ',');
	*p = 0;
	strcpy(id, line);
    }
    if (fgets(readseq, 1000, fp) == NULL) return 0;
    return 1;
}

int outputpair(node *n, node *n2, int flag_u, char s, int len1, int len2, int &pairDist)
{
    if (n->pos <0 && n2->pos>0) return 0;

    pairDist = n2->pos - n->pos;
    pairDist += len2;

    if(discont_input) {
      if(n->chrom!=n2->chrom) pairDist = PAIR_DIST_NA;
    }

    int nc1,  nm1, nc2, nm2;
    long long np1, np2;
    node *nn1 = n, *nn2 = n2;
    int L1 = len1;
    if (flag_u == 2) {
	nc1 = n2->chrom+1; nc2 = n->chrom+1;
	np1 = n2->pos;     np2 = n->pos;
	nm1 = n2->numMis;  nm2 = n->numMis;
	L1 = len2; nn1 = n2; nn2 = n;
    } else if (flag_u == 3) {
        nc1 = n->chrom+1; nc2 = n2->chrom+1;
        np1 = n->pos;     np2 = opp(n2->pos, len2);
        nm1 = n->numMis;  nm2 = n2->numMis;
    } else if (flag_u == 4) {
        nc1 = n2->chrom+1; nc2 = n->chrom+1;
        np1 = n2->pos;     np2 = opp(n->pos, len1);
        nm1 = n2->numMis;  nm2 = n->numMis;
	L1 = len2; nn1 = n2; nn2 = n;
    } else {
        nc1 = n->chrom+1; nc2 = n2->chrom+1;
        np1 = n->pos;     np2 = n2->pos;
        nm1 = n->numMis;  nm2 = n2->numMis;
    }
    if (mode == 1) np2 = opp(np2,  L1);

    int outsize;
    np1 += nn1->r_start;
    np2 += nn2->r_start;

    if (flag_u == 1) {
	if (mis_penalty != 0 && nn1->score >= 0 && nn2->score >= 0) {
	    int s = localscore(nn1->len+nn2->len, 10, mis_penalty, nn1->score+nn2->score); 
	    if (s < bestscore) { 
		if (s > secondbest) secondbest = s;
		return 1;
	    } 
	    secondbest = bestscore;
	    bestscore = s; 
	}
	if (discont_input)
	    sprintf(temp_out, ",%d_%lld.%d%s%c%d_%lld.%d%s", nc1,np1, nm1,nn1->ext, s, nc2, np2, nm2, nn2->ext);
	else
	    sprintf(temp_out, ",%lld.%d%s%c%lld.%d%s", np1, nm1, nn1->ext,s, np2, nm2, nn2->ext);
    } else {
	if (discont_input)
	    outsize = printf(",%d_%lld.%d%s%c%d_%lld.%d%s", nc1,np1, nm1, nn1->ext, s, nc2, np2, nm2, nn2->ext);
	else
	    outsize = printf(",%lld.%d%s%c%lld.%d%s", np1, nm1, nn1->ext,s, np2, nm2, nn2->ext);
	check_printf;
    }
    return 1;
}

int outputpair(node *n, node *n2, int flag_u, char s, int len1, int len2)
{
  int dummy;
  return outputpair(n,n2,flag_u,s,len1,len2,dummy);
}

int merge(node *list1, int n1, node *list2, int n2, int lower, int upper, int &i_index, int &j_index, char sep, int flag_unique, int len1, int len2,
	  vector<int> &nonResc_pairDist)
{
    int i, j, k, count = 0;
    i = j = 0;
    i_index = j_index = -1;
     bestscore = secondbest = 0;

    if (mode ==0) { 
	if (lower < len1+len2) lower = len1+len2;
	if (upper < lower) return 0;
    }
    while (j < n2) {
	while (i < n1 && list1[i].chrom< list2[j].chrom) i++;
	while (i < n1 && list1[i].chrom ==  list2[j].chrom &&
		list1[i].pos < list2[j].pos+len2-upper) {
		if (same_sign(list1[i].pos, list2[j].pos)) {
		    i_index = i;
		    j_index = j;
		}
		i++;
	}
	k = i;
	while (k < n1 && list1[k].chrom == list2[j].chrom && list1[k].pos <= list2[j].pos+len2-lower) {
	  int pairDist;
	    if (mnum > 0)  {
		int x = -list1[k].pos+list2[j].pos+len2;
		int a;
	        for (a = 0; a < mnum; a++) {
		    if (x > LL[a] && x <= UU[a]) break;
	        }
	        if (a < mnum && outputpair(list1+k, list2+j, flag_unique, sep, len1, len2, pairDist)) {
		     count++;
		     if (mis_penalty == 0) {
		     	if (flag_unique > 1) return count;
		     	if (count > max_out && mis_penalty == 0) goto label_exit;
		    }
		}
	    } else if (outputpair(list1+k, list2+j, flag_unique, sep, len1, len2, pairDist)) {
		count++;
		if (mis_penalty == 0) {
		    if (flag_unique > 1 && mis_penalty == 0) return count;
		    if (count > max_out) goto label_exit;
		}
	    }
	    k++;
	    nonResc_pairDist.push_back(pairDist);

	}
	if (i_index < 0 && k < n1 &&list1[k].chrom == list2[j].chrom &&  list1[k].pos < list2[j].pos && same_sign(list1[k].pos, list2[j].pos)) {
	    i_index = k; j_index = j;
	}
	j++;
    }
label_exit:
    if (flag_unique) {
	if (count == 1 || (count > 1 && bestscore-secondbest > scoredrop)) 
      	    print_line(temp_out);
    } else if (flag_unique) {
      nonResc_pairDist.clear();
    }

    return count;
}

int merge(node *list1, int n1, node *list2, int n2, int lower, int upper, int &i_index, int &j_index, char sep, int flag_unique, int len1, int len2) {

  vector<int> dummy;  //maintain backwards call compatibility.
  return merge(list1, n1, list2, n2, lower, upper, i_index, j_index, sep, flag_unique,len1,len2, dummy);
}

void copy_n(node *n, node *p)
{
    if (n == p) return;
    memcpy(n, p, sizeof(node));
	/*
    n->chrom = p->chrom;
    n->pos = p->pos;
    n->numMis = p->numMis;
	*/
}

void mini_score(node *list, int &n)
{
    if (n <= 0) return;
    int nn =1, i;
    if (mis_penalty < 0) {
	int max_s = localscore(list+0, mis_penalty);
	list[0].score = max_s;
	for (i = 1; i < n; i++) {
	    int sc = localscore(list+i, mis_penalty);
	    list[i].score = sc;
	    if (sc > max_s) max_s = sc;
	}
	nn = 0;
	max_s -= scoredrop/2;
	for (i = 0; i < n; i++) {
	    if (list[i].score >= max_s) {
		copy_n(list+nn, list+i);
		nn++;
	    }
	}
	n = nn;
	return;
    } 
    int min_s = list[0].numMis;
    for (i = 1; i < n; i++) {
	if (list[i].numMis < min_s) {
	    copy_n(list, list+i);
	    nn = 1;
	    min_s = list[i].numMis;
	} else if (list[i].numMis == min_s) {
	    copy_n(list+nn, list+i);
	    nn++;
	}
    }
    n = nn;
}

int min_mis(node *list, int start, int n, int &last)
{
    last = start+1;
    int c = list[start].chrom;
    int minm = list[start].numMis;
    while (last < n) {
	if (list[last].chrom > c) break;
	if (list[last].numMis < minm) minm = list[last].numMis;
	last++;
    }
    return minm;
}

void mov(node *list, int &to, int from, int end, int mis)
{
    int i;
    for (i = from; i < end; i++) {
	if (list[i].numMis == mis) {
	    copy_n(list+to, list+i);
	    to++;
	}
    }
}

void mini_score(node *list1, int & n1, node *list2, int &n2)
{
    int nn1, nn2, i, j;
    int min_s = 100000;
    nn1 = nn2 = i = j= 0;
    while (1)  {
    	while (i < n1 && list1[i].chrom < list2[j].chrom) i++;
	if (i >= n1) break;
	while (j < n2 && list1[i].chrom > list2[j].chrom) j++;
	if (i >= n1 || j >= n2) break;
	if (list1[i].chrom != list2[j].chrom) continue;
	int m1, m2;
	int x = min_mis(list1, i, n1, m1);
	int y = min_mis(list2, j, n2, m2);
	if (x+y <= min_s) {
	    if (x+y < min_s) {
	    	min_s = x+y;
	    	nn1 = nn2 = 0;
	    }
	    mov(list1, nn1, i, m1, x);
	    mov(list2, nn2, j, m2, y);
	}
	i = m1;
	j = m2;
	if (i >= n1 || j >= n2) break;
    }
    n1 = nn1; n2 = nn2;
}

static void reverse_list_opp(node *list, int n, int len)
{
    int i, j;
    for (i = 0, j= n-1; i<j; i++, j--) {
	node n = list[i];
	list[i] = list[j];
	list[j] = n;
	list[i].pos = opp(list[i].pos, len);
	list[j].pos = opp(list[j].pos, len);
    }
    if (i == j) list[i].pos = opp(list[i].pos, len);
}

int annot(node *list1, int n1, node *list2, int n2, int lower, int upper, int len1, int len2)
{
    int ii, jj;
    //mini_score(list1, n1);
    //mini_score(list2, n2);
    if (mis_penalty < 0) {
	mini_score(list1, n1);
	mini_score(list2, n2);
    } 
    if (n1 == 1 && n2 == 1 && list1[0].chrom != list2[0].chrom) {
        outputpair(list1, list2, 0, '*', len1, len2); // C**
        return 1;
    }
    if (n1+n2 == 1) {// D** class
	printf(",");
	if (n1 == 1) {list1[0].output(discont_input); printf("*");}
	else { printf("*");list2[0].output(discont_input);}
	return 1;
    }
    if (mis_penalty >= 0)  mini_score(list1, n1, list2, n2);
    int i = merge(list1, n1, list2, n2, lower, upper, ii, jj, '*', 1, len1, len2);
    if (ii >= 0) {
	outputpair(list1+ii, list2+jj, 0, '*', len1, len2);
	return 1;
    }
    i = merge(list2, n2, list1, n1, lower, upper, ii, jj, '*', 2, len1, len2); // same strand wrong orientation
    if (i > 0) return 1;
    reverse_list_opp(list2, n2, len2);
    //for (i = 0; i < n2; i++) list2[i].pos = opp(list2[i].pos, len2);
    int ii1, jj1;
    i = merge(list1, n1, list2, n2, lower, upper, ii1, jj1, '*', 3, len1, len2); // same orientation wrong strand
    if (i > 0) return 1;
    int ii2, jj2;
    i =  merge(list2, n2, list1, n1, lower, upper, ii2, jj2, '*', 4, len1, len2); // wrong oreintation and wrong strands
    if (i > 0) return 1;
    if (ii >=0) {
        //list2[ii].pos = opp(list2[i].pos, len2);  // need to change to original number
	reverse_list_opp(list2, n2, len2);
        outputpair(list2+ii,list1+jj, 2, '*', len2, len1);
        return 1;
    }
    if (ii1 >= 0) {
        outputpair(list1+ii1, list2+jj1, 3, '*', len1, len2);
        return 1;
    }
    if (ii2 >= 0) {
        outputpair(list2+ii2, list1+jj2, 4, '*', len1, len2);
        return 1;
    }
    return 0;
}

int compare_id(char *id, char *id1, int order)
{
    if (order) return 0;
    return compare_id(id, id1);
}

static int check_for_good_pair(char *line)
{
    char *p = strchr(line, ',');
    if (p == NULL) return 0;
    if (mis_penalty >= 0 || indelcut <= 0) return 1;
    do {
	p++;
	char *q = strchr(p, ',');
	char *qq = strchr(p, ':'); 
	if (qq == NULL) return 1;
	int len;
	p = q;
	if (qq && (!q || qq < q)) {
	    len = -1;
	    sscanf(qq+2, "%*d.%*d.%d", &len);
	    if (len<indelcut) continue;
	    qq = strchr(qq+2, ':');
	    if (qq && (!q || qq < q)) {
            	len = -1;
            	sscanf(qq+2, "%*d.%*d.%d", &len);
            	if (len<indelcut) continue;
		return 1; // both length long enough, a good pair is found
	    }
	} else return 1; //regular pairing, full length.
    } while (p);
    return 0;
}

int main(int argc, char *argv[])
{
    if (argc<5) {
	fprintf(stderr, "pairing(V2.2.4:05/26/2009) R3file F3file lower upper [M=mode][f=flen][r=rlen][Z=(L1,U1);(L2,U2)...][S=refseq] [O=0/1] [U=0/1] [T=##] [X=matrix] [F=maskingF] [R=maskingR] [G=0/1][I=##][D=##][i=##][d=##][P=0/1][q=##][V=tmp_res_file][z=maxrescue][a=TreshR][b=TreshF][B=##][E=##][u=snprate][H=iubmapping][A=VA][m=mis_p][L=min_len][h=histogramFile][x=dropscore][y=indelcutoff]\n");
	fprintf(stderr, "M=mode (default 0): 0 R3-F3 pairing;  1: F5-F3 pairing\n");
	fprintf(stderr, "O=0 sorted order; O=1 pairing by order(1st w 1st, 2nd w/ 2nd & so on)\n");
	fprintf(stderr, "U=0 all non-rescued pairs (default), U=1 only unique non-rescued pair\n");
	fprintf(stderr, "UR=0 all rescued pairs (default), UR=1 only unique rescued pair (default)\n");
	fprintf(stderr, "Uh=0 for histograms, include all pairs, Uh=1 only include unique pairs in histogram");
	fprintf(stderr, "T=## is the cutoff for rescue output\n P=0/1 i/o format 0 pos.mis, 1: chrom_pos.mis\nq= number of rescued alignments\n B,E specify the beginning and ending penal numbers;\n A=0/1 whether treat Valid adjacent as 1\n H=0/1 whether allow IUB mapping\n");
	exit(1);
    }
    char *file1 = argv[1];
    char *file2 = argv[2];
    int lower=atoi(argv[3]);
    int upper = atoi(argv[4]);
    pairingDistances.init(lower,upper);
    char *refseq = NULL;
    rescue_basic *re;
    int order = 0;
    int thres = 500;
    char *matr = NULL;
    char *maskingF = NULL;
    char *maskingR = NULL;
    int flen = 100000, rlen =100000;
    int gapalign = 2;
    int ano = 0;
    int numRes = 2000;
    char *tempfile = NULL;
    int maxrescue=100000000;
    int noaddition = 0;
    int TR = 100000, TF= 100000;
    float SNPrate = 0.0;


    int amballow = 0;
    int adjas1 = 0;
    char *qfile = NULL;
    string histFile;

    int i, ins = 0, del = 0, hi = 0, hd = 0;
    printf("#");
    for (i = 0; i <argc; i++) {
      printf("%s ", argv[i]);
    }
    printf("\n");

    pairingDistances.addHeaderArgLine(argc, argv);  //information added to header line of hist files
    pairingDistances.setReqUniqPairing(true);  //Histogram will only have unique pairings

    for (i=5; i<argc; ++i)
    {
    if (argv[i][2] == '=') {
    	if (argv[i][0] == 'U') {
    		if (argv[i][1] == 'R') {
    			flag_rescue_unique = atoi(argv[i]+3);
    			flag_rescue_was_set = 1;
    		}
    		if (argv[i][1] == 'h') {
    			int setting = atoi(argv[i]+3);
    			if (setting!=0)
    				pairingDistances.setReqUniqPairing(true);
    			else
    				pairingDistances.setReqUniqPairing(false);
    		}
    	}
    }
	if (argv[i][1] == '=')         // X=value
	{
	    if (argv[i][2] == 0)
		fprintf(stderr,"Missing value for %s\n",argv[i]);
	    else if (argv[i][0] == 'S') {
		refseq = argv[i]+2; //discont_input = 0;
	    } else if (argv[i][0] == 'O')
		order = atoi(argv[i]+2);
	    else if (argv[i][0] == 'U') {
	    	flag_unique = atoi(argv[i]+2);
	    	//if(flag_rescue_was_set == 0)
	    	//	flag_rescue_unique = flag_unique;
	    }
	    else if (argv[i][0] == 'T')
		thres = atoi(argv[i]+2);
	    else if (argv[i][0] =='M')
		mode = atoi(argv[i]+2);
            else if (argv[i][0] == 'X')
                matr = argv[i]+2;
            else if (argv[i][0] == 'F')
                maskingF = argv[i]+2;
            else if (argv[i][0] == 'R')
                maskingR = argv[i]+2;
	    else if (argv[i][0] == 'G')
		gapalign = atoi(argv[i]+2);
            else if (argv[i][0] == 'I')
                ins = atoi(argv[i]+2);
            else if (argv[i][0] == 'D')
                del = atoi(argv[i]+2);
	    else if (argv[i][0] == 'x')
                scoredrop = (int) (atof(argv[i]+2)*10);
            else if (argv[i][0] == 'd')
                hd = atoi(argv[i]+2);
            else if (argv[i][0] == 'i')
                hi = atoi(argv[i]+2);
            else if (argv[i][0] == 'f')
                flen = atoi(argv[i]+2);
            else if (argv[i][0] == 'r')
                rlen = atoi(argv[i]+2);
	    else if (argv[i][0] == 'P')
		discont_input = atoi(argv[i]+2);
	    else if (argv[i][0] == 'Z')
		getLR(argv[i]+2, lower, upper);
	    else if (argv[i][0] == 'Q')
		ano = atoi(argv[i]+2);
	    else if (argv[i][0] == 'q')
		numRes = atoi(argv[i]+2);
	    else if (argv[i][0] == 'a')
		TR = atoi(argv[i]+2);
            else if (argv[i][0] == 'b')
		TF = atoi(argv[i]+2);
	    else if (argv[i][0] == 'y')
		indelcut = atoi(argv[i]+2);
	    else if (argv[i][0] == 'V') {
		tempfile = argv[i]+2;
		if (tempfile[0] == '-') {
		    tempfile++;
		    noaddition = 1;
		}
	    } else if (argv[i][0] == 'z')
		maxrescue = atoi(argv[i]+2);
	    else if (argv[i][0] == 'L') 
		min_length = atoi(argv[i]+2);
	    else if (argv[i][0] == 'm') 
		mis_penalty = (int)(atof(argv[i]+2)*10.0-0.5);
	    else if (argv[i][0] == 'B') 
		startPanelNum = atoi(argv[i]+2);
	    else if (argv[i][0] == 'E')
                endPanelNum = atoi(argv[i]+2);
	    else if (argv[i][0] == 'u') {
                SNPrate = atof(argv[i]+2);
            } else if (argv[i][0] == 'H') {
                amballow = atoi(argv[i]+2);
            } else if (argv[i][0] == 'A') {
                adjas1 = atoi(argv[i]+2);
            } else if (argv[i][0] == 'h') {
	      vector<string> histVect; // can take in hist,bin_size
	      histFile.assign(argv[i]+2); //histogram file
	      histVect = split_string(histFile,",");
	      histFile.assign(*histVect.begin());
	      if(histVect.size()>1) { //if there's a bin_size specified, set it.
		pairingDistances.set_bin_size((size_t)atoi((histVect.begin()+1)->c_str()));
	      }
	    } else if (argv[i][0] == 's') {
	      usePanelList = true;
	      string inputStr(argv[i]+2);
	      vector<string> strVect;
	      strVect = split_string(inputStr,",");
	      //panelsToDo.clear();

	      vector<string>::const_iterator argPartItr = strVect.begin();
	      vector<int> panelVect;
	      for(;argPartItr!=strVect.end();argPartItr++)
		panelVect.push_back(atoi(argPartItr->c_str()));
	      sort(panelVect.begin(),panelVect.end()); //get into assending order

	      vector<int>::const_iterator panelVectItr = panelVect.begin();
	      for(;panelVectItr!=panelVect.end();panelVectItr++)
		panelsToDo.push(*panelVectItr);
	      if(panelsToDo.empty()) {
		fprintf(stderr,"ERROR, '%s' is not valid for list of panel ids (s arg).\n",
			inputStr.c_str());
		exit(1);
	      }
	      startPanelNum = panelsToDo.front();
	    }
	}
    }
    if (refseq == NULL) {
	rescue_no *rn = new rescue_no();
	re = rn;
    } else {
	FastaFile fafile(SEQTYPE_NT);
	FastaSeq faseq;

	if (!fafile.Open(refseq,"r")) {
	  //rescue_no *rn = new rescue_no();
	  // re = rn;
	  fprintf(stderr,"ERROR: Cannot find ref at %s.\n", refseq);
	  exit(1);
	} else {
	    rescue *r;
            if (gapalign == 1) {
            	r = new rescue();
            } else if (gapalign == 2  && (del >0 || ins > 0)) {
            	rescue_oneindel *r1 = new rescue_oneindel();
            	r1->setMaxIndel(ins, del, SSIZE, thres, hi, hd);
            	r = r1;
            } else {
            	rescue_nogap *r1 = new rescue_nogap();
            	r = r1;
            }
	    fafile.Close();
	    r->setNum(numRes);
	    r->ProcessFile(refseq);
	    r->setMatrix(matr);
	    r->Fmask(maskingF);
	    r->setLocal(mis_penalty, min_length);
	    r->Rmask(maskingR);
	    r->setmatching(adjas1, SNPrate, amballow, qfile, mis_penalty);
	    r->setMode(mode);
	    re = r;
	}
    }

    int scale = 10;
    if (qfile) scale = 100;
    re->setThreshold(thres, TR, TF, scale);
    re->set_unique_only(flag_rescue_unique);
    if  (flag_unique > 1) {
	max_out = flag_unique-1;
	flag_unique = 0;
    }
    if (flag_unique == 1) max_out = 1;
    node *list1 = new node[MAXSIZE], *list2= new node[MAXSIZE];
    if (list1 == NULL || list2 == NULL) fatal("Can not allocate memory\n");
    char id[1000], id1[1000], read1[1000], read2[1000], line3[100000];
    int n1, n2;
    genFile gf1, gf2;
    genFile *gF_fp1, *gF_fp2, gF_fp3;
    gF_fp1 = &gf1;
    if (strcmp(file1, file2)==0) gF_fp2 = gF_fp1;
    else {
   	gf2.setUsePipe();
    	gf2.open(file2,"r");
	gF_fp2 = &gf2;
    }
    gF_fp1->setUsePipe();
    gF_fp1->open(file1,"r");

    if(startPanelNum>0) {
      gF_fp1->seekToPanelId(startPanelNum);
      if (gF_fp2 != gF_fp1) gF_fp2->seekToPanelId(startPanelNum);
    }

    FILE *fp1 = gF_fp1->getFILE(); //ckopen(file1, "r");
    FILE *fp2 = NULL;
    if (gF_fp2 != gF_fp1) fp2 = gF_fp2->getFILE(); //ckopen(file2, "r");
    FILE *fp3 = NULL;

    line3[0] ='#';
    if (tempfile) { //this is the exclusion list
      gF_fp3.open(tempfile, "r");
      fp3 = gF_fp3.getFILE();  //ckopen(tempfile, "r");
    }
    readin(list1, n1, gF_fp1, id, read1, 1);
    readin(list2, n2, gF_fp2, id1, read2, 1);
    do {
	int a = compare_id(id, id1, order);
	if (a != 0 && fp2 == NULL) {
	    if (ano) {// E** class
		printf("%s", id);
		if (n1 == 1) {
		    printf(",");
		    if (id[strlen(id)-2] == 'F') printf("+");
		    list1[0].output(discont_input);
		    if (id[strlen(id)-2] != 'F') printf("+");
		}
                printf("\n");
	    }
	    node *s = list1;
	    list1 = list2;
	    list2 = s;
	    n1 = n2;
	    strcpy(id, id1);
	    strcpy(read1, read2);
	    if (!readin(list2, n2, gF_fp2, id1, read2, order)) break;
	    continue;
	}
	if (a < 0) {
	    if (ano) {
		printf("%s", id);
		if (n1 == 1) {
		    printf(",");
		    list1[0].output(discont_input);
		    printf("+");
		}
		printn();
	    }
	    if (!readin(list1, n1, gF_fp1, id, read1, order)) break;
	} else if (a > 0) {
            if (ano) {
                printf("%s", id1);
                if (n2 == 1) {
		    printf(",+");
                    list2[0].output(discont_input);
                }
                printn();
            }
	    if (!readin(list2, n2, gF_fp2, id1, read2, order)) break;
	} else {
	    if (fp3) {
		char *x = NULL;
		int b = 1;
		do {
		    if (line3[0] != '>') continue;
		    b = compare_id(id, line3, order);
		    if (b <= 0) break;
                } while (x= fgets(line3, sizeof line3, fp3));

		if (b == 0) {
		    if (check_for_good_pair(line3)) {
			if (noaddition == 0) print_line(line3);
            		if (!readin(list1, n1, gF_fp1, id, read1, order)) break;
            		if (!readin(list2, n2, gF_fp2, id1, read2, order)) break;
			continue;
		    }
		}
	    }

	    qsort((void *) list1, n1, sizeof(node), compare);
	    qsort((void *) list2, n2, sizeof(node), compare);
            print_line(id);
	    int ii, jj;
	    int len1 = lenofread(read1);
	    if (len1 > rlen) len1 = rlen;
	    int len2 = lenofread(read2);
	    if (len2 > flen) len2 = flen;

            if (mode ==1) reverse_list_opp(list2, n2, len2);
	    vector<int> nonResc_pairDists, r3Resc_pairDists, f3Resc_pairDists;
	    if (merge(list1, n1, list2, n2, lower, upper, ii, jj, '=', flag_unique, len1, len2,nonResc_pairDists)== 0) {
		if (mode == 1) reverse_list_opp(list2, n2, len2);
		if (mnum == 0)
		  re->pair_rescue(list1, n1, list2, n2, lower, upper, read1, read2, discont_input, len1, len2, maxrescue,r3Resc_pairDists, f3Resc_pairDists);
		else {
		    int i;
		    for (i = 0; i < mnum; i++)
		      re->pair_rescue(list1, n1, list2, n2, LL[i], UU[i], read1, read2, discont_input, len1, len2, maxrescue,r3Resc_pairDists, f3Resc_pairDists);
		}
		//if (ii >= 0) {
		//    outputpair(list1+ii, list2+jj, 0, '*', len1, len2);
		//} else {
		if (ano)
		    annot(list1, n1, list2, n2, lower, upper, len1, len2);
		//}
	    }
	    /*
	    std::cout << "{";
	    printVector(nonResc_pairDists,"+",std::cout);
	    std::cout << "_";
	    printVector(r3Resc_pairDists,"+",std::cout);
	    std::cout << "_";
	    printVector(f3Resc_pairDists,"+",std::cout);
	    std::cout << "}";
	    */

	    pairingDistances.add_3types_distances(r3Resc_pairDists,
						  f3Resc_pairDists,
						  nonResc_pairDists);

	    printn();
            if (!readin(list1, n1, gF_fp1, id, read1, order)) break;
            if (!readin(list2, n2, gF_fp2, id1, read2, order)) break;
	}
    } while(1);
    fclose(fp1);
    if (fp2) fclose(fp2);

    //Now make the histogram file if argument specified.
    //pairingDistances.bin_histograms(1);
    if(histFile.size()!=0) {
      ofstream histStrm;
      histStrm.open(histFile.c_str());
      pairingDistances.printHistograms(histStrm);
      histStrm.close();

      ofstream binStrm;
      histStrm.open((histFile+".binned").c_str());
      pairingDistances.bin_histograms();
      pairingDistances.printHistograms(histStrm);
      binStrm.close();
    }
    return 0;
}

