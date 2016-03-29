#include "map.h"

#include "zcompress.h"


#define ePCR_WDSIZE_DEFAULT   8
#define ePCR_WDSIZE_MIN       3
#define ePCR_WDSIZE_MAX       14
#define ePCR_THRESHOLD        180
#define MAXPROBELEN  1000
#define gpenalty 30
#define Aadiffer 'a'-'A'
static char adjclass[4][4] = {{0,1,2,3},{1,0,3,2},{2,3,0,1},{3,2,1,0}};
static char adjclassbig[128][128];
static char reorder[128];
static void sets2c(char s[][128], char a, char b, char c)
{
    s[a][b] = s[a+Aadiffer][b] = s[a][b+Aadiffer]= s[a+Aadiffer][b+Aadiffer] = c;
}

mapping_machine::mapping_machine ()
{
    init_scode(_scode);
    init_compl(_compl);
    
    mat = NULL;
    hash = NULL;
    psi = NULL;
    qfile = NULL;
    benum = 0; 
    mismatch = 0;
    do_classification = 0;
    outp = stdout;
    outputfile = NULL;
    //SetWordSize("11111111");
    SetThreshold(ePCR_THRESHOLD);
    amb_use_worst = 1;
    prefix_len = 0; 
    setgapsize(0);
    tempout = NULL;
    filter = 0;
    colorspace = 0;
    averageQV = 20.0; 
    dir = 1;
    both = 0;
    hlimit = 100000;
    selfscore = 0;
    reference = 0;
    firstr = 0; lastr = ((long long) (-1));
    int i, j;
   for (i = 0; i < 128; i++) {
	reorder[i] = '.';
        for (j = 0; j < 128; j++) {
	    seq2color[i][j] = 'N';
	    adjclassbig[i][j] = 4;
	    color2seq[i][j] = 'N';
	}
	seq2color[i]['.'] = 'X';
    }
    reorder['A'] = reorder['a'] = '0';
    reorder['C'] = reorder['c'] = '1';
    reorder['G'] = reorder['g'] = '2';
    reorder['T'] = reorder['t'] = '3';

    char order[] = "ACGT";
    for (i = 0; i < 4; i++) {
	for (j = 0; j < 4; j++) 
	    adjclassbig[order[i]][order[j]] = adjclass[i][j];
    }
    for (i = 0; i < 128; i++) {
	seq2color['.'][i] = 'X';
    }
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
    xdrop = 1000;
    uc = new unique_hits();
    int len = MAXPROBELEN;
    psi = new int*[len];
    //gap = new int[len];
    psi_qvalue = new int*[len+1];
    psi_qvalue[0] = new int[(len+1)*128];
    for (i = 1; i <= len; i++) psi_qvalue[i] = psi_qvalue[i-1]+128;
    for (i = 0; i < (len+1)*128; i++) psi_qvalue[0][i] = 1000;
}

mapping_machine_color_adja::mapping_machine_color_adja()
{
    int i;
    for (i = 0; i < 128; i++) amb_this[i] = i;
}	

mapping_machine_color_adj::mapping_machine_color_adj()
{
    score_snp = 100000;
}

void mapping_machine_color_adj::setSNPrate(double a, int amb)
{
    if (a >= 1.0) {
	score_adj = 10;
    } else if (a <= 0.000000001)  {
	score_adj = 100000000;
    } else {
	double x = -log(a)/log(10.0)*10.0/averageQV;
	x *=10;
	score_adj = (int) x;
    }
    if (amb == 1) {
	score_snp = 0;
	char a[]="KSWMYR"; //VBHD";
	char b[]="GCACTAACTA"; 
	char d[]="TGTACGGTCT";
	char c[]="ACGT";
	int n = strlen(a);
	int i,j;
	for (i = 0; i < n; i++) {
	    //sets2c(seq2color, a[i], 'N', seq2color[b[i]]['A']);
	    //sets2c(seq2color, 'N', a[i], seq2color['A'][b[i]]);
	    for (j = 0; j < 4; j++) {
		sets2c(seq2color,a[i], c[j], seq2color[b[i]][c[j]]);
		sets2c(seq2color, c[j],a[i],seq2color[c[j]][b[i]]); 
	    }
	    for (j = 0; j < n; j++) {
		sets2c(seq2color,a[i], a[j], seq2color[b[i]][b[j]]);
	    }
	}
	/*
	for (i = 0; i < 4; i++) {
            sets2c(seq2color, c[i], 'N', seq2color[c[i]]['A']);
            sets2c(seq2color, 'N', c[i], seq2color['A'][c[i]]);
	}
	*/
	/*
	for (i = 0; i < 10; i++) {
	    for (j = 0; j < 10; j++) {
		sets2c(seq2color,a[i], a[j], seq2color[b[i]][b[j]]);  
	    }
	}
	*/
	for (i = 0; i < 128; i++) amb_other[i] = 'N';
	for (i = 0; i < n; i++) amb_other[a[i]] = d[i];
	for (i = 0; i < n; i++) amb_this[a[i]] = b[i];
    }
}
 
void mapping_machine::setcolorcode(char *cs)
{
    char *p = cs;
    char a, b, c;
    char code[] = "ACGT";
    do {
	a = toupper(p[0]); b= toupper(p[1]); c= p[2];
	sets2c(seq2color, a, b, code[c-'0']);
	color2seq[a][c] = b;
	p = strchr(p, ',');
	if (p) p++; else break;
    } while (1);
}

mapping_machine_color::mapping_machine_color() 
{
    probe_first = NULL;
    colorspace = 1; 
}

//set matrix using default matrices
void mapping_machine::setMatrix() 
{

    int s =1;
    if (qfile) s=10;
    if (mat) delete mat;
    mat = new matrix(m_wsize, amb_use_worst, s);
    score_adj *= mat->Scale;
}

//set matrix reading from a file with filename given by mfile
void mapping_machine::setMatrix(char *mfile)
{
    int s =1;
    if (qfile) s=10;
    if (mat) delete mat;
    if (mfile) {
	mat = new matrix(mfile, m_wsize, amb_use_worst, s);
    } else mat = new matrix(m_wsize, amb_use_worst, s);
    score_adj *= mat->Scale;
}

//set matrix by reading a file already opened
void mapping_machine::setMatrix(FILE *mfp)
{
    if (mfp == NULL) {
	fatal("file is NULL");
    }
    mat = new matrix(mfp, m_wsize, amb_use_worst);
    score_adj *= mat->Scale;
}

// 
mapping_machine::~mapping_machine ()
{
	if (mat) delete mat;
	if (hash) delete hash;
	if (psi) delete psi;
	if (qfile) fclose(qfile);
}

//
void mapping_machine::SetThreshold (int t)
{
    threshold = t;
    if (qfile) threshold *=10;
}

// set wordsize m_wsize, m_asize, m_mask
void mapping_machine::SetWordSize (char *pat)
{
    if (pat == NULL) return;
    char *q = strchr(pat, ',');
    if (q) {
	if (hash == NULL) {HHHm *m = new HHHm(); hash = m;}
	hash->setpattern(pat, m_wsize, tailingzero);
    } else {
	if ( hash == NULL) {HHH1 *m = new HHH1();hash = m;}
	if (length == strlen(pat)) {
          int i, j = 0;
          inc[0] = 0;
          for (i = 1; i < strlen(pat); i++) {
            if (pat[i] != pat[i-1])  {
                if (pat[i] == '1') {
                    end[j] = i;
                    j++;
                } else {
                    if (j == 0) inc[j] = i;
                    else inc[j] = i-end[j-1];
                }
            }
          }
          if (pat[i-1] == '0') {
            end[j] = i;
            j++;
          }
          benum = j;
	}

	hash->setpattern(pat, m_wsize, tailingzero);
    } 
    if (mat) mat->setwindow(m_wsize);
}

static char code[] = "ACGT";

void mapping_machine::probe(char *p, int nh)
{
    int i;
    /* no translation in sequence space
    p[0] = color2seq[p[0]][p[1]-'0'];
    for (i = 1; i < length; i++) {
	p[i] = color2seq[p[i-1]][p[i+1]-'0'];
    }
    */
    if (nh == -1) return;
    if (prefix_len > 0) {
	reverse(p, prefix_len,probe_prefix);
	p+=prefix_len;
    }
    nh = probe_main(p, mat, nh);
    if (both) {
	reverse_seq(p, length, _compl);
	dir = -1;
	probe_main(p, mat, nh);
	dir = 1;
    }
    printn();
}
void mapping_machine_color::probe(char *p, int nh)
{
    int i;
    char a;
    if (nh >= hlimit ||  p[1] == 'X') goto rlabel;
    if (prefix_len > 0) {
	first_nuc = 'A';
	probe_first = mat->ma(0, '8', length);
    } else {
        first_nuc = color2seq[p[0]][p[1]];
        probe_first = mat->ma(0, first_nuc, length);
    }
    for (i = 0; p[i+2]; i++) {
	if (p[i+2] == 'X') goto rlabel;
        if (p[i+2]== '.' || p[i+2] == 'N') {
            p[i] = 'N';
        } else if (p[i+2] == '8') {
	    p[i] = '8';
	} else {
            p[i] = code[p[i+2]-'0'];
        }
    }
    p[i] = 0;
    if (prefix_len > 0) {
        reverse(p, prefix_len, probe_prefix);
        p+=prefix_len;
    }
    if (nh == -1) return;
    nh = probe_main(p, mat, nh);
    
    if (both) {
	int j;
	for (i = 0, j = length-1; i < j; i++, j--) {
	    char t = p[j];
	    p[j] = p[i];
	    p[i] = t;
	}
	first_nuc = _compl[first_nuc];
	if (prefix_len > 0) probe_first = mat->ma(0, '8', length);
	else probe_first = mat->ma(0,first_nuc, length);
	dir = -1;
	probe_main(p, mat, nh);
	dir = 1;
    }
rlabel:
    printn();
}

void mapping_machine_color_adj::probe(char *p,int nh) 
{
    int i;
    char a;
    if (nh >= hlimit || p[1] == 'X') goto exit_label;
    if (prefix_len > 0) {
        first_nuc = probe_second = 'A'; // when the probe_first is a wild card, both are not used.
        probe_first = mat->ma(0, '8', length);
    } else {
        a = first_nuc = color2seq[p[0]][p[1]];
        probe_first = mat->ma(0, first_nuc, length);
        probe_second = color2seq[a][p[2]];
    }
    for (i = 0; p[i+2]; i++) {
	if (p[i+2] == 'X') goto exit_label;
	int cla = 5;
	if (p[i+2] == '.' || p[i+2] == 'N') {
	    p[i] = 'N';
	} else {
            p[i] = code[p[i+2]-'0'];
            if (p[i+3]) cla = adjclass[p[i+2]-'0'][p[i+3]-'0'];
	    else cla = 0;
	}
	int j = i-prefix_len;
	if (j >= 0) probe_class[j] = cla;
    }
    p[i+1] = 0;
    if (prefix_len > 0) {
        reverse(p, prefix_len, probe_prefix);
        p+=prefix_len;
    }
    if (nh < 0) return;
    nh = probe_main(p, mat, nh);
    if (both) {
	int j;
	for (i = 0, j = length-1; i < j; i++, j--) {
	    char t = p[j];
	    p[j] = p[i];
	    p[i] = t;
	}
	for (i = 0, j = length-2; i < j; i++,j--) {
	    char t = probe_class[j];
	    probe_class[j] = probe_class[i];
	    probe_class[i] = t;
	}
	first_nuc = _compl[first_nuc];
        if (prefix_len > 0) probe_first = mat->ma(0, '8', length);
        else probe_first = mat->ma(0,first_nuc, length);
	probe_second = _compl[probe_second];
	dir = -1;
	probe_main(p, mat, nh);
	dir = 1;
    }
 exit_label:
    printn();
}


int mapping_machine::ReadprobeFile (const char *fname)
{
  //FILE *m_file = fopen(fname,"r");
        genFile m_file;
	m_file.setUsePipe();
	m_file.open(fname,"r");
	FILE *tpout = NULL;
	int adjlen = colorspace*2;
	if (tempout) {
	    if (strcmp(tempout,"stdin")==0) tpout =stdin; 
	    else tpout = fopen(tempout, "r");
	}
	if (m_file.getIsNull())
	{
		fprintf(stderr,"Unable to open file [%s]\n",fname);
		return 0;
	}
	long long countr = 0;
	char line[500000], p[1000];
	if (firstr > 0) {
	    while (m_file.gets(line,sizeof line)) {
		if (line[0] != '>') continue;
		countr++;
		if (countr == firstr) {
		    m_file.gets(line, sizeof line);
		    break;
		}
	    }
	}
	uc->setOffset(offset);	
	if (tpout == NULL) {
	    delete uc;
	    uc = new unique_hits();
	    while (m_file.gets(line,sizeof line))
	    {
		if (line[0] == '#') {
		    while (!(strchr(line, '\n'))) {
			if (m_file.gets(line, sizeof line) == NULL) break;
		    } 
		    continue;  // comment
		}
		if (line[0] != '>') {
		    if (line[0] == ' ' || line[0] == '\n') {
			fprintf(stderr, "Warning a line is skipped: empty line or start with space\n");
			continue ;
		    } 
		    fatal("reads file format is wrong, expecting > sign\n"); 
		}
		if (!(strchr(line, '\n'))) {
		    line[100] = 0;
		    fprintf(stderr, "read defline too long, check read file\n %s\n", line); 
		    exit(1);
		}
		//char seq[1000];
		int i = -1;
		p[0] = 0;
		sscanf(line, ">%s", probeid);
		char *q = strchr(probeid, ',');
		if (q) *q = 0;
		m_file.gets(line,sizeof line);
		sscanf(line, "%s", p);
		if (p[0] == 0) {
		    continue; 
		}
		real_length = strlen(p)-adjlen-prefix_len;
		if (shortout==0) { print_line(probeid); }
		if (real_length < length) {printn();continue;}
		p[0] = toupper(p[0]);
		probe(p, 0);
		if (countr == lastr) break;
		countr++;
	    }
	} else {
	    char line1[100000];
	    while (m_file.gets(line, sizeof line) && fgets(line1, sizeof line1, tpout))
	    {
		while (line[0] == '#') {
		    if (!m_file.gets(line, sizeof line)) {
			fclose(tpout);
			m_file.close();
			return 0;
		    } 
		}
		if (line[0] == ' ' || line[0] == '\n') {
		    fprintf(stderr, "warning, a line skipped start with space or end of line\n");
		}
		uc->renew();
		//char seq[1000];
		int i = -1;
		p[0] = 0;
		char *q = strchr(line1, '\n');
		while (q == NULL) {
		    print_line(line1);
		    if (fgets(line1, sizeof line1, tpout)  == NULL) return 0;
		    q = strchr(line1, '\n');
		}
		if (q) *q = 0;
		print_line(line1);
		int mini_mis = 1000;
		int nh = uc->addline(line1, mini_mis, prefix_len);
		sscanf(line, ">%s", probeid);
		q = strchr(probeid, ',');
		if (q) *q = 0;
		m_file.gets(line,sizeof line);
		sscanf(line, "%s", p);
		if (p[0] == 0) {
		    continue; 
		}
		if (mini_mis <= threshold-xdrop*10*mat->Scale) {printn(); continue;}
		real_length = strlen(p)-adjlen-prefix_len;
                if (real_length < length) {printn();continue;}
		if (filter) {
		    char *qq = strchr(line1, ',');
		    if (qq) {
			if (filter == 2) {printn(); continue;}
			int score;
			sscanf(qq+1, "%*d.%d", &score);
			if (score == 0) {printn(); continue;}
		    }
		}
		p[0] = toupper(p[0]);
		probe(p, nh);
		if (countr == lastr) break;
		countr++; 
	    }
	    if (tpout != stdin) fclose(tpout);
	}
	m_file.close();
     	if (outp != stdout) fclose(outp);
	return 1;
}

void mapping_machine::run_classification(const char *seq_data, const char *fname, const char *out_repeat, const char *outUniq, long long beg, long long end)
{
    do_classification = 1;
    int adjlen = colorspace*2;
    ProcessSeq("name", seq_data);
    FILE *frepeat = ckopen(out_repeat, "w");
    FILE *fother = NULL;
    if (outUniq) fother = ckopen(outUniq, "w");
    int count = 0;

    genFile m_file;
    m_file.open(fname,"r");
    fprintf(frepeat, "#max count for this 14me is %d\n", hash->max_count());
    char pp[1000], *p, line[10000];

        if (m_file.getIsNull())
        {
                fprintf(stderr,"Unable to open file [%s]\n",fname);
                exit(1);
        }

    while (m_file.gets(line,sizeof line))
    {
                if (line[0] == '#') {
                    while (!(strchr(line, '\n'))) {
                        if (m_file.gets(line, sizeof line) == NULL) break;
                    }
                    continue;  // comment
                }
		count++;
		if (count <= beg)  continue;
		if (count > end) break;

                //char seq[1000];
                pp[0] = 0;
                sscanf(line, ">%s", probeid);
                m_file.gets(line,sizeof line);
                sscanf(line, "%s", pp);
                if (pp[0] == 0) {
                    continue;
                }
                real_length = strlen(pp)-adjlen;
                if (real_length < length) {continue;}
	int i;
        for (i=0, p = pp+2; i<m_wsize; i++)
        {
            unsigned int j=(*p++)-'0';
            hash->addone(j);
        }
	int x, y;
	x = hash->findfirst();
	if (fother && x > threshold) {
	    fprintf(frepeat, ">%s\n%s", probeid, line);
	} else {
	    for (i = 0, p = pp+length+1; i < m_wsize; i++) {
		unsigned int j=(*p--)-'0';
		hash->addone(j);
	    }
	    y = hash->findfirst();
	    if (fother == NULL || y > threshold) {
		if (fother) 
            	    fprintf(frepeat, ">%s\n%s", probeid, line);
		else 
		    fprintf(frepeat, ">%s:%d\n%s", probeid, ((x>y)?x:y)+1, line);
 	    } else {
	    	fprintf(fother, ">%s\n%s", probeid, line);
	    }
        }
    }
}

int mapping_machine::header(unique_hits *uc, FILE *tpout, int shortout, int pos1, int filter)
{
    char line[1000000];
    uc->renew();
    if (tpout == NULL || fgets(line, sizeof line, tpout) == NULL) {
	if (shortout) { int outsize = fprintf(outp, "%d",pos1); check_printf;}
	return 1;
    }
    char *q = strchr(line, '\n');
    int mini_mis = 1000;
    uc->addline(line, mini_mis, prefix_len); 
    if (mini_mis <= threshold-xdrop*10*mat->Scale) return 0;
    if (q) *q = 0;
    print_line(line);
    if (filter) {
	char *qq = strchr(line, ',');
	if (qq) {
	    if (filter == 2) {printn(); return 0;}
	    int score;
	    sscanf(qq+1, "%*d.%d", &score);
	    if (score == 0) {printn(); return 0;}
	}
    }
    return 1;
}

int mapping_machine::buildprobe(const char *seq1)
{

    FILE *tpout = NULL;
    if (tempout) {
        if (strcmp(tempout,"stdin")==0) tpout =stdin;
        else tpout = fopen(tempout, "r");
    }
    long long len;
    const char *seq_data;
    char *seqcolor = NULL;
    if (seq1 && (len = strlen(seq1)) >= length)
    {
	seqc* seq;
	if (colorspace == 0) {
	    seq_data = seq1;
	} else {
	    seqcolor = new char[len];
	    int i;
	    for (i = 0; i < len-1; i++) {
		int a = _scode[seq2color[seq1[i]][seq1[i+1]]];
		if (a == AMBIG) seqcolor[i] = 'X';
		else seqcolor[i] = a+'0';
	    }
	    seqcolor[i] = 0;
	    seq_data = seqcolor;
	}
	int d = 1;
	char pr[1000];
	const char *p = seq_data;
	long long i;
	int j, pos, pos1,  N = 0;

	if (colorspace==0) {
	    for (pos1=0; pos1 <= len-length; ++pos1, p++)
	    {
		strncpy(pr, p, length);
		pr[length] = 0;
		if (header(uc,tpout, shortout, pos1, filter)) 
		    probe(pr, 0);
	    }
	} else {
	    for (pos1=0; pos1 < len-length; ++pos1, p++) {
		strncpy(pr+2, p, length);
		pr[length+2] = 0;
		pr[0] = toupper(seq1[pos1]);
		char a = _scode[seq2color[pr[0]][pr[0]]];
		if (a == AMBIG) pr[1] = 'X';
		else pr[1] = a+'0';
		if (header(uc, tpout, shortout, pos1, filter)) 
		    probe(pr, 0);
	    }
	    if (both)  {
	    for (i = 0; i < len-1; i++) {
		char a = _scode[seq2color[seq1[len-1-i]][seq1[len-i-2]]];
                if (a == AMBIG) seqcolor[i] = 'X';
                else seqcolor[i] = a+'0';
	    }
	    seqcolor[i] = 0;
	    for (pos1=0, p= seq_data; pos1 < len-length; ++pos1, p++) {
		strncpy(pr+2, p, length);
		pr[length+2] = 0;
		pr[0] = _compl[seq1[len-pos1-1]];
		char a = _scode[seq2color[pr[0]][pr[0]]];
                if (a == AMBIG) pr[1] = 'X';
                else pr[1] = a+'0';
		if (header(uc, tpout, shortout, pos1-len+1, filter)) 
		    probe(pr, 0);
	    } 
	    } //if both 
	}
	if (seqcolor) delete [] seqcolor;
	return 1;
    }
    return 0;
}

long long mapping_machine::init_hash(const char *seq_label, const char *seq1)
{
    long long len;
    if (reversed) both = 0;
    if (seq1 && (len = strlen(seq1)) >= m_wsize)
    {
	//seqc* seq;
	sname = strsave(seq_label); 
	seqseq = strsave(seq1, len);
	long long j;
	for (j = 0; j < len; j++) {
	    if (_compl[seqseq[j]] == 'X') seqseq[j] = 'N';
	}
	if (reversed) reverse_seq(seqseq, len, _compl);
	//fprintf(stderr, "length is %lld\n", len); 
	if (!do_classification) hash->setseq(len);
        totallen = len;
        if (colorspace == 0) {
            seqcolor = NULL;
        } else {
            seqcolor = new char[len];
            int i;
            char *s = seqcolor, *t = seqseq;
            while (s < seqcolor+len-1) {
                *s = seq2color[*t][*(t+1)];
                *t = toupper(*t);
                s++; t++;
                //seqcolor[i] = seq2color[seqseq[i]][seqseq[i+1]];
            }
            *s = 0;
            //seqcolor[i] = 0;
            //translate to color space;
        }

	return len;
    } 
    return 0;
}

int mapping_machine::build_hashing(long long beg, long long end)
{
    const char *seq_data;

	if (colorspace == 0) {
	    seq_data = seqseq;
	    seqcolor = NULL;
	} else {
	    seq_data = seqcolor;
 	}	
	int d = 1;
	    const char *p = seq_data+beg;
	    int  i;
	    int  N = 0;
	    unsigned int pos1;
	    //hashitem->renew();
	    for (i=0; i<m_wsize-1; i++)
	    {
		unsigned int j=_scode[*p++];
		hash->addone(j);
	    }
	    
	    for (pos1= (unsigned int) beg; *p && pos1 <= end; ++pos1)
	    {
		//// If N > 0 it means there is an N within the the
		//   first m_wsize characters and we just ignore it.
		hash->addone(_scode[*p++]);
		hash->insert(pos1);
	    }
	return 1;
}

int mapping_machine::ProcessSeq (const char *seq_label, const char *seq1)
{
     long long len;
    if ((len =init_hash(seq_label, seq1))== 0) return 0;
    return build_hashing(0, len);
}

// Display the alignment
void mapping_machine::display_align_core(const char *seq_data, int pos, int pos_e, int left, int right)
{
    char line[1000];
    char line1[1000];
    strncpy(line, seq_data+pos, pos_e-pos+1);
    strncpy(line1, probe_seq+left, length-left-right);
    line[pos_e-pos+1] = '\0';
    line1[length-left-right] = '\0';
    
    if (length-left-right == pos_e-pos+1) {
	int i;
	for (i = 0; i < pos_e-pos+1; i++) {
	    if (line1[i] == line[i]) line[i] = '.';
	}
    }
    int outsize = fprintf(outp, "%s\n%s\n", line1, line);
    check_printf;
}

void mapping_machine::display_align(const char *seq_data,int pos, int pos_e)
{
    if (gapsize > 0) {
	return;
    } else { 
	display_align_core(seq_data, pos,pos_e, 0,0);
    }
}
static void build_seq(char *det, const char *src, int b, int e, int r)
{
    int i;
    char *a = det;
    if (r) {
	for (i = e; i >= b; i--) {
            *a = reorder[src[i]];
            a++;
        }
	*a = 0;
	return;
    }
    for (i = b; i <= e; i++) {
	*a = reorder[src[i]];
	a++;
    }
    *a = 0;
}
// This funtion output each hit to  probe.

static long long position_convert(long long pos, long long pos_e, int reversed, int dir, long long slen, long long off)
{
        if (reversed) {
	    return pos-slen+1-off;
        } else if (dir == -1) {
	    return -pos_e-1-off;
        } else {
	    return pos+off;
        }
}

//seq1 is probe, seq2 is reference
int mapping_machine::left_local_align(const char *seq1, const char *seq2, int match, int mismatch, int &len, int dir, int lstate)
{
    len = 0;
    int max_score = 0;
    const char *s = seq1; 
    int cur_score = 0;
    while (*s && *seq2 && *seq2 != 'X') { 
	if (*s == *seq2) {
	    cur_score += match;
	    if (cur_score > max_score) {
		max_score = cur_score;
		len = s-seq1+1;
	    }
	} else cur_score += mismatch;
	s++; seq2 += dir;
    }
    return max_score; 
}

int mapping_machine::left_local_align(const char *seq1, const char *seq2, int match, int mismatch, int &len, int dir, int s, int t)
{
    return left_local_align(seq1, seq2, match, mismatch, len, dir, adjclassbig[s][t]);
}

int mapping_machine_color_adj::left_local_align(const char *seq1, const char *seq2, int match, int mismatch, int &len, int dir, int lstate)
{
    len = 0;
    int max_score = 0;
    const char *s = seq1;
    int cur_score = 0;
    while (*s && *seq2 && *seq2 != 'X') {
        if (*s == *seq2) {
            cur_score += match;
	    lstate = 0;
            if (cur_score > max_score) {
                max_score = cur_score;
                len = s-seq1+1;
            }
        } else {
	    int state = adjclassbig[*s][*seq2];
	    if (lstate == state) {
		cur_score += match;
		if (cur_score > max_score) {
                    max_score = cur_score;
                    len = s-seq1+1;
                }
		lstate = 0; // can not start a VA again.
	    } else {
		lstate = state;  
		cur_score += mismatch;
	    }
	}
        s++; seq2 += dir;
    }
    return max_score;
}

long long mapping_machine::ReportHit(const char *seq_label, const char *seq_data, long long  pos, long long pos_e, int sum1, long long slen)
{
    long long x = pos+offset;
    if (reversed) {
        x= pos-slen+1-offset;
    } else if (dir == -1) {
        if (colorspace == 0) {
     	    x=-pos_e-offset;
        } else {
            x=-pos_e-1-offset;
        }
    }

    if (sum1 == 0 && noperfect) return 0;
    int outsize;
    int len, lscore, lscore1;
    if (mismatch != 0) {
	if (dir==1) 
		lscore = left_local_align(probe_seq+pos_e-pos+1, seq_data+pos_e+1, 10, mismatch, len, dir, adjclassbig[probe_seq[pos_e-pos]][seq_data[pos_e]]);
	else  lscore = left_local_align(probe_seq+pos_e-pos+1, seq_data+pos-1, 10, mismatch, len, dir, adjclassbig[probe_seq[0]][seq_data[pos]]);
	lscore += ((pos_e-pos+1)*10+sum1/10*(mismatch-10))*mat->Scale;
	if (dir == 1) pos_e += len;
	else pos -= len;
	if (prefix_len > 0) {
	    if (dir == 1) {
		lscore1 = left_local_align(probe_prefix, seq_data+pos-1, 10, mismatch, len, -1, adjclassbig[probe_seq[0]][seq_data[pos]]); 
		pos -= len;
	    } else {
		lscore1 = left_local_align(probe_prefix, seq_data+pos_e+1, 10, mismatch, len, 1, adjclassbig[probe_seq[length-1]][seq_data[pos_e]]);
		pos_e+= len;
	    }
	    lscore += lscore1;
	}
	lscore = ((pos_e-pos+1)*10-lscore)/(10-mismatch); 
    }
    if (display_alig) {
	outsize = fprintf(outp, "=%s %s %lld %lld %d\n", probeid, seq_label, pos, pos_e, sum1);
	check_printf;
	display_align(seq_data, pos, pos_e);
    } else {
	if (reversed) {
		outsize = fprintf(outp, ",%lld.%d",  pos-slen+1-offset, sum1);
	} else if (dir == -1) {
	    if (colorspace == 0) {
		outsize = fprintf(outp, ",%lld.%d", -pos_e-offset, /*-pos,*/ sum1);
	    } else {
		outsize = fprintf(outp,",%lld.%d",-pos_e-1-offset, /*pos_e-slen,*/ sum1);
	    }
	} else {
	    outsize = fprintf(outp, ",%lld.%d", pos+offset, /*pos_e,*/ sum1);
	}
	check_printf;
    }
    if (mismatch != 0) {
	if (prefix_len > 0)  outsize = fprintf(outp, ":(%lld.%d.%d)", pos_e-pos+1,lscore, prefix_len-len);
	else outsize = fprintf(outp, ":(%lld.%d.0)", pos_e-pos+1, lscore);
	check_printf;
    }
    if (reference) {
	char tmp[1000];
        if (colorspace == 1) {
            if (dir == 1)
		build_seq(tmp, seq_data, pos, pos+real_length-1, reversed);
	    else 
		build_seq(tmp, seq_data, pos_e-real_length+1, pos_e, 0);
	    if (!(reversed || dir == -1)){
		outsize = fprintf(outp, ":%c%s",seqseq[pos], tmp);
	    } else {
		outsize = fprintf(outp, ":%s%c", tmp, (dir == -1) ?  seqseq[pos_e+1] : _compl[seqseq[pos]]);  
	    }
        } else {
	    strncpy(tmp, seq_data+pos, pos_e-pos+1);
	    tmp[pos_e-pos+1] = 0;
            outsize = fprintf(outp,":%s", tmp);
	}
	check_printf;
    }
    return x;
}



int mapping_machine::match(long long seed, int index, int wscore, int threshold, int &sum1, long long slen) 
{
    if (dir == 1 && !uc->is_new(seed, seed+length-1, reversed, dir, slen)) return 0;
    if (dir == -1 && !uc->is_new(seed, seed+length, reversed, dir, slen)) return 0;
    const char *seq_label = sname;
    const char *seq = seqseq;
    int lindel, rindel;
    int limit;
    int flag;
    limit = threshold;
    sum1 = 0;
    if (gapsize == 0) {
	flag = ungappedmatch(seq, seed, limit, sum1, lindel, rindel, 0);
	if (flag > 0){
	    int x = ReportHit(seq_label, seq, seed-lindel, seed-1+length-rindel, sum1, slen);
	    uc->addhits(x);
	}
    } else {
	int score = gapmatch(seq, seed, limit, gapsize, 0);
	if (score >= 0){
	    flag = 1;
	    sum1 = score;
	    ReportHit(seq_label, seq, seed+1, seed+length, score, slen);
	    int i;
	    /*
	    for (i = seed-gapsize; i <= seed+gapsize; i++) {
		p->newhit_test(i+seq_count);
	    }
	    */
	    if (display_alig) {
		gapmatch(seq, seed, limit, gapsize, 1);
	    }
	} else flag = -1;
    }
    return flag;
}

static int cross_dot(const char *seq, long long start, int len, long long slen)
{
    long long end = start+len;
    if (end > slen) return 1;
    if (start < 0) return 1;
    const char *s = seq+start; 
    while (s < seq+end) {
	if (*s == 'X') return 1;
	s++;
    }
    return 0;
}

int mapping_machine_color::match(long long seed, int index, int wscore, int threshold, int &sum1, long long slen)
{
    if (!uc->is_new(seed, seed+length-1, reversed, dir, slen)) return 0;
    const char *seq_label = sname;
    const char *seq = seqcolor;
    const char *seq2 = seqseq;
    //dir = hn->seq->dir;
    int lindel = 0, rindel = 0;
    int limit;
    int flag;
    limit = threshold;
    // color space no gap;
    {
	sum1 = selfscore;
    	int **p;
	const char *s;
	if (seed < 0) return -1;
	p = psi;
	s = seq+seed;
	if (benum == 0) {
	    while (p < psi+length) {
	        sum1 += (*p)[*s];
                if (sum1 > threshold) return -1;
                p++; s++;
            }
	} else {

        int i;
        for (i = 0; i < benum; i++) {
            p += inc[i];
            s += inc[i];
            while (p < psi+end[i]) {
                sum1 += (*p)[*s];
                if (sum1 > threshold) return -1;
                p++; s++;
            }
        }

	}
	/*
        if (dir == 1)
            sum1 += probe_first[seq2[seed]];
        else sum1 += probe_first[seq2[seed+length]];
	if (sum1 > threshold ) return -1;
	*/ 
	//flag = ungappedmatch(seq2, seed, limit, sum1, lindel, rindel, 0);
	//if (flag > 0){
	    if (real_length>length) {
		int i = real_length-length; 
		if (dir == 1) {
		    if (cross_dot(seq,seed+length, i, slen)) return -1;
		} else {
		    if (cross_dot(seq,seed-i, i, slen)) return -1;
		}
	    }

        if (dir == 1)
            sum1 += probe_first[seq2[seed]];
        else sum1 += probe_first[seq2[seed+length]];
        if (sum1 > threshold ) return -1;

	int x= ReportHit(seq_label, seq, seed-lindel, seed+length-1-rindel, sum1, slen);
	uc->addhits(x);
	//}
    }
    return 1;
}

int mapping_machine_color_adja::match(long long seed, int index, int wscore, int threshold, int &sum1, long long slen)
{
    if (!uc->is_new(seed, seed+length-1, reversed, dir, slen)) return 0;
    const char *seq_label = sname;
    const char *seq = seqcolor;
    const char *seq2 =seqseq;
    int lindel, rindel;
    int limit;
    int flag;
    limit = threshold;
    sum1 = 0;
        if (dir == 1)
            sum1 = probe_first[amb_this[seq2[seed]]];
        else sum1 = probe_first[amb_this[seq2[seed+length]]];

        flag = ungappedmatch(seq, seed, limit, sum1, lindel, rindel, 0);
        if (flag > 0){
            if (real_length>length) {
                int i = real_length-length;
                if (dir == 1) {
                    if (cross_dot(seq,seed+length, i, slen)) return -1;
                } else {
                    if (cross_dot(seq,seed-i, i, slen)) return -1;
                }
            }
            int x= ReportHit(seq_label, seq, seed-lindel, seed-1+length-rindel, sum1, slen);
	    uc->addhits(x);
	
        }
 	return flag;
}


// Ungapped alignment for a probe. In the ungapped case, if not local (global),
// one just need to walk along the diagonal in the alignment graph for a fiexd
// length (the length of the probe). In the local alignment case, one need to
// keep a running optimum so in the end we can truncate alignment to get the
// best local alignment.
int mapping_machine::ungappedmatch(const char *seq, long long seed, int threshold, int &sum1, int &ll, int &rr, char islocal) 
{
    int **p;
    const char *s;
    int sum = sum1+selfscore, maxi = -10000;
    if (seed < 0) return -1;
    p = psi;
    s = seq+seed;
    if (benum == 0) {
    	while (p < psi+length) { 
	    sum += (*p)[*s];
	    if (sum > threshold) return -1;
	    p++; s++;
       	}
    } else {
	int i;
	for (i = 0; i < benum; i++) {
	    p += inc[i];
	    s += inc[i];
	    while (p < psi+end[i]) {
            	sum += (*p)[*s];
            	if (sum > threshold) return -1;
            	p++; s++;
            }
	}
    }

    //if (sum > threshold) return -1;
    sum1 = sum;
    ll = rr = 0;
    return 1;
    
}

 
int mapping_machine_color_adja::ungappedmatch(const char *seq, long long seed, int threshold, int &sum1, int &ll, int &rr, char islocal)
{
    int **p;
    const char *s;
    int sum = sum1, maxi = -10000;
    if (seed < 0) return -1;
    p = psi;
    s = seq+seed;
    int len1 = length;
    if (sum > 0) {
        if (sum > threshold ) return -1;
        if (dir == 1) {
                s++;p++;
        } else {
                len1--;
        }
    }
   sum += selfscore;
    while (p < psi+len1) {
        int a = (*p)[*s];
        if (a > 0) {
            sum+= a;
            if (sum > threshold) return -1;
            p++;s++;
        }
        p++; s++;
    }
    //if (sum > threshold) return -1;
    sum1 = sum;
    ll = rr = 0;
    return 1;
}

int mapping_machine_color_adj::ungappedmatch(const char *seq, long long seed, int threshold, int &sum1, int &ll, int &rr, char islocal)
{
    int **p;
    const char *s;
    int sum = sum1, maxi = -10000;
    if (seed < 0) return -1;
    p = psi;
    s = seq+seed;
    int len1 = length;
    if (sum > 0) {
	if (sum > threshold ) return -1;
	if (dir == 1) {
	    if (probe_second == seqseq[seed+1]) {
		if (sum > score_snp && first_nuc==amb_other[seqseq[seed]]) {
		    sum = score_snp;
		} else {
		    sum+= (*p)[*s];
		    if (sum > score_adj) sum = score_adj; 
		}
                s++;p++;
	    } 
	} else {
	    if (probe_second == seqseq[seed+length-1]) {
		len1--;
                if (sum > score_snp && first_nuc ==amb_other[seqseq[seed+length]]) {
                   sum = score_snp;
                } else {
		    sum += psi[len1][seqseq[seed+len1]];
		    if (sum > score_adj) sum = score_adj;
		} 
	    }
	}
    }
    sum += selfscore;
    while (p < psi+len1-1) {
	int a = (*p)[*s];
	if (a > 0) {
	    if (probe_class[p-psi] == adjclassbig[*s][*(s+1)]) { 
	    	p++;s++;
		int d = p-psi;
		if (a > score_snp && probe_seq[d] == seq2color[amb_other[seqseq[seed+d]]][seqseq[seed+d+1]]){
		    a = score_snp;
		} else {
		    a+= (*p)[*s];
		    if (a > score_adj) a = score_adj;
		}
	    }
	    sum += a;
	    if (sum > threshold) return -1;
	}
        p++; s++;
    }
    if (p < psi+len1) {
	sum+=(*p)[*s];
	if (sum > threshold) return -1; 
    }	
    sum1 = sum;
    ll = rr = 0;
    return 1;
}

// Initialize probe. Psi stores the profile for each base (like in
// profile search). Each positionof the probe has an array telling
// use what score will be if it is matched to different nucleotides.
void mapping_machine::probe_init(int len) 
{
    //psi = new int*[len];
    length = len;
    if (len > MAXPROBELEN) fatal("probe too long\n");
/*	
    psi_qvalue = new int*[len+1];
    psi_qvalue[0] = new int[(len+1)*128];
    int i;
    for (i = 1; i <= len; i++) psi_qvalue[i] = psi_qvalue[i-1]+128; 
    for (i = 0; i < (len+1)*128; i++) psi_qvalue[0][i] = 1000;
*/
}
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
/*
    0.75,
    0.75,
0.6309573444801930,
0.5011872336272720,
0.3981071705534970,
0.3162277660168380,
0.2511886431509580,
0.1995262314968880,
0.1584893192461110,
0.1258925411794170,
0.1000000000000000,
0.0794328234724281,
0.0630957344480193,
0.0501187233627272,
0.0398107170553497,
0.0316227766016838,
0.0251188643150958,
0.0199526231496888,
0.0158489319246111,
0.0125892541179417,
0.0100000000000000,
0.0079432823472428,
0.0063095734448019
*/
};
 
void mapping_machine::modify(int x, int i, int index)
{
    char a[]="ACGTacgt", c = probe_seq[i-1];
    int j, *old = psi[i-1], sum;
    double r, p;
    int *nw  = psi_qvalue[index]; 
    if (i == 0) {
	old = probe_first;
	c = first_nuc;
	probe_first = nw;
    } else {
	psi[i-1] = nw;
    }
    sum = old['A']+old['C']+old['G']+old['T'];
    if (x == 0) x=10;
    p = q2p[x]/averageQV;//  prob of wrong call
    r = ((double)(x+5))/averageQV;
    int s = (int) (p*sum/3);
    for (i = 0; i < 8; i++) {
	nw[a[i]] = (int) (r*old[a[i]])-s;
    }
    nw[c] = 0; // match value
    selfscore += s;
    nw['N'] = old['N'];
    nw['X'] = old['X'];
}


char line_map[10000];

int mapping_machine::quality_value(char *qline)
{
    int i, index_tempmem = 0;
    char *q = qline-1;
    selfscore = 0;
        for (i = 0; i <= length; i++) {
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

void mapping_machine::setPsi(const char *seq)
{
    
    strcpy(probe_seq, seq);
    int i;
    if (dir == 1)
        for (i = 0; i < length; i++) {
            psi[i] = mat->ma(colorspace+i, probe_seq[i], length);
        }
    else
        for (i = 0; i < length; i++) {
            psi[i] = mat->ma(length-i-(1-colorspace), probe_seq[i], length);
        }
}

void mapping_machine_color::setPsiMask(const char *seq, int *allzero, const char *mask)
{
    strcpy(probe_seq, seq);
    int i;
    for (i = 0; i < length; i++) {
	if (mask[i+1] == '1') {
            psi[i] = mat->ma(i+1, probe_seq[i], length);
        } else psi[i] = allzero;
    }
	if (mask[0] == '0') probe_first = allzero; 
}

void mapping_machine::getQ()
{

    if (qfile) {
        selfscore = 0;
        if (!both || dir == 1) {
            line_map[0] = 0;
            while (fgets(line_map, sizeof line_map, qfile)) {
                if (line_map[0] != '>' && line_map[0] != '#') break;
            }
        }
        if (line_map[0] != 0) {
            quality_value(line_map);
            /*
            int x, index_tempmem = 0;
            char *q = line_map-1;
            for (i = 0; i <= length; i++) {
                if (q == NULL) break;
                q++;
                x = atoi(q);
                q = strchr(q, ' ');
                if (x > 20) {x=20;}
                modify(x, i, index_tempmem);
                index_tempmem+
	    }
	   */
	}
    }
}

int mapping_machine::probe_main(char *seq, matrix *pm, int nh) 
{
    int i;
    if (nh >= hlimit) return nh;
    int last_map=-1;
    setPsi(seq);
    if (qfile) {
	selfscore = 0;
	if (!both || dir == 1) {
	    line_map[0] = 0;
	    while (fgets(line_map, sizeof line_map, qfile)) {
		if (line_map[0] != '>' && line_map[0] != '#') break;
	    }
	}
	if (line_map[0] != 0) {
	    quality_value(line_map);
	    /*
	    int x, index_tempmem = 0;
	    char *q = line_map-1;
	    for (i = 0; i <= length; i++) {
		if (q == NULL) break;
		q++;
		x = atoi(q);
		q = strchr(q, ' ');
		if (x > 20) {x=20;}
		modify(x, i, index_tempmem);
		index_tempmem++;
	    }
	    */
	} 
    }
    /* new code , commented out for now
    
    if (has_qvalue) {
	if (dir == 1) {
	    for (i = 0; i < length; i++) {
		if (n2s[i]<qthresh) continue;
		int *second = pm->ma(colorspace+i, probe_seq_sec[i], length);
		average(psi_save[i], psi[i], second, n2s[i]);
	    }
	} else {
            for (i = 0; i < length; i++) {
                if (n2s[i]<qthresh) continue;
                int *second = pm->ma(length-i-(1-colorspace), probe_seq_sec[i], length);
                average(psi_save[i], psi[i], second, n2s[i]);
            }
	}
    }
    */
    /*
    for (i = 0; i < p->length; i++) {
	gap_p[i] = gap[i]-psi[i][p->probe_seq[i]];
    }
    */
    if (length >= m_wsize) {
	char *p = probe_seq;
	//hashitem->renew();
	for (i=0; i<m_wsize-1; i++)
	{
	    unsigned int j=_scode[*p++];
	    hash->addone(j);
	}
	int pos1 = 0;
	if (tailingzero) {
	  int N = 0;
	  for (; pos1 <= length-m_wsize; ++pos1)
	  {
	    hash->addone(_scode[*p++]);
	    //int h = hashitem->hash();
	    if (N > 0) {
		N--;
		continue;
	    } else N=tailingzero;
	    {
		int sum, sum1;
		long long hn;
		int j = 0;
		for (hn = hash->findfirst(); hn >= 0; hn = hash->findnext()) {
		    long long pos = hn-pos1;
		    //if (pos != last_map) {
			if (match(pos, pos1, 0, threshold, sum1, totallen)>0) {
			    //last_map = pos;
			    nh++; j= 0;
			    if (nh >= hlimit) return nh;
			} else {
			    j++;
			    if (num_ext_allowed && j >= num_ext_allowed) return nh;
			}
		
		    //}
		}
	    }
	  }
	  return nh;
	}
	

	for (; pos1 <= length-m_wsize; ++pos1)
	{
	    hash->addone(_scode[*p++]);
	    /*
	    int h = hashitem->hash();
	    if (h >= 0)
	    */
	    {
		int sum, sum1;
		long long hn;
		for (hn = hash->findfirst(); hn >= 0; hn = hash->findnext()) {
		    long long pos = hn-pos1;
		    //if (pos != last_map) {
			int j = 0;
			if (match(pos, pos1, 0, threshold, sum1, totallen)>0) {
			    nh++;j = 0;
			    if (nh >= hlimit) break;	
                            //last_map = pos;
			}else {
                            j++;
                            if (num_ext_allowed && j >= num_ext_allowed) return
nh;
                        }

		    //}
		}
	    }
	}
    }
    return nh;
}

void makehash_sdiscont::addone(int a)
{
	addbase(a);
        h1 <<=2;
        h1 &=mask1;
        h1 |= (h2 >> shift2);
        h2 <<=2;
        h2 &= mask2;
	h2 |= (h3 >>shift3);
	h3 <<= 2;
	h3 &=mask3;
        if (a == AMBIG) a = 0;
        h3 |= a;
        hashv = (int) ((h1 << shift1) | (h3));
}

void makehash_sdiscont::set_param(int zero1, int one1, int zero2, int one2)
{
    shift2 = 2*zero2;
    mask1 = (1 << (2*one1))-1;
    mask2 = (((long long)1) <<(shift2))-1;
    shift2 -= 2;
    shift1 = one2*2;
    mask3 = (1 << shift1)-1;
    shift3 = shift1-2;
    wd = one1+one2+zero1+zero2;
}
makehash_sdiscont::makehash_sdiscont()
{
    h1 = h2 = h3=0;
    mask1 = mask2 = mask3 = 0;
    shift1 = shift2 = shift3 = 0;
}

void makehash_simple_discont::addone(int a) 
{
	addbase(a);
	h1 <<=2;
	h1 &=mask1;
	h1 |= (h2 >> shift2);
	h2 <<=2;
	h2 &= mask2;
	if (a == AMBIG) a = 0;
	h2 |= a;
	hashv = (int) ((h1 << shift1) | (mask & h2));
}
    
void makehash_simple_discont::set_param(int zero1, int one1, int zero2, int one2) 
{
    shift2 = 2*zero2+2*one2;
    mask1 = (1 << (2*one1))-1;
    mask2 = (((long long)1) <<(shift2))-1;
    shift2 -= 2;
    shift1 = one2*2;
    mask = (1 << shift1)-1;
    wd = one1+one2+zero1+zero2;
} 

makehash_simple_discont::makehash_simple_discont() 
{
    h1 = h2 = 0;
    mask1 = mask2 = mask = 0;
    shift1 = shift2 = 0;
}

makehash_3one_discont::makehash_3one_discont() 
{
    h1 = h2 = 0;
    mask1 = mask2 = mask3 = mask2n= mask3n = 0;
    shift1 = shift2 = shift3 = shift2n=0;
}

void makehash_3one_discont::addone(int a) 
{
	addbase(a);
	h1 <<=2;
	h1 &=mask1;
	h1 |= (h2 >> shift2);
	h2 <<=2;
	h2 &= mask2;
	h2 |= (h3 >> shift3);
	h3 <<= 2;
	h3 &= mask3;
	if (a == AMBIG) a = 0;
	h3 |= a;
	hashv = (int) ( (h1 << shift1) | (mask2n & (h2 << shift2n)) | (mask3n & h3));
}


    
void makehash_3one_discont::set_param(int zero1, int one1, int zero2, int one2, int zero3, int one3) 
{
    shift2 = 2*zero2+2*one2;
    shift3 = 2*zero3+2*one3;
    mask1 = (1 << (2*one1))-1;
    mask2 = (((long long) 1) <<(shift2))-1;
    mask3 = (((long long) 1) <<(shift3))-1;
    shift2 -= 2;
    shift3 -= 2;
    shift1 = one2*2+one3*2;
    mask2n = (1 << shift1)-1;
    shift2n = one3*2;
    mask3n = (1 << shift2n)-1;
    wd = one1+one2+zero1+zero2+zero3+one3;
} 

void makehash_gen_discont::addone(int a)
{
    addbase(a);
    int i;
    hashv = 0;
    for (i = 0; i < n_ones-1; i++) {
	h[i] <<=2;
	h[i] &= mask[i];
	h[i] |= (h[i+1] >> shift[i+1]);
	hashv |= (int) ((h[i] & maskn[i]) << shiftn[i]);
    }
    h[i] <<=2;
    h[i] &= mask[i];
    if (a == AMBIG) a = 0;
    h[i] |= a;
    hashv |= (int) (h[i] & maskn[i]);
}

makehash_gen_discont::makehash_gen_discont()
{
    memset(shift, 0, 100*sizeof(int));
    memset(h, 0, 100*sizeof(long long));
    memset(shiftn, 0, 100*sizeof(int));
    memset(mask, 0, 100*sizeof(long long));
    memset(maskn, 0, 100*sizeof(long long));
    n_ones = wd=0;
}

void makehash_gen_discont::set_param(char *pat)
{
    char *p = pat;
    int zero, one, i;
    wd = n_ones = 0;
    do {
	zero = one = 0;
	for (; *p; p++) {
	    if (*p != '0') break;
	    zero++;
	}
	for (; *p; p++) {
	    if (*p != '1') break;
	    one++;
	}
	if (zero == 0 && one == 0) break;
	if (n_ones == 0) {
	    mask[0] = maskn[0] = (((long long) 1) << (one*2))-1;
	    shift[0] = 0;
	} else {
	    while (one+zero > 30) {
		int s = 30;
		if (s > zero-3) s = zero-3; 
		shift[n_ones] = 2*s;
		mask[n_ones] = (((long long) 1) << (shift[n_ones]))-1;
		shift[n_ones] -= 2;
		maskn[n_ones] = ((long long) 0);
		shiftn[n_ones] = 0;
		n_ones++;
		zero -= s;
	    } 
	    shift[n_ones] = 2*(one+zero);
	    mask[n_ones] = (((long long) 1) << (shift[n_ones]))-1;
	    shift[n_ones] -= 2;
	    maskn[n_ones] = (((long long)1) << (2*one))-1;
	}
	shiftn[n_ones] = 0;
	int j;
	for (j = 0; j < n_ones; j++) {
	    shiftn[j] += 2*one;
	}
	wd += zero+one;
	n_ones++;
    } while (1);
}

make_hash *HHH::setonepattern(const char *ppat, _hashtable* &htb, int &m_wsize, int &tailingzero)
{
    if (ppat == NULL) return NULL;
    char *pat = strsave(ppat);
    int z, one;
    z=one=0;
    char *p;
    tailingzero = 0;
    for (p = pat+strlen(pat)-1; p >=pat && *p == '0'; p--) {
	*p = 0;
	tailingzero++; // control the number of skips.
    }
    p = pat;
    while (*p) {
	if (*p == '0') z++;
	else if (*p == '1') one++;
	else fatal("pat has none 0/1 letters\n");
	p++;
    }
    m_wsize = one+z;
    int m_real_wsize = one;
    if (htb == NULL) htb = new _hashtable((1<<(2*m_real_wsize)));
    int z1, z2, o1, o2, z3, o3;
    z1 = z2=o1=o2 =z3=o3= 0;
    for (p = pat; *p; p++) {
	if (*p != '0') break;
	z1++;
    }
    for (; *p; p++) {
	if (*p != '1') break;
	o1++;
    }
    for (; *p; p++) {
	if (*p != '0') break;
	z2++;
    }
    for (; *p; p++) {
	if (*p != '1') break;
	o2++;
    }
    for (; *p; p++) {
	if (*p != '0') break;
	z3++;
    }
    for (; *p; p++) {
	if (*p != '1') break;
	o3++;
    }
    make_hash *hashitem;
    if (*p || (z3==0 && z2>30) || (z3+o3)>30 || (z3>0 && z2+o2>30)) {
	makehash_gen_discont *n = new makehash_gen_discont();
	n->set_param(pat);
	hashitem = n;
    } else if (z2 == 0) {
	makehash_cont *n = new makehash_cont();
	n->setwordsize(o1);
	hashitem = n;
    } else if (z3 ==0){
	if (z2+o2<=30) {
	    makehash_simple_discont *n = new makehash_simple_discont();
	    n->set_param(z1, o1, z2, o2);
	    hashitem = n;
	} else {
	    makehash_sdiscont *n = new makehash_sdiscont();
	    n->set_param(z1, o1, z2, o2);
	    hashitem = n;
	}
    } else {
	makehash_3one_discont *n = new makehash_3one_discont();
	n->set_param(z1, o1, z2,o2, z3,o3);
	hashitem = n;
    }
    hashitem->setpat(pat, m_wsize);
    return hashitem;
}


void HHHm::setpattern(char *pat, int &m_wsize, int &tzero) 
{
    N = 0;
    char *p, *q = pat;
    int m, t;
    while (q) {
	p = strchr(q, ',');
	if (p) {*p = 0; p++;}
	hashitem[N] = setonepattern(q, htb[N], m, t);
	if (N==0) {m_wsize = m; tzero = t; }
	else {
	    if (m_wsize != m || tzero != t) {
		fatal("multiple pattern not same length or tailing zeros\n");
	    }
	}
	N++;
	q = p;
    }
}

long long HHH1::findfirst()
{
    int h = hashitem->hash();
    return htb->findfirst(h);
}

long long HHH1::findnext()
{
    return htb->findnext();
}

long long  HHHm::findfirst()
{
    currentN = 0;
    long long x;
    do {
	int h =  hashitem[currentN]->hash();
	x = htb[currentN++]->findfirst(h);
    } while (x < 0 && currentN < N);
    return x;
}

long long  HHHm::findnext()
{
    long long  x = htb[currentN-1]->findnext();
    while (x < 0 && currentN < N) {
	int h =  hashitem[currentN]->hash();
	x = htb[currentN++]->findfirst(h);
    }
    return x;
}

#define mm 1

int probe_array_size = mm;

static int pm[mm*16] = {
    0, 10, 10, 10, 10, 0, 10, 10, 10, 10, 0, 10, 10, 10, 10, 0,
};

static int *probe_matrix = pm;
static int gap_def[mm]={20};
static int *gap_value = gap_def;

static void set(int *a, int i, int j, int x)
{
    int lowi = i+'a'-'A';
    int lowj = j+'a'-'A';
    if (i < 'A') lowi = i;
    if (j < 'A') lowj = j;
    a[128*i+j] = a[128*lowi+j] = a[128*i+lowj] = a[128*lowi+lowj] = x;
}

static int minimal(int amb1, int amb2, int *b)
{
    int i, j;
    int m = 0;
    for (i = 0; i < 4; i++) {
	if ((amb1 & (1 << i)) == 0) continue;
	for (j = 0; j < 4; j++) {
	    if ((amb2 & (1 << j)) == 0) continue;
	    if (m < b[i*4+j]) m = b[i*4+j];
	}
    }
    return m;
}

static int average(int amb1, int amb2, int *b)
{
    int i, j;
    int m = 0, c= 0;
    for (i = 0; i < 4; i++) {
	if ((amb1 & (1 << i)) == 0) continue;
	for (j = 0; j < 4; j++) {
	    if ((amb2 & (1 << j)) == 0) continue;
	    m += b[i*4+j];
	    c++;
	}
    }
    return m/c;
}

static void setarray(int *a, int *ar,  int *b, int aow)
{
    char tmp[] = "ACGTNVBHDKSWMYR";
    char tmp1[] ="TGCANBVDHMSWKRY";
    int ambcode[] = {1,2,4,8,15,7,14,11,13,12,6,9,3,10,5};
    int i;
    for (i = 0; i < 128*128; i++) a[i] = ar[i] = 1000; 
    for (i = 0; i < 15; i++) {
	int j;
	for (j = 0; j < 15; j++) {
	    int m;
	    if (aow) {
		m = minimal(ambcode[i], ambcode[j], b);
	    } else {
		m = average(ambcode[i], ambcode[j], b);
	    }
	    set(a, tmp[i], tmp[j], m);
	    set(ar, tmp1[i], tmp1[j], m);
	}
	int m = average(ambcode[i], 15, b);
	set(a, 'X', tmp[i],10000);
	set(a, tmp[i], 'X', 10000);
	set(ar, 'X', tmp1[i], 10000);
	set(ar, tmp1[i], 'X', 10000);
	set(a, '8', tmp[i], 0);
	set(a, tmp[i], '8', 0);
	set(ar, '8', tmp1[i], 0);
	set(ar, tmp1[i], '8', 0);
    }
}

void matrix::matrix_main(int wdsize)
{
    num = probe_array_size;
    int each = 128;
    int total = num*each*each;
    m = new int[total];
    mr = new int[total];
    gap_value = gap_def;
    int i;
    if (Scale != 1) {
	for (i = 0; i < num*16; i++) probe_matrix[i] *= Scale;
    } 
    for (i = 0; i < num; i++) {
	setarray(m+i*each*each, mr+i*each*each, probe_matrix+i*16, AOW);
    }   
}

matrix::~matrix()
{
    delete [] m;
    delete [] mr;
}

void matrix::matrix_read(FILE *fp, int wdsize)
{
    char line[1000];
    if (!find_string(fp, "Header", line)) {
	fatal("Does not find header in matrix file");
    }
    char order[128];
    int i;
    init_scode(order);
    char find = 0;
    num = 0;
    while (fgets(line, 1000, fp)) {
	if (strncmp(line, "End", 3) == 0) {
	    find = 1; break;
	} else if (strncmp(line, "NumMatrices", 11) == 0) {
	    char tmp[100];
	    sscanf(line, "%s %d", tmp, &probe_array_size);
	    num = probe_array_size;
	} else {
	    fprintf(stderr, "Matrix file ignore line \"%s\"", line);
	}
    }
    if (!find) {fatal("End of header not found");}
    if (num <= 0) {
	fatal("number of matrices wrong\n");
    }
    probe_matrix = new int[num*16];
    gap_value = new int[num];
    memset(gap_value, 0, sizeof(int)*num);
    find = 0;
    int cur = 0;
    while (fgets(line, 1000, fp)) {
	if (strncmp(line, "Matrix", 6) == 0) {
	    int n;
	    sscanf(line+6, "%d", &n);
	    n--;
	    if (n < cur) fatal("Matrix Order not increasing");
	    int *a = probe_matrix+16*n;
	    gap_value[n] = read_one_matrix(a, cur, cur, 'L', 0, n, 1, order, line, fp);
	    //a[0] = a[5] = a[10] = a[15] = 0;
	}
    }
    for (i = 1; i < num; i++) {
	if (gap_value[i] == 0) gap_value[i] = gap_value[i-1];	
    }
    matrix_main(wdsize);
}


void matrix::revising()
{
    zero_diagnal(m, mr, num);
}

int *matrix::ma(int a, int b, int len)
{
    if (a >= num) a = num-1;
    return m+(a*128+b)*128;
}

int *matrix::mar(int a, int b, int len)
{
    if (a >= num) a = num-1;
    return mr+(a*128+b)*128;
}

int matrix::gap(int a)
{
    if (a >=num) a = num-1;
    return gap_value[a];
}

void unique_hits_check::renew() {
    memset(num, 0, sizeof(int)*16);
}

void unique_hits_check::addhits(long long pos) 
{
    long long w = ( llabs(pos) & ((long long)  15)); 
    if (num[w] == MMLIST) fatal("outof memory for hit list\n"); 
    list[w][num[w]] = pos; num[w]++;
}

int unique_hits_check::is_new(long long pos, long long pos_e, int reversed, int dir, long long slen) 
{
    pos = position_convert(pos, pos_e, reversed, dir, slen, offset); 
    long long w = ( llabs(pos) & ((long long)  15));
    int i; 
    long long *l = list[w];
    for (i = 0; i < num[w]; i++) {
	if (l[i] == pos) return 0;
    }
    return 1;
}

int unique_hits_check::addline(char *line, int &mini_mis, int left)
{
    char *p = line;
    int c = 0;
    mini_mis = 1000;
    while (p = strchr(p, ',')) {
	p++;
	long long pos; int m; int st =0; 
	sscanf(p, "%lld.%d:(%*d.%*d.%d", &pos, &m, &st);
	if (m < mini_mis) mini_mis = m;
	addhits(pos-st+left);
	c++;
    }
    return c;
}

int localscore(int len, int mis, int mat, int numM)
{
    return len*mat+numM*(mis-mat);
}

