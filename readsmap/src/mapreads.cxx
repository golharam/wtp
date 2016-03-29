/* Notes
    5/21/09  Fix a bug in map that build hash.   
    5/19/09  Fix some bug in multi-run memory efficient mode . V2.3.1
    5/15/09  First with k-best schema setting V2.3
    5/7/09 Fix bigs on VA=1 in extension V2.1
    5/5/09   Add masking feature to automatically finding the right schema (only when using DBschema). V1.3.2
    4/20/09  Fix some bugs in multi-pass V1.3.1
    4/17/09  Add multiple pass when memory is limited. V1.3
    4/10/09  Add feature that allow local  to be anchored any place of the read  V 1.2.4
    2/16/09  Fflush(stdout) added after headerlines output
    2/13/09   Fix a bug in map multi-thread that turned off VA
    11/14/08 New m= option for local alignment 
    11/3/08 Add   r= option for both strands. 
    8/15/08 Fix setting V-option for all the runs.
	    Also fix return value for system call when multiple processes
		are run.
    8/8/08  Add a feature to allow a schema database
	    Also Fix a bug in dealing with dot.
    	    Fix a bug in dealing with VA at the start of the read	
    7/30/08 Add V-option V1.1.4
    1/28/09  Multi-thread version with t= option. V1.2.1 
*/

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <unistd.h>
#include <libgen.h>
#include <limits.h>
#include <string.h>
#include <ctype.h>
#include "zutil.h"
#include <sys/types.h>
#include <sys/time.h>
#include <time.h>
#include "fasta-io.h"
#include <unistd.h>

#include "zcompress.h"
#define MAX_SEQ 5000000
#define MAXSPLIT 100
char *path = NULL;
char tmpin[1000]; 
char tmpout[1000];
char tmpmatrix[1000];
char tmpseqsts[1000];
char tmpindex[1000];
char tmpreference[1000];
int doseq=0;
int compress_A = 0;  //compress is a global in zlib, so added _A
static int conv = 0, ref = 0;
static char *matrix = NULL;
static long long  firstr = 0, lastr = ((long long) -1);
static char *templatefile = NULL;
static FILE *tempfp = NULL;
static int adj = 0;
static char headerline[1000000];
static int pertaglimit = 1000;
static long long  starting = 0;
static int reportP = 1;
static char *codestring = NULL;
static int cleanflag = 1;
static int incremental = 0;
static int mul = 0;
static char *qfile = NULL;
static float snprate = 0.00000001;
static int timing = 0;
static int amballow = 0;
static int unique = 1;
static char tmpstsfile[1000];
static float AQV = 15.0; 
static int both = 1;
static char *extra_options = " ";
static int mismatch = 0;
static char *outputfile = NULL; 
static int Mlimit = 0; // in GBs 0 means no limit
static int num_of_run = 1;
static char split_genome[MAXSPLIT][1000];
static long long split_place[MAXSPLIT];
static char *masking = NULL;
static int numZero = 0;
static int kbest = 0;

int Usage()
{
	fprintf(stderr,"\nmapreads (2.3.2 5/21/2009) utility program\n");

	fprintf(stderr,"USAGE:  mapreads readsC seqfile [options]\n\n");
	fprintf(stderr,"OPTIONS:\n");
	fprintf(stderr, "\tM=##    Number of mismatches allowed (default 0)\n");
	fprintf(stderr, "\tP=str   Path whether to find program map (Def same as mapreads)\n");

	fprintf(stderr, "\tS=##    0 color space, 1 seq space, 2 both, 3 analyze seqs in color, 4 ana in seq space(0)\n");
	fprintf(stderr, "\tF=0/1   Wether to filter reads progressively when hits are found(0)\n");
	fprintf(stderr, "\tL=##    length of reads (Def 15)\n");
	fprintf(stderr, "\tX=str   matrix file (default empty)\n");
	fprintf(stderr, "\tC=0/1   clean up temp file or not (default 1)\n");
	fprintf(stderr, "\tT=str   template file (default empty)\n");
	fprintf(stderr, "\tA=0-2   2: count ajacent mismatches as 1, 1:only consist ones; 0: not counting them as 1 (Def 0)\n");
	fprintf(stderr, "\tO=##    offset of reference sequence(def 0)\n"); 
	fprintf(stderr, "\tZ=##    Maximum number of hits per tag(def 1000)\n");
	fprintf(stderr, "\tY=str   To change color code\n");
	fprintf(stderr, "\tI=0/1   Input reference sequence is single sequence (0) or multiple sequences (1)\n");
	fprintf(stderr, "\tb=##    first read to process (0-based)\n");
	fprintf(stderr, "\te=##    last read to process (0-based)\n");
        fprintf(stderr, "\tU=.##   SNP rate in reference\n");
        fprintf(stderr, "\tH=0/1   whether allow amb letter to be treated as known SNP(no penalty)\n");
        fprintf(stderr, "\tQ=qualfile   Optional quality values file for quality-weighted matching\n");
	fprintf(stderr, "\tB=##    The first base of reference sequence to be matched against(0-based)\n"); 
	fprintf(stderr, "\tE=##    The last starting base of reference sequence to be matched against\n");
	fprintf(stderr, "\tV=##.## Average QV value for the run\n");
	//fprintf(stderr, "\tm=#.##  mismatch penalty for local alignment, matching is 1. If m=0, no local alignment\n");
	fprintf(stderr, "\tr=0/1   Whether mapping to both strands(default 1)\n");
	fprintf(stderr, "\tout=file output file name, only effective when multi-thread\n"); 
	fprintf(stderr, "\tt=0/1   Turn on(1)/off time stamp\n");
	fprintf(stderr, "\tm=-#.## mismatch score for local alignment option\n");
	fprintf(stderr, "\tML=##   Total memory available in GBs\n"); 
	fprintf(stderr, "\tK=0/1   Whether to use k-best or not,if not use k-first\n"); 
	fprintf(stderr,"\n");
	return 1;
}

static void clean_up()
{
    char command[2000];
    if (cleanflag  == 1) {
        sprintf(command, "rm -f %s*", tmpin);
        system(command);
        sprintf(command, "rm -f %s", tmpseqsts);
        system(command);
        sprintf(command, "rm -f %s", tmpmatrix);
        system(command);
        sprintf(command, "rm -f %s", tmpout);
        system(command);
	if (mul) {
	    sprintf(command, "rm -f %s", tmpindex);
	    system(command);
	    sprintf(command, "rm -f %s*", split_genome[0]);
	    system(command);
	}
	if (firstr >0 || lastr != -1 || incremental) {
	    sprintf(command, "rm -f %s", tmpstsfile); 
	    system(command);
	}
	tmpin[1]='t';
        sprintf(command, "rm -f %s", tmpin);
        system(command);
    }
}
    
static int run_no = 0; 

static void runmap_core(char *command, char *pat, int compr)
{
    char com[20000], com1[1000];
    char m[1000];
    run_no++;
    fprintf(stderr, "map start run No. %d\n", run_no);
    if (timing) system("date");
    if (matrix) {
	sprintf(m, "X=%s B=0 m=%d %s", matrix, mismatch, extra_options);
    } else sprintf(m, "B=%d m=%d %s", both, mismatch, extra_options);
    if (codestring) {
	char n[1000];
	sprintf(n, " Y=\"%s\"", codestring);
	strcat(m, n);
    }
    if (qfile) {
        char n[1000];
        sprintf(n, " Q=\"%s\"", qfile);
	strcat(m, n);
    }
    char *qq;
    if (qq=strchr(pat, ',')) {
	char n[1000];
	*qq = 0;
	int x = atoi(qq+1);
	sprintf(n, " T=%d", x*10);
	strcat(m, n);
    }

    if (compr == 0 ) {
	if (timing) {
	    sprintf(com, "time %s u=%d r=%d n=%d Z=%d P=\"%s\" M=%d U=%f H=%d %s > %s.%d", command, unique, ref, reportP,  pertaglimit,pat, doseq, snprate, amballow, m,  tmpout, run_no);
	} else {
	    sprintf(com, "%s u=%d r=%d n=%d Z=%d P=\"%s\" M=%d U=%f H=%d %s > %s.%d", command, unique, ref, reportP,  pertaglimit,pat, doseq, snprate, amballow, m,  tmpout, run_no);
	}
    } else 
	sprintf(com, "%s u=%d r=%d n=%d Z=%d P=\"%s\" M=%d U=%f H=%d %s | gzip -3 -c -f >  %s.%d ; exit ${PIPESTATUS[%d]}", command, unique, ref, reportP,  pertaglimit,pat, doseq, snprate, amballow, m, tmpout, run_no, (run_no == 1)? 0:1);
    //reportP = 0;
    //fprintf(stderr, "%s\n", com);
    strcat(headerline, "#");
    strcat(headerline, com);
    strcat(headerline, "\n");
    int i = system(com);
    if (i) {
 	fprintf(stderr, "fail to execute command:\n\t%s\n", com); 
	clean_up();
	exit(1);
    }
    if (timing) system("date");
    system("sleep 1");
    sprintf(com1, "mv -f %s.%d %s", tmpout, run_no, tmpin);
    system(com1);
    if (!both || (!matrix && !codestring)) return; 
    if (compr == 0 ) {
        char *f = strchr(com, '>');
        sprintf(f, "R=1 > %s.%d", tmpout, run_no);
    } else {
        char *f = strstr(com, "| gzip ");
        sprintf(f, "R=1 E=stdin | gzip -3 -c -f > %s.%d ; exit ${PIPESTATUS[%d]}", tmpout, run_no, (run_no==1)? 0:1);
        if (run_no == 1 ) {
           char g[10000];
           sprintf(g, "gzip -dc %s | %s", tmpin, com);
	   /*
           char *h = strstr(g,"E=.Tmpfile");
           char *i = strchr(h, ' ');
           sprintf(h, "E=stdin%s", i);
	   */
           strcpy(com,g);
        } 
    }
    //fprintf(stderr, "%s\n", com);
    strcat(headerline, "#");
    strcat(headerline, com);
    strcat(headerline, "\n");
    i = system(com);
    if (i) {
        fprintf(stderr, "fail to execute command:\n\t%s\n", com);
        clean_up();
        exit(1);
    }
    system(com1);
}

static void runmap(char *command, char *pat, int compr)
{
    if (num_of_run == 1){runmap_core(command, pat, compr); return;}
    int i; 
    int first = 0;
    if (run_no == 0) first = 1;
    char comm[100000];
    char *atsts = strstr(command, "map"), *rest;
    //fprintf(stderr, "%s\n", command); //debug
    if (atsts == NULL) fatal("command wrong\n");
    //fprintf(stderr, "%s\n", atsts);// debug
    atsts = strchr(atsts, ' '); 
    if (atsts == NULL) fatal("command wrong\n");
    atsts = strchr(atsts+1, ' ');
    if (atsts == NULL) fatal("command wrong\n");
    *atsts = 0;
    if ((rest = strchr(atsts+1, ' ')) == NULL) fatal("command wrong\n");
    //fprintf(stderr, "%s---%s\n", command, rest); //debug
    for (i = 0; i < num_of_run; i++) {
	if (first && i > 0 && compr)  {
	    sprintf(comm, "gzip -dc %s | %s %s %s q=%lld E=stdin", tmpin, command, split_genome[i], rest+1, split_place[i]);
	} else { 
	    sprintf(comm, "%s %s %s q=%lld", command, split_genome[i], rest+1, split_place[i]);
	}
	// adjust start
	//fprintf(stderr, "runcore %s\n", comm); //debug
	runmap_core(comm, pat, compr);
    }
    *atsts = ' ';
}

static void runmap(char *command, char *pat)
{
    runmap(command, pat, 0);
} 

static void mapping_0mis(char *seq, char *stsfile, int cspace, int len, int dir, int numproc)
{
    char command[10000];
    char pat[100];
    int i;
    int n = 0;
    //if (qfile) n+=1;
    for (i = 0; i < len-13; i++) 
	pat[i] = '0';
    for (; i < len; i++) pat[i] = '1';
    pat[i] = 0;
    sprintf(command, "%s/map %s %s T=%d L=%d C=%d D=%d np=%d V=%f", path, stsfile, seq,n*10, len, cspace, dir, numproc, AQV);
    runmap(command, pat);
}

static void zero(char *p, int n)
{
    int i;
    for (i = 0; i < n; i++) p[i] = '0';
    p[i] = 0;
}


static void mapping_1mis(char *seq, char *stsfile, int cspace, int len, int filter, int dir, int numproc)
{
    char command[10000];
    int n = 1;
    //if (qfile) n+=1;
    sprintf(command, "%s/map %s %s T=%d L=%d C=%d E=%s F=%d D=%d np=%d V=%f", path, stsfile, seq, n*10, len, cspace, tmpin, filter, dir, numproc, AQV);
    if (len == 15) {
	runmap(command, "1111110111111");
	runmap(command, "110000111111111");
	runmap(command, "111111111000011");
    } else if (len == 14) {
	if (adj >= 1) {
	    runmap(command, "1111111111000");//skip 4, first and last
	    runmap(command, "11100001111111");
	    runmap(command, "11111100000111");
	} else {
            runmap(command, "1111111111000");//skip 4, first and last
            runmap(command, "11110001111111");
            runmap(command, "11111110001111");
	}
    } else if (len == 19) {
	runmap(command, "111111111111000000");
	runmap(command, "1111110000000111111");
    } else if (len == 20) {
	runmap(command, "0111101111111111");
	runmap(command, "11111111110000001111");
    } else if (len == 24) {
	runmap(command, "11111110001111110"); // skip every two bases, eq to 5 patterns; still cover all adj mismatches
    } else if (len >= 22) {
	char pat[400];
	zero(pat, len-22);
	strcat(pat, "111111101111111");
	runmap(command, pat);
    } else {
	fatal("No internal scheme for this read length and number of mismatches\n");
    }
}

static void mapping_2mis(char *seq, char *stsfile, int cspace, int len, int filter, int dir, int numproc)
{

    char command[10000];
    int n = 2;
    //if (qfile) n+=1;
    sprintf(command, "%s/map %s %s T=%d L=%d C=%d E=%s F=%d D=%d np=%d V=%f", path, stsfile, seq, n*10,len, cspace, tmpin, filter, dir, numproc, AQV);
    if (len == 24) {
	runmap(command, "11111111111100000");//same as 111111111111000000000000, 000000111111111111000000, and 000000000000111111111111 three patterns
	runmap(command, "11111100000011111100000"); // same as 111111000000111111000000, 000000111111000000111111 two patterns.
	runmap(command, "111111000000000000111111");
    } else if (len == 20) {
	runmap(command, "111111111111000"); //get 3,0^4_1^12; 0^8_1^12
	runmap(command, "1111000011111111000");// skip 4, so also 00001111000011111111
	runmap(command, "1111111100001111000"); //same, get 2;
	runmap(command, "11110000111100001111");
	runmap(command, "11110000000011111111");
	runmap(command, "11111111000000001111");
    } else if (len == 19) {
	runmap(command, "11111111111000"); //get 3,0^4_1^11_0^4; 0^8_1^11 as well
	runmap(command, "111100001111111000");// same as 2
	runmap(command, "111111100001111000"); // same as 2
	//runmap(command, "1111000011111110000,1111111100011110000");  
	runmap(command, "1111000011100001111");
	runmap(command, "1111000000011111111");
	runmap(command, "1111111100000001111");
	//runmap(command, "1111000011100001111,1111000000011111111");
	//runmap(command, "1111111100000001111");
    } else if (len == 14) {
	/*
	runmap(command, "00001111111111");
	runmap(command, "01110001111111");
	runmap(command, "01111110001111");
	runmap(command, "01111111110001");
	runmap(command, "10110110110111");
	runmap(command, "10111011011011");
	runmap(command, "10111101101101");
	runmap(command, "11010111011101");

	runmap(command, "11100111101011");
	runmap(command, "11101101010111");
	runmap(command, "11011011100111");
	runmap(command, "11011100111011");
	runmap(command, "11101010111101");
	*/	
	runmap(command, "00001111111111,01110001111111,01111110001111,01111111110001,10110110110111,10111011011011,10111101101101,11010111011101,11100111101011,11101101010111,11011011100111,11101010111101,11011100111011");

	/*
	runmap(command, "11011010111110");
	runmap(command, "00101111111110");
	runmap(command, "11110101011110");
	runmap(command, "11111111100010");
	runmap(command, "00111111111100");
	*/

	runmap(command, "11011010111110,00101111111110,11110101011110,11111111100010");
	runmap(command, "00111111111100");
    }  else if (len == 15) {
	runmap(command, "100001111111111");
	runmap(command, "111111110001110");
	runmap(command, "110111101110011");
	runmap(command, "010111011101111");
	runmap(command, "101111011110101");
	runmap(command, "110110111111100");
	runmap(command, "011110100111111");
	runmap(command, "111101001011111");
	runmap(command, "111011011111010");
	runmap(command, "011101111111001");
	runmap(command, "011011111010111");
	runmap(command, "101110111011011");
	runmap(command, "111100111100111");
	runmap(command, "111011101101101");
	runmap(command, "001111101111110");
	runmap(command, "111010010111111");
	runmap(command, "111101110110110");
	runmap(command, "101111110101011");
	runmap(command, "110111110011101");
    } else {
	fatal("No internal scheme for this read length and number of mismatches\n");
    }
}

static void mapping_3mis(char *seq, char *stsfile, int cspace, int len, int filter, int dir, int numproc)
{
    int n = 3;
    //if (qfile) n++;
    char command[10000];
    sprintf(command, "%s/map %s %s T=%d L=%d C=%d E=%s F=%d D=%d np=%d V=%f ", path, stsfile, seq, n*10, len, cspace, tmpin, filter, dir, numproc, AQV);
    if (len == 24) {
	runmap(command, "111111111111000"); // same as 4 patterns;
	runmap(command, "1111111100001111000"); // same as 3;
	runmap(command, "1111000011111111000"); // same as 3;
	runmap(command, "11111111000000001111000"); // same as 2;
	runmap(command, "11110000000011111111000") ; // same as 2;
	runmap(command, "11110000111100001111000"); // same as 2;
	runmap(command, "111100000000000011111111");
	runmap(command, "111100000000111100001111");
	runmap(command, "111100001111000000001111");
	runmap(command, "111111110000000000001111");
    } else if (len == 20) {
	//runmap(command, "00000000111111111111");
	runmap(command, "00001111000011111111");
	runmap(command, "00001111111100001111");
	runmap(command, "00001111111111110000"); //take first as well
	runmap(command, "00110011001111001111");
	runmap(command, "00110011110011111100");
	runmap(command, "00110011111100110011");
	runmap(command, "00111100001111110011");
	runmap(command, "00111100110000111111");
	runmap(command, "00111100111111001100");
	runmap(command, "00111111001100111100");
	runmap(command, "00111111110011000011");
	runmap(command, "11010101010101011011");
	runmap(command, "11011010101010100111");
	runmap(command, "11100110100110011101");
	runmap(command, "11101001011001101110");
	runmap(command, "11010101101010111010");
	runmap(command, "11011010010101110101");
	runmap(command, "11100110011011010110");
	runmap(command, "11101001100111101001");
	runmap(command, "11010101011110100101");
	runmap(command, "11011010101101011010");
	runmap(command, "11100110111001101001");
	runmap(command, "11101001110110010110");
	runmap(command, "11100111100101100110");
	runmap(command, "11101011011010011001");
	runmap(command, "11011101101001010101");
	runmap(command, "11011110010110101010");
    } else if (len == 19) {
	runmap(command, "0000111011111101110");
	runmap(command, "0011011101110111001");
	runmap(command, "0011101110001011111");
	runmap(command, "0011110111101110100");
	runmap(command, "0100001111011110111");
	runmap(command, "0100110110110111011");
	runmap(command, "0100111101101011101");
	runmap(command, "0111010011101011110");
	runmap(command, "0111011100010101111");
	runmap(command, "0111100101011101101");
	runmap(command, "0111101011010111010");
	runmap(command, "0111111000111110100");
	runmap(command, "0111111111100000011");
	runmap(command, "1001001011110111011");
	runmap(command, "1001010100101111111");
	runmap(command, "1001111111010010101");
	runmap(command, "1010011101011011110");
	runmap(command, "1010100111110111100");
	runmap(command, "1010111110101110010");
	runmap(command, "1011001001111100111");
	runmap(command, "1011110011011001011");
	runmap(command, "1011111010100101101");
	runmap(command, "1101011111111001000");
	runmap(command, "1101100010111011101");
	runmap(command, "1101101101100101110");
	runmap(command, "1101110110011100110");
	runmap(command, "1101111001001110011");
	runmap(command, "1110010111000101111");
	runmap(command, "1110011010011111001");
	runmap(command, "1110101011101100101");
	runmap(command, "1110101100111001011");
	runmap(command, "1110110001110010111");
	runmap(command, "1111000111101110001");
	runmap(command, "1111001110110010110");
	runmap(command, "1111110101001111000");
    } else if (len == 15) {
	runmap(command, "011011111101001");
	runmap(command, "111111100100101");
	runmap(command, "110101011101011");
	runmap(command, "111100111100011");
	runmap(command, "101111010001111");
	runmap(command, "101101101101110");
	runmap(command, "101101011010111");
	runmap(command, "111110110011001");
	runmap(command, "100011111011110");
	runmap(command, "110111011110100");
	runmap(command, "101110011101101");
	runmap(command, "001111111011100");
	runmap(command, "101111011011010");
	runmap(command, "110010011111011");
	runmap(command, "101011101011101");
	runmap(command, "001111101110101");
	runmap(command, "100100111011111");
	runmap(command, "000111111101011");
	runmap(command, "011011111010011");
	runmap(command, "110011110010111");
	runmap(command, "111001001011111");
	runmap(command, "010111011011101");
	runmap(command, "010110110111101");
	runmap(command, "010011100111111");
	runmap(command, "001101110111101");
	runmap(command, "100111000111111");
	runmap(command, "011111001101110");
	runmap(command, "110111101001011");
	runmap(command, "110110101111100");
	runmap(command, "111011010111100");
	runmap(command, "011101110001111");
	runmap(command, "110111110101010");
	runmap(command, "011100101111011");
	runmap(command, "101110100111011");
	runmap(command, "010110111100111");
	runmap(command, "110001111101101");
	runmap(command, "111011101110010");
	runmap(command, "111101001111001");
	runmap(command, "111111111000100");
	runmap(command, "010101011111110");
	runmap(command, "101110110110110");
	runmap(command, "011110111011010");
	runmap(command, "011001111110110");
	runmap(command, "111000111110101");
	runmap(command, "001010011111111");
	runmap(command, "011111010110011");
	runmap(command, "111100010111110");
	runmap(command, "111010110001111");
	runmap(command, "111010101101110");
	runmap(command, "101011011100111");
	runmap(command, "110101100110111");
	runmap(command, "011111100011110");
	runmap(command, "101001110111011");
	runmap(command, "111110001010111");
	runmap(command, "100111111110001");
	runmap(command, "111101111011000");
    } else if (len == 14) {
	runmap(command, "0011011111111");// 1
	runmap(command, "01011111110110");
	runmap(command, "01101011111110");
	runmap(command, "01111100111110");
	runmap(command, "01111111011100");
	runmap(command, "01111111101010");
	runmap(command, "1000111111111");//2
	runmap(command, "10111010111110");
	runmap(command, "1011110101111");//3
	runmap(command, "10111111101100");
	runmap(command, "10111111110010");
	runmap(command, "11010011111110");
	runmap(command, "11011101111010");
	runmap(command, "11011110111100");
	runmap(command, "11011111001110");
	runmap(command, "11100111011110");
	runmap(command, "11101101110110");
	runmap(command, "11101110101110");
	runmap(command, "11101111111000");
	runmap(command, "11110101111100");
	runmap(command, "11110110111010");
	runmap(command, "11110111100110");
	runmap(command, "11111001101110");
	runmap(command, "11111011011010");
	runmap(command, "11111011110100");
	runmap(command, "11111110010110");
	//runmap(command, "00011011111111"); 1
	runmap(command, "00101101111111");
	runmap(command, "00111110110111");
	runmap(command, "00111111001111");
	runmap(command, "00111111111001");
	// runmap(command, "01000111111111"); 2
	runmap(command, "01011101111101");
	//runmap(command, "01011110101111");
	runmap(command, "01011111011011");
	runmap(command, "01101110111011");
	runmap(command, "01101111010111");
	runmap(command, "01101111101101");
	runmap(command, "01110011101111");
	runmap(command, "01110101111011");
	runmap(command, "01110110011111");
	runmap(command, "01110111110101");
	runmap(command, "01111001011111");
	runmap(command, "01111010111101");
	runmap(command, "01111011110011");
	runmap(command, "01111101100111");
	runmap(command, "10010111111101");
	runmap(command, "10011101101111");
	runmap(command, "10011110111011");
	runmap(command, "10011111010111");
	runmap(command, "10100011111111");
	runmap(command, "10101110011111");
	runmap(command, "10101111101011");
	runmap(command, "10101111110101");
	runmap(command, "10110101110111");
	runmap(command, "10110110101111");
	runmap(command, "10110111011011");
	runmap(command, "10111001111011");
	runmap(command, "10111011011101");
	runmap(command, "10111011100111");
	runmap(command, "10111100111101");
	runmap(command, "11001011111011");
	runmap(command, "11001100111111");
	runmap(command, "11001111011101");
	runmap(command, "11001111100111");
	runmap(command, "11010101011111");
	runmap(command, "11010110110111");
	runmap(command, "11010111101011");
	runmap(command, "11011001110111");
	runmap(command, "11011010011111");
	runmap(command, "11011011101101");
	runmap(command, "11011111110001");
	runmap(command, "11100101101111");
	runmap(command, "11100110111101");
	runmap(command, "11100111110011");
	runmap(command, "11101001111101");
	runmap(command, "11101010110111");
	runmap(command, "11101011001111");
	runmap(command, "11101101011011");
	runmap(command, "11110000111111");
	runmap(command, "11110011010111");
	runmap(command, "11110011111001");
	runmap(command, "11110111001101");
	runmap(command, "11111010101011");
	runmap(command, "11111100001111");
	runmap(command, "11111100110011");
	runmap(command, "11111101010101");
	runmap(command, "11111101101001");
	runmap(command, "11111110011001");
	runmap(command, "11111110100101");
	runmap(command, "11111111000011");
    } else {
	fatal("No internal scheme for this read length and number of mismatches\n");
    }
}
#define  MAXSLINES 5000 
char *schema_line[MAXSLINES];
int num_of_lines;

static int findschema(FILE *tmpfile, int len, int mis)
{
    char line[10000];
    len -= numZero;
    int use_adj = adj;
    if (amballow) use_adj = 1;
    while (fgets(line, sizeof line, tmpfile)) {
	if (line[0] == '$') {
	    int a, b, c;
	    a = b = c = -1 ;
	    sscanf(line+1, "%d %d %d", &a, &b, &c);
	    if (c < 0) continue;
	    if (a == len && b == mis && (b<=1 || c == use_adj))
		return 1;
	}
    }
}

static void masking_schema_line(char *line)
{
    int len = strlen(masking);
    int i, j;
    char tmp[1000];
    strcpy(tmp, line);    
    for (i = 0, j = 0; i < len-1; i++) {
	if (masking[i+1] == '0') line[i] = '0'; 
	else {
	    line[i] = tmp[j];
	    j++;
	}
    }
    line[i] = 0;
}

static int findkbestschemas(FILE *tempfile, int len, int mis)
{
    int i, j;
    num_of_lines = 0;
    char line[10000];
    for (i = 0; i <= mis; i++) {
	if (findschema(tempfile, len,  i) != 1) fatal("do not find schema\n");
	while (fgets(line, sizeof(line), tempfile)) {
	    if (line[0] == '-'  || line[0] == '$') break;
	    if (line[0] == '#') continue; 
            char *p = strchr(line, '\n');
            *p = 0;
            if (num_of_lines >= MAXSLINES) fatal("too many lines of schema\n");
            if (masking) masking_schema_line(line);
	    char line1[10000];
	    sprintf(line1, "%s,%d", line, i);
            schema_line[num_of_lines++] = strsave(line1);
	}
	rewind(tempfile);
    } 
}

static void map(char *seq, char *stsfile, int n, int cspace, int len, int f, int dir, int numproc)
{
    char command[1000];
    if (adj && cspace == 1) cspace += adj ;
    int zero = 0;
    if (masking) {
	int i; 
	for (i = 0 ; i < len; i++) if (masking[i] == '0') zero++;  
    }
    if (templatefile /*&& !doseq*/) {
	//if (qfile) n += 1; 
        char line[10000];
	/*
	if (n > 0) {
	    sprintf(command, "%s/map %s %s T=%d L=%d C=%d E=%s F=%d D=%d", path, stsfile, seq, 0, len, cspace, tmpin, f, dir);
	    int i;
	    for (i = 0; i < 14; i++) line[i] = '1';
	    for (; i < len; i++) line[i] = '0'; 
	    line[i] = 0;
	    runmap(command, line, compress);
	}
	*/
	sprintf(command, "%s/map %s %s T=%d L=%d C=%d E=%s F=%d D=%d np=%d V=%f ", path, stsfile, seq, 10*n, len, cspace, tmpin, f, dir, numproc, AQV );
	if (tempfp == NULL) tempfp = ckopen(templatefile, "r");
	int count = 3;
	int run = 0; 
	num_of_lines = 0;
	while (fgets(line, sizeof line, tempfp)) {
	    if (line[0] == '#') {
		if (run == 0) {
		    if (strncmp(line, "#Schema data base", 17) == 0) {
			if (kbest == 1) {findkbestschemas(tempfp, len+1-zero, n);break;}
			int x;
			if ((x=findschema(tempfp, len+1-zero, n)) == 0) {
			    fatal("Does not find the right schema\n");
			} else if (x == 2) { // future, build internally success
			    break;
			}	
		    }
		    run = 1;
		} 
		continue;
	    }
	    if (line[0] == '-') break;
	    if (line[0] == '$') break;
	    run = 1;
	    char *p = strchr(line, '\n'); 
	    *p = 0;
	    if (num_of_lines >= MAXSLINES) fatal("too many lines of schema\n");
	    if (masking) masking_schema_line(line); 
	    schema_line[num_of_lines++] = strsave(line);
	}
	int i;
	for (i = 0; i < num_of_lines; i++) {
	    char *line = schema_line[i];
	    runmap(command, line, compress_A);
	    /*
	    if (count == 0 && unique == 0) {
		char  command1[1000];
		sprintf(command1, "%s/remduphits %lld %d %s %d %d %s > %s.out2", path, (long long) 0, pertaglimit,tmpin, 0, -1, "rd", tmpout);
		int r = system(command1);
		if (r != 0) {fatal("remduphits fail\n");}
		sprintf(command1, "mv -f %s.out2 %s", tmpout, tmpin);
		//fprintf(stderr, "%s\n", command1);
		system(command1);
		count = 3; 
	    } else count--;
	    */
	    if (compress_A) {
	      sprintf(command, "gzip -dc < %s | %s/map %s %s T=%d L=%d C=%d E=stdin F=%d D=%d np=%d V=%f ",tmpin, path, stsfile, seq, 10*n, len, cspace, f, dir, numproc, AQV);
            }
	}
	return;
    }
    if (n == 0) mapping_0mis(seq, stsfile, cspace, len, dir, numproc);
    else if (n == 1) mapping_1mis(seq, stsfile, cspace, len, f, dir, numproc);
    else if (n == 2) mapping_2mis(seq, stsfile, cspace, len, f, dir, numproc);
    else if (n == 3) mapping_3mis(seq, stsfile, cspace, len, f, dir, numproc);
    else {
	fatal("Internal schemes support only up to 3 mismatches\n");
    }
}

static void convert(char *file, char *name)
{
    if (conv == 0) {
	strcpy(name, file);
	return;
    }
    char color2seq[128][128];
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
   
    sprintf(tmpseqsts, "%s.tmpseqsts", tmpin); 
    sprintf(name, "%s", tmpseqsts);
    FILE *fw = ckopen(name, "w");
    FILE *fr = ckopen(file, "r");
    char line[100000];
    int outsize;
    while (fgets(line, sizeof line, fr)) {
	if (line[0] == '#' || line[0] != '>') continue;
	outsize = fprintf(fw, "%s", line);
	check_printf;
	if (fgets(line, sizeof line, fr) == NULL) break;
	char line2[1000];
	sscanf(line, "%s", line2);
	int l = strlen(line2);
	int i;
	char c = color2seq[line2[0]][line2[1]];
	outsize = fprintf(fw, "%c", c);
	check_printf;
	for (i = 2; i < l; i++) {
	    c = color2seq[c][line2[i]];
	    outsize = fprintf(fw, "%c", c);
	    check_printf;
	}
	outsize = fprintf(fw, "\n");
	check_printf;
    }
    fclose(fr);
    fclose(fw);
}

static int secondoftime()
{

    time_t clock, c1;
    struct tm *ttm;
    clock = time(&c1);
    ttm = localtime(&clock);
    char str[100];
    strftime(str, 100, "%s", ttm);
    return atoi(str);
}

static void matrixunit(FILE *f, int i)
{
    int outsize = fprintf(f, "Matrix %d\nDefault 0 10\nEnd\n", i);
    check_printf;
}

static void matrixzero(FILE *f, int i)
{
    int outsize = fprintf(f, "Matrix %d\nDefault 0 0\nEnd\n", i);
    check_printf;
}

void process_matrix(char * &a, int sec)
{
    if (a == NULL) return;
    char *p = a;
    int zero = 0;
    while (*p) {
	if (*p != '0' && *p != '1') return;
	if (*p == '0') zero++; 
	p++;
    }
    masking = a;
    numZero = zero;
    sprintf(tmpmatrix, ".tmpfile.matrix.%dXXXXXX", sec);
    int fd  = mkstemp(tmpmatrix);
    close(fd);
    FILE *fp = ckopen(tmpmatrix, "w");
    int len = strlen(a);
    int outsize = fprintf(fp,  "Header\nNumMatrices %d\nEnd\n", len);
    check_printf;
    int i;
    for (i = 0; i < len; i++) {
        if (a[i] == '0') {
            matrixzero(fp, i+1);
        } else {
            matrixunit(fp, i+1);
        }
    }
    fclose(fp);
    a = tmpmatrix;
}

long long make_ref_bydot(char * &seqfile, int sec, long long be, long long end, int readlen, long long each )
{
    sprintf(tmpindex, ".tmpfile.index%dXXXXXX", sec);
    int fd = mkstemp(tmpindex);
    long long en = end;
    long long beg = be;
    if (en > be + each) en = be+each;
    close(fd);
    FILE *fp = ckopen(tmpindex, "w");
    char line[1000];
    long long start = -1;
    FastaFile fafile(SEQTYPE_NT);
    en += readlen;

    if (!fafile.Open(seqfile,"r"))
        fatal("Fail to open reference file\n");
    FastaSeq faseq;
    seqfile = split_genome[0];
    sprintf(seqfile, ".tmpfile.reference%dXXXXXX", sec);
    fd = mkstemp(seqfile);
    close(fd);
    FILE *fpw = ckopen(seqfile, "w");
    long long  total = 0;
    int num_seq = 0;
    long long all_total = 0;
    int outsize = fprintf(fpw, ">dot file\n");
    check_printf;
    int npass = 0;
    split_place[0] = 0;
    while (fafile.Read(faseq)) {
	long long start1 = be-all_total; 
	all_total += faseq.Length();
	if (all_total < be) continue;
        total += faseq.Length()+1;
	if (num_seq > MAX_SEQ) fatal("too many sequences in reference file\n"); 
        outsize = fprintf(fp, "%lld\n", total);
        num_seq++;
        check_printf;
	while (all_total >= be) {
	    if (start < 0) start = start1;
	    if (all_total >= en) {
                long long len = faseq.Length()-start1+en-all_total;
                char *temp = new char[len+1];
                strncpy(temp, faseq.Sequence()+start1,len);
                temp[len] = 0;
                outsize = fprintf(fpw, "%s\n", temp);
		delete [] temp;
		start1 += en-readlen-be;
		be = en-readlen; 
		en = en+each;
		if (en > end) en = end; 
		npass++;
		split_place[npass] = be-beg;
		fclose(fpw);
		sprintf(split_genome[npass], "%s.%d", seqfile, npass);
		fpw = ckopen(split_genome[npass], "w"); 
		fprintf(fpw, ">seq%d\n", npass);
	    } else {
		outsize = fprintf(fpw, "%s.\n", faseq.Sequence()+start1);
		be = all_total;
		check_printf;
		break;
	    } 
	    check_printf;
	    if (all_total >= end) goto exitpoint;
	}
    }
/*
    while (fafile.Read(faseq)) {
	all_total += faseq.Length();
	if (all_total >= en) {
	    int len = faseq.Length()+en-all_total;
	    char *temp = new char[len+1]; 
	    strncpy(temp, faseq.Sequence(),len);  
	    temp[len] = 0;
            outsize = fprintf(fpw, "%s.\n", temp);
            check_printf;
	    delete [] temp;
            total += len+1;
            outsize = fprintf(fp, "%lld\n", total);
            check_printf;
	    num_seq++;
	    if (num_seq > MAX_SEQ) fatal("Too many sequences in reference file for this run of\n");
	    break;
	}
	outsize = fprintf(fpw, "%s.\n", faseq.Sequence());
	check_printf;
	total += faseq.Length()+1;
	outsize = fprintf(fp, "%lld\n", total);
	check_printf;
    }
*/
  exitpoint:
    fclose(fpw);
    fclose(fp);
    fafile.Close();
    return start;
}
static long long  make_refs(char * &seqfile, int sec, long long be, long long en, int readlen)
{
    //fprintf(stderr, "Mlimit =%d\n", Mlimit);
    if (Mlimit > 25) Mlimit = 25; // enough space to process 4G genome, this will limit each run to 4G genome 
    if (Mlimit <= 0) return make_ref_bydot(seqfile, sec, be, en, readlen, en-be);
    long long s = ((en-be) >> 20)+1;
    s = s*6;
    //fprintf(stderr, "s=%lld\n", s);
    if (s < (Mlimit-1)*1024) return make_ref_bydot(seqfile, sec, be, en, readlen, en-be);
    int fs = MBofFile(seqfile)+1;
    //fprintf(stderr, "Input read file is %dMB\n", fs);
    fs = fs-(be >> 20)+1;
    int ffs = fs*6;
    if (ffs < (Mlimit-1)*1024) return make_ref_bydot(seqfile, sec, be, en, readlen, en-be);
    if (ffs > s) ffs = s;
    int n = ffs/(Mlimit-1)/1024+1;
    if (n > MAXSPLIT) fatal("the genome File too big to process\n"); 
    num_of_run = n;
    long long each = (((long long ) fs) << 20) / n+1;
    //fprintf(stderr, "size of each %lld\n", each);
    return make_ref_bydot(seqfile, sec, be, en, readlen, each); 
}

static int count(char *line)
{
    int i = 0;
    char *s = line;
    while (s = strchr(s, ',')) {
	i++;
	s++;
    }
    return i;
}

static void process_line(char *line, char *line1)
{
    int maxs = 0;
    char *s = line;
    s = strchr(line, '\n'); *s = 0;
    printf("%s", line);
    s = line;
    while (s = strchr(s, '.')) {
	s++;
	int x = atoi(s);
	if (x > maxs) maxs = x;
    }
    char *t = strchr(line1, ',');
    while (t) {
	s = strchr(t, '.');
	if (!s) break;
	int x = atoi(s+1);
	if (x > maxs) {
	    char *tt = strchr(t, ',');
	    int j = strspn(s, ",\n");
	    s[j] = 0;
	    printf("%s", t);
	    t = tt;
	} else {
	    t = strchr(t, ',');
	}
    }
}

int split(char *res, char *outsts, int Zlimit, char *insts)
{
    FILE *fres = ckopen(res, "r");
    FILE *fout = ckopen(outsts, "w");
    FILE *fp = NULL;
    int rv = 0;
    if (insts) fp = ckopen(insts, "r");
    char line[1000000], line1[1000000];
    while (fgets(line, sizeof line, fres)) {
	if (line[0] == '#') {
	    if (!fp) printf("%s", line);
	    continue;
	}
	int i = 0;
	if (fp) {
	    while (fgets(line1, sizeof line1, fp)) {
		if (line1[0] != '#') {i = 1; break;}
	    }
	}
	if (count(line) >= Zlimit) {
	    fprintf(fout, "%s", line);
	    if (fgets(line, sizeof line, fres)==NULL) break;
            fprintf(fout, "%s", line);
	    rv = 1;
	} else {
	    if (i) process_line(line, line1);
	    else  printf( "%s", line);
            if (fgets(line, sizeof line, fres)==NULL) break;
            printf("%s", line);
	}
	if (fp) fgets(line1, sizeof line1, fp);
    }
    return rv;
}

main(int argc, char *argv[])
{

    char *stsfileIn = NULL;
    char *stsfile = new char[10000];  //fn without gz
    char *seqfile = NULL;
    int numMis = 0;
    int readlen = 15;
    int filter = 0;
    int mode = 0;
    int chrom_start = 0;
    int i;
    long long seqbeg=0, seqend=(((long long) 1) << 60);
    int compact = 0;
    int numproc = 1;
    for (i=1; i<argc; ++i)
    {
	if (argv[i][1] == '=')         // X=value
	{
	    if (argv[i][2] == 0)
		fprintf(stderr,"Missing value for %s\n",argv[i]);
	    else if (argv[i][0] == 'M')
		numMis = atoi(argv[i]+2);
	    else if (argv[i][0] == 'P') 
		path = argv[i]+2;
	    else if (argv[i][0] == 'S') 
		mode = atoi(argv[i]+2);
	    else if (argv[i][0] == 'F') 
		filter = atoi(argv[i]+2);
	    else if (argv[i][0] == 'L') 
		readlen = atoi(argv[i]+2);
	    else if (argv[i][0] == 'X') 
		matrix = argv[i]+2;
	    else if (argv[i][0] == 'T') 
		templatefile = argv[i]+2;
	    else if (argv[i][0] == 'C') 
		cleanflag = atoi(argv[i]+2);
	    else if (argv[i][0] == 'A') 
		adj = atoi(argv[i]+2);
	    else if (argv[i][0] == 'Z') 
		pertaglimit = atoi(argv[i]+2);
	    else if (argv[i][0] == 'O') {
		starting = atoll(argv[i]+2);
		/*
		int x = 1, y=-7;
		printf("%lld %lld\n", x+starting, y-starting);
		*/
	    }
	    else if (argv[i][0] == 'Y')
                codestring = argv[i]+2;
	    else if (argv[i][0] == 'c') 
		conv = atoi(argv[i]+2); 
	    else if (argv[i][0] == 'R') 
		ref = atoi(argv[i]+2);
	    else if (argv[i][0] == 'I') 
		mul = atoi(argv[i]+2);
	    else if (argv[i][0] == 'b')
                firstr = atoll(argv[i]+2);
            else if (argv[i][0] == 'e')
                lastr = atoll(argv[i]+2);
	    else if (argv[i][0] == 'Q') 
		qfile = argv[i]+2;
	    else if (argv[i][0] == 'U') 
		snprate = atof(argv[i]+2);
	    else if (argv[i][0] == 'H') 
		amballow = atoi(argv[i]+2);
	    else if (argv[i][0] == 'i')
		incremental = atoi(argv[i]+2);
	    else if (argv[i][0] == 'B') 
		seqbeg = atoll(argv[i]+2);
	    else if (argv[i][0] == 'E')
		seqend = atoll(argv[i]+2); 
	    else if (argv[i][0] == 'q') 
		compact = atoi(argv[i]+2);
	    else if (argv[i][0] == 'u') 
		unique = atoi(argv[i]+2);
	    else if (argv[i][0] == 'V') 
		AQV =  atof(argv[i]+2);
	    else if (argv[i][0] == 'K')
		kbest = atoi(argv[i]+2);
	    else if (argv[i][0] == 'm') 
		mismatch = (int)(atof(argv[i]+2)*10.0-0.5);
	    else if (argv[i][0] == 'x') 
		extra_options = argv[i]+2;
	    else if (argv[i][0] == 'r') 
		both = atoi(argv[i]+2); 
	    else if (argv[i][0] == 't')
		timing = atoi(argv[i]+2);
	} else if (argv[i][0] == 'n' && argv[i][1] == 'p' && argv[i][2] == '=') {
	  numproc = atoi(argv[i]+3);
	} else if (strncmp(argv[i], "out=", 4)== 0) {
	   	outputfile = argv[i]+4;
	} else if (strncmp(argv[i], "ML=",3) == 0) {
		Mlimit = atoi(argv[i]+3);  
	} else   // filename
	{
	    if (stsfileIn == NULL) {
		stsfileIn = argv[i];
		//remove_gz_ext(stsfileIn,stsfile);
		strcpy(stsfile,stsfileIn);  //It works now when file has .gz in it
	    }
	    else if (seqfile ==NULL)
		seqfile = argv[i];
	    else 
		fprintf(stderr,"Argument \"%s\" ignored\n",argv[i]);
	}
    }
    if (unique >= 2)  { unique = 1; compress_A =1;} 
    else unique = 1;
    if (numproc > 1) {compress_A = 0; unique = 1;} // not using compress in multi-thread
    if (seqfile == NULL) {
	return Usage();
    }
    if (path == NULL) {
      char tmp[64];
      /* Get the executable name */
      if(sprintf(tmp, "/proc/%d/exe", (int) getpid()) < 0)
	abort();
      path = (char *) malloc(PATH_MAX + 1);
      i = readlink(tmp, path, PATH_MAX);
      if(i < 0)
	{
	  fprintf(stderr,"Failed to readlink(%s): %s",tmp,strerror(errno));
	  exit(1);
	}
      path[i] = 0; // readlink does not null terminate returned string
      path = dirname(path);
    }
    
    //system("rm -f out1");
    sprintf(tmpindex, "rd");
    int sec = secondoftime();
    sprintf(tmpin, ".tmpfile%dXXXXXX", sec);
    int fd =mkstemp(tmpin);
    close(fd);
    tmpin[1]='T';
    sprintf(tmpout, "%s.out", tmpin);
    headerline[0] = 0;
    if (firstr >0 || lastr > 0 || incremental) {
	sprintf(tmpstsfile, "%s.sts", tmpin);
	char endp[1000];
	if (lastr > 0) {
	    sprintf(endp, "No. %d", lastr);
	} else {
	    sprintf(endp, "end");
	}
	if (incremental) {
	    if (firstr > 0 || lastr > 0)
	    	sprintf(headerline, "#tempfile %s contains reads %d to %s from file %s that are not mapped previously\n", tmpstsfile, firstr, endp, stsfile);
	    else sprintf(headerline, "#tempfile %s contains all reads from file %s that are not mapped previously\n", tmpstsfile, stsfile);
	} else { 
	    sprintf(headerline, "#tempfile %s contains reads No. %d to %s of read file %s\n",  tmpstsfile, firstr, endp, stsfile);
	}
	FILE *fp = ckopen(tmpstsfile, "w");
	//FILE *fpr = ckopen(stsfile, "r");
	genFile fpr;
	fpr.setUsePipe();
	fpr.open(stsfile, "r");
	char line[100000];
	long long c = 0;
	if (firstr > 0) {
	    while (fpr.gets(line, sizeof line)) {
		if (line[0] != '>') continue;
		c++;
		if (c == firstr) {
		    fpr.gets(line, sizeof line);
		    break;
		}
	    }
	}
	if (incremental) {
	    while (fpr.gets(line, sizeof line)) {
		if (line[0] == '#') {
		    fprintf(fp, "%s", line);
		    continue;
		}
		if (strchr(line, ',')) {
		    if (!fpr.gets(line, sizeof line)) break;
		} else {
		    fprintf(fp, "%s", line);
            	    if (fpr.gets(line, sizeof line)) {
                	fprintf(fp, "%s", line);
            	    } else break;
		}
                if (c == lastr) break;
		c++;
	    }
	} else {
	while (fpr.gets(line, sizeof line)) {
  	    fprintf(fp, "%s", line);
	    if (line[0] == '#') continue;
	    if (fpr.gets(line, sizeof line)) {
            	fprintf(fp, "%s", line);
	    } else break;
	    if (c == lastr) break;
	    c++;
	}
	}
	fclose(fp);
	fpr.close();
	stsfile = tmpstsfile;
    }
    char command[10000];
    sprintf(command, "rm -f %s", tmpin);
    system(command);
    process_matrix(matrix, sec); 
    if (mul) {
	chrom_start = starting;
	starting =  make_refs(seqfile, sec, seqbeg, seqend, readlen) ;//change the file and replace thename 
    }
    if (mode > 2) {
	doseq = 1;
	if (mode == 3) 
	    map(seqfile, stsfile, numMis, 1, readlen-1, filter, 1, numproc);
	else 
	    map(seqfile, stsfile, numMis, 0, readlen, filter, 1, numproc);
	if (compress_A == 0) 
	    sprintf(command, "%s/remduphits %lld %d %s %d %d %s", path, starting, pertaglimit, tmpin, mul, chrom_start, tmpindex);
	else 
	    sprintf(command, "gzip -dc < %s | %s/remduphits %lld %d stdin %d %d %s", tmpin, path, starting, pertaglimit, mul, chrom_start, tmpindex);
    } else {
	if (/*numMis > 0 && mode ==*/ 0) {
            char name[10000], tmpsts[10000], *old_sts = NULL, *new_sts = NULL;
            sprintf(name, "%s.sts", tmpout);
	    sprintf(tmpsts, "%s.sts1", tmpout);
	    new_sts = stsfile;
	    do {
	        map(seqfile, new_sts, numMis, 1, readlen-1, filter, 1, numproc);
     	    	sprintf(command, "%s/remduphits %lld %d %s %d %d %s %s > %s", path, starting, pertaglimit,tmpin, mul, chrom_start, tmpindex, new_sts, tmpout);
		if (old_sts == NULL) {
	    	    if (split(tmpout, name, pertaglimit, old_sts)== 0) break; 
		    new_sts = name;
		    old_sts = tmpsts;
		    if (numMis > 1) numMis = 1; else numMis--;
		} else {
		    if (split(tmpout, old_sts, pertaglimit, new_sts) == 0) break;
		    char *t = new_sts;
		    new_sts = old_sts; 
		    old_sts = t;
		    numMis--;
		}
	    } while( numMis >= 0);
	} else {
	  if (numMis == 0) {
	    map(seqfile, stsfile, 0, 1, readlen-1, filter, 1, numproc);
	  } else {
	    if (mode == 0 || mode == 2) {
	        map(seqfile, stsfile, numMis, 1, readlen-1, filter, 1, numproc);
	    }
	    if (mode == 1 || mode == 2) {
		char name[1000];
		convert(stsfile, name);
		map(seqfile, name, numMis, 0, readlen, filter, 1, numproc);
	    }
	  }
 	  if (compact) pertaglimit *= -1;

	  if (compress_A == 0) {
	    if (numproc > 1) {
		int i;
    		int outsize = printf("%s", headerline);
    		check_printf;
		fflush(stdout);
		char *p = stsfile, *q, *out, *qq; 
		out = outputfile;
		for (i = 0; i < numproc; i++) {
		    if (q = strchr(p, ':')) {
			*q = 0;
		    }
	   	    if (outputfile && (qq = strchr(out, ':'))) {
			*qq = 0;
		    } 
		    if (outputfile) {
			FILE *fp = ckopen(out, "w");
			fprintf(fp, "%s", headerline);
			fclose(fp);
			sprintf(command, "%s/remduphits %lld %d %s_%d %d %d %s %s >> %s", path, starting, pertaglimit,tmpin, i, mul, chrom_start, tmpindex, p, out); 
		    } else 
			sprintf(command, "%s/remduphits %lld %d %s_%d %d %d %s %s", path, starting, pertaglimit,tmpin, i, mul, chrom_start, tmpindex, p); 
		    int sts = system(command);
		    if (sts) fatal("remove redundant hits fails\n");
		    if (q) {
			p = q+1;
			if (outputfile && !qq) fatal("Number of output files fewer than the number of input files\n");
			out = qq+1; 
		    } else break;
		}
		clean_up();
		exit(0);
	    } else {
		sprintf(command, "%s/remduphits %lld %d %s %d %d %s %s", path, starting, pertaglimit,tmpin, mul, chrom_start, tmpindex, stsfile);
	    }
	  } else { 
		sprintf(command, "gzip -dc < %s | %s/remduphits %lld %d stdin %d %d %s %s", tmpin, path, starting, pertaglimit, mul, chrom_start, tmpindex, stsfile);
	  }
	}
    }
    int outsize = printf("%s", headerline);
    check_printf;
    fflush(stdout);
    int sts = system(command);
    if (sts) fatal("remove redundant hits fails\n");
    clean_up();
}
