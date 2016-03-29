#include "map.h"
#include "fasta-io.h"
// V1.5 Fix a bug in hash building
// V1.4     VA=1 in extension in local
// V 1.3    Add a offset parameter option
// Version 1.2 First version of seeded local 
// Version 1.1 fix some bugs in IUB mapping. Leading base is IUB code, also 
//       some cases of adjacent IUBs
// Version 1.0 mostly multi-threads workin

int Usage()
{
	fprintf(stderr,"\nmap (1.5 05/21/09) utility program\n");

	fprintf(stderr,"USAGE:  map readsfile seqfile [options]\n\n");
	fprintf(stderr,"OPTIONS:\n");
	fprintf(stderr,"\tT=##     Threshold Score\n");
	fprintf(stderr,"\tW=##     Word size \n");
	fprintf(stderr, "\tA=0/1   whether output alignment (default 0)\n");
	fprintf(stderr, "\tX=file  this gives the input matrix file\n");
	fprintf(stderr, "\tG=#     this gives the maximum number of indel bases allowed (default 0)\n");
	fprintf(stderr, "\tP=str   discontig word pattern(default11111111)\n");

	fprintf(stderr, "\tL=##    lenght of reads(dafault 15)\n");
	fprintf(stderr, "\tE=file  Intermediate outfile, add the hits to it\n");
	fprintf(stderr, "\tF=0/1   Wether to filter reads that already have hits in E-file\n");
	fprintf(stderr, "\tC=0/1   Whether it is in color space\n");
	//fprintf(stderr, "\tD=+/-1  Direction of the search for color space\n");
	fprintf(stderr, "\tB=0/1   Whether to search for both direction in seqeunce space\n");
	fprintf(stderr, "\tM=0/1   Default 0=mapping reads, 1=compare two sequences (maybe 1) at length L\n");
	fprintf(stderr, "\tZ=##    maximum number of hits reported per tag (default 1000) \n"); 
	fprintf(stderr, "\tn=0/1   whether report perfect match(1)\n");
	fprintf(stderr, "\tR=0/1   whether reversed complement ref seq(0)\n");
	fprintf(stderr, "\tQ=file  quality value file (NULL)\n");
	fprintf(stderr, "\tY=str   Color code string\n");
	fprintf(stderr, "\tb=##    first read to process\n");
	fprintf(stderr, "\te=##    last read to process\n"); 
	fprintf(stderr, "\tU=.##   SNP rate in reference\n");
	fprintf(stderr, "\tH=0/1   whether allow amb letter to be treated as known SNP(no penalty)\n");
	fprintf(stderr, "\tx=##    not report hits with x+m mismatches where m is minimum mismatch in any hit reported so far\n");
	fprintf(stderr, "\tm=-##   mismatch penalty for local alignment m=0, default no local alignment\n"); 
	fprintf(stderr, "\ts=##    start of the seeded area for local\n"); 

	fprintf(stderr,"\n");
	return 1;
}

long long  buildlookup(mapping_machine *e_PCR, const char *seqfile, int mode)
{
 FastaFile fafile(SEQTYPE_NT);
long long len= (long long) 1;
    if (!fafile.Open(seqfile,"r"))
        return 0;
    FastaSeq faseq;
        if (fafile.Read(faseq))
        {
            //len = e_PCR->init_hash(faseq.Label(),faseq.Sequence());
	    e_PCR->ProcessSeq(faseq.Label(),faseq.Sequence());
        }
        fafile.Close();
	return len;
}

int run_pcr(mapping_machine *e_PCR, const char *stsfile, const char *seqfile, int mode)
{
    FastaFile fafile(SEQTYPE_NT);

    if (!fafile.Open(seqfile,"r"))
	return 0;
    FastaSeq faseq;
    if (mode == 0) {
	if (fafile.Read(faseq))
	{
	    e_PCR->ProcessSeq(faseq.Label(),faseq.Sequence());
	}
	fafile.Close();

	return e_PCR->ReadprobeFile(stsfile);
    } else {
	if (fafile.Read(faseq))
	{
	    e_PCR->ProcessSeq(faseq.Label(),faseq.Sequence());
	} else return 0;
	if (strcmp(stsfile, seqfile) != 0) {
	    fafile.Close();
	    if (!fafile.Open(stsfile,"r"))
		return 0;
	    if (!fafile.Read(faseq)) return 0;
	}
	return e_PCR->buildprobe(faseq.Sequence());
	fafile.Close();
    }
}

mapping_machine **listmachine;
char **filename;
const char *refseq;
long long *begp; 

void *process_run_buildhash(void *targ)
{
    long long process_num = (long long) targ;
    listmachine[process_num]->build_hashing(begp[process_num], begp[process_num+1]-1);
}

void *process_run(void *targ)
{
    long long process_num = (long long) targ;
    if (listmachine[process_num]->ReadprobeFile(filename[process_num])) {
	listmachine[process_num]->moveOuttoTemp();
    } else exit(1);
}
#define num_mux 1024 
int multi_thread(const char *stsfile, const char *seqfile, int mode, int numproc, long long len)
{
    long long i;
/*  // multi-threading the hash building. Turn off for now, may be worth it when have 16 core 
    pthread_t th_id[numproc];
    pthread_mutex_t mutex[num_mux];
    for (i = 0; i < num_mux; i++) pthread_mutex_init(&mutex[i], NULL); 
    begp = new long long [numproc+1];
    begp[0] = 0; begp[1] = len/numproc;
    for (i =2; i < numproc; i++) {
	begp[i] = begp[i-1]+len/numproc;
    }
    begp[i] = len;
    for (i = 0; i < numproc; i++) {
	listmachine[i]->getHash()->setMutex(mutex, num_mux);
        pthread_create(&th_id[i], NULL, &process_run_buildhash, (void *) i);
    } 

    for (i = 0; i < numproc; i++) {
        pthread_join(th_id[i], NULL);
    }
*/
    filename = new char*[numproc];
    refseq = seqfile;
    char *s = strsave(stsfile), *t; 
    filename[0] = s;
    i = 1;
    while (s = strchr(s, ':')) {
	*s = 0; 
	s++;
	if (i < numproc)
	    filename[i++] = s;
	else break;
    }
    if (i < numproc) numproc = i;
    // launch multi-thread
    fprintf(stderr, "start %d threads\n", numproc);
    pthread_t thread_id[numproc];
    for (i = 0; i < numproc; i++) {
	pthread_create( &thread_id[i], NULL, &process_run, (void *)i );
    }
    int return_state = 1;
    for (i = 0; i < numproc; i++) {
	int *x;
	pthread_join( thread_id[i], NULL);
    }
    //if (return_state != 1) return 0; 
    //fprintf(stderr, "all threads done\n");
    return 1;
}

        const char *matrixfile = NULL;
        int threshold = 0, wdsize= 0;
        int program = 0;
        int gap = 0;
        char *mfile = NULL;
        char islocal = 0, dis=0;
        int dw = 0;
        char *pat = NULL;
        int len = 15;
        char *tempout = NULL;
        int filter = 0;
        int colorspace = 0;
        int dir = 1;
        int both = 1;
        int shortout = 1;
        int mode = 0;
        int hitlimit = 1000;
        int reportP = 1;
        char *qfile = NULL;
        char *cc = NULL;
        double AQV = 15.0;
        int Re = 0;
        int ref = 0;
        float SNPrate = 0;
        int amballow = 0;
        int uniqueH = 0;
        int mismatch = 0;
        char *outfile = NULL;
        int pack = 0;
        long long  first_read = 0;
        long long  last_read = ((long long) 1) << 62;
        int numproc = 1;
        int xdrop = 1000;
	int prefix_len = 0;
	long long offset = 0;
	int num_ext = 0;

static void buildmachine(mapping_machine * &e_PCR, mapping_machine *copy)
{
        if (colorspace >= 1) {
            delete e_PCR;
            if (colorspace == 1 && SNPrate <= 0.00001 && amballow == 0) {
                mapping_machine_color *mp = new mapping_machine_color();
                e_PCR = mp;
            } else {
                if (colorspace<=2) {
                    mapping_machine_color_adj *mp = new mapping_machine_color_adj();

                    mp->setAverageQV(AQV);
                    if (colorspace == 2) {
                        mp->setSNPrate(2.0, amballow);
                    } else {
                        mp->setSNPrate(SNPrate, amballow);
                    }
                    e_PCR = mp;
                } else {
                    mapping_machine_color_adja *mp = new mapping_machine_color_adja();
                    e_PCR = mp;

                }
            }
        }
        if (Re == 1) both = 0;
        if (cc) e_PCR->setcolorcode(cc);
        if (qfile) e_PCR->setQfile(qfile);
        e_PCR->setgapsize(gap);
        e_PCR->probe_init(len);

	if (copy) {
            e_PCR->share(copy);
        }
        e_PCR->SetWordSize(pat);
        e_PCR->setMatrix(mfile);
	e_PCR->setdisplay(dis);
        e_PCR->setNoperfect(!reportP);
        // fprintf(stderr, "matrix built\n");
        e_PCR->SetThreshold(threshold);
        e_PCR->setTempout(tempout);
        e_PCR->setFilter(filter);
        e_PCR->setColorspace(colorspace);
        e_PCR->setDir(dir);
        e_PCR->setBoth(both);
	e_PCR->num_ext_allowed = num_ext;
        e_PCR->setShort(shortout);
        e_PCR->setHitlimit(hitlimit);
        e_PCR->setReversed(Re);
        e_PCR->setReference(ref);
        e_PCR->setfirstlast(first_read, last_read);
        e_PCR->setAverageQV(AQV);
        e_PCR->setXdrop(xdrop);
        e_PCR->setMismatch(mismatch);
	e_PCR->setOffset(offset);
	if (mismatch != 0) e_PCR->setPrefixLen(prefix_len);
        if (uniqueH) e_PCR->setUniqueCheck();
}


int epcr_main(char *title, int argc, char **argv, mapping_machine *e_PCR)
{
	//// Gather arguments

	const char *stsfile = NULL;
	const char *seqfile = NULL;
	int i;

	for (i=1; i<argc; ++i)
	{
		if (argv[i][1] == '=')         // X=value
		{
			if (argv[i][2] == 0)
				fprintf(stderr,"Missing value for %s\n",argv[i]);
			else if (argv[i][0] == 'T')
				threshold = atoi(argv[i]+2);
			else if (argv[i][0] == 'V') 
				AQV = (double) atof(argv[i]+2);
			else if (argv[i][0] == 'W')
				wdsize = atoi(argv[i]+2);
			else if (argv[i][0] == 'G') {
			    gap = atoi(argv[i]+2);
			} else if (argv[i][0] == 'X') {
			    mfile = argv[i]+2;
			}  else if (argv[i][0] == 'A') {
			    dis = atoi(argv[i]+2);
			} else if (argv[i][0] == 'Y') {
                            cc = argv[i]+2;
			} else if (argv[i][0] == 'd') {
			    dw = atoi(argv[i]+2);
			} else if (argv[i][0] == 'P') {
			    pat = argv[i]+2;
			} else if (argv[i][0] == 'L') {
			    len = atoi(argv[i]+2);
			} else if (argv[i][0] == 'E') {
			    tempout = argv[i]+2;
			} else if (argv[i][0] == 'F') {
			    filter = atoi(argv[i]+2);
			} else if (argv[i][0] == 'C') {
			    colorspace = atoi(argv[i]+2);
			/*
			} else if (argv[i][0] == 'D') {
			    dir = atoi(argv[i]+2);
			*/
			} else if (argv[i][0] == 'B') {
			    both = atoi(argv[i]+2);
			} else if (argv[i][0] == 'S') {
			    shortout = atoi(argv[i]+2);
			} else if (argv[i][0] == 'M') {
			    mode = atoi(argv[i]+2);
			} else if (argv[i][0] == 'Z') {
			    hitlimit = atoi(argv[i]+2);
			} else if (argv[i][0] == 'n') {
			    reportP = atoi(argv[i]+2);
			} else if (argv[i][0] == 'R') {
                            Re = atoi(argv[i]+2);
			} else if (argv[i][0] == 'r') {
                            ref = atoi(argv[i]+2);
                        } else if (argv[i][0] == 'Q') {
			    qfile = argv[i]+2;
			} else if (argv[i][0] == 'b') {
                            first_read= atoll(argv[i]+2);
                        } else if (argv[i][0] == 'e') {
                            last_read= atoll(argv[i]+2);
                        }  else if (argv[i][0] == 'U') {
			    SNPrate = atof(argv[i]+2);
			} else if (argv[i][0] == 'H') { 
			    amballow = atoi(argv[i]+2);
			} else if (argv[i][0] == 'u') {
			    uniqueH = atoi(argv[i]+2);
			} else if (argv[i][0] == 'x') {
			    xdrop = atoi(argv[i]+2);
                        }else if (argv[i][0] == 'm') {
			    mismatch = atoi(argv[i]+2);
			} else if (argv[i][0] == 's') {
			    prefix_len = atoi(argv[i]+2);
                        } else if (argv[i][0] == 'q') {
			    offset = atoll(argv[i]+2);
			} 
		} else if (argv[i][0] == 'n' && argv[i][1] == 'p' && argv[i][2] == '=') {
		    numproc = atoi(argv[i]+3);
		} else if (argv[i][0] == 'n' && argv[i][1] == 'x' && argv[i][2] == '=') {
		    num_ext = atoi(argv[i]+3);
                }
		else   // filename
		{
			if (stsfile == NULL)
				stsfile = argv[i];
			else if (seqfile ==NULL)
				seqfile = argv[i];
			else 
				fprintf(stderr,"Argument \"%s\" ignored\n",argv[i]);
		}
	}

	if (stsfile==NULL || seqfile==NULL) {
	    return Usage();
	}
	if (numproc > 1) {mode = 0;}
	buildmachine(e_PCR, NULL);
	if (numproc > 1) {
	    long long llen = buildlookup(e_PCR, seqfile, mode);
	    listmachine = new mapping_machine*[numproc];
	    char *qfileend = NULL;
	    char *tp = tempout, file[1000]; 
	    tempout = file; /// set the tempfile pointer to this new bariable file name
	    for (i = 0; i < numproc; i++) {
		if (qfile) {
		    qfileend = strchr(qfile, ':');
		    if (qfileend) {*qfileend = 0;}
		}
		sprintf(file, "%s_%d", tp, i); // each machine has its own tempfile
		mapping_machine *mp = new mapping_machine();
		buildmachine(mp, e_PCR);
		sprintf(file, "%s_%d_tmpout", tp, i);
		mp->setOutputfile(file);
		listmachine[i] = mp;
		if (qfileend) qfile = qfileend+1;
		else qfile = NULL;
	    } 
	    return !multi_thread(stsfile, seqfile, mode, numproc, llen);
	} 
	return !run_pcr(e_PCR, stsfile, seqfile, mode);
}

int main (int argc, char **argv)
{
	mapping_machine *probem = new mapping_machine();
	
	return epcr_main("pbm", argc, argv,probem);
}
