#include "map.h"
#include "fasta-io.h"

static long long beg = 0;
static long long end = ((long long) 1) << 60;

int Usage()
{
	fprintf(stderr,"\nclassify (0.1) utility program\n");

	fprintf(stderr,"USAGE:  classify readsfile seqfile outputRepeat [outputNoneRepeat] [options]\n\n");
	fprintf(stderr,"OPTIONS:\n");
	fprintf(stderr,"\tT=##     Threshold Score\n");
	fprintf(stderr, "\tL=##    lenght of reads(dafault 15)\n");
	fprintf(stderr, "\tB=beg\n\tE=end\n");
	//fprintf(stderr, "\tD=+/-1  Direction of the search for color space\n");

	fprintf(stderr,"\n");
	return 1;
}

int run_pcr(mapping_machine *e_PCR, const char *stsfile, const char *seqfile, const char *outR, const char *outN)
{
    FastaFile fafile(SEQTYPE_NT);

    if (!fafile.Open(seqfile,"r"))
	return 0;
    FastaSeq faseq;
	if (!fafile.Read(faseq))
	{
	    fatal("Sequence file is empty\n");
	}
	fafile.Close();

	e_PCR->run_classification(faseq.Sequence(), stsfile, outR, outN, beg, end);
}

int epcr_main(char *title, int argc, char **argv, mapping_machine *e_PCR)
{
	//// Gather arguments

	const char *stsfile = NULL;
	const char *seqfile = NULL;
	const char *outR = NULL, *outN = NULL;
	int threshold = 0;
	char *pat = NULL;
	char pat1[1000];
	int len = 15;
	int colorspace = 1;

	int i;

	for (i=1; i<argc; ++i)
	{
		if (argv[i][1] == '=')         // X=value
		{
			if (argv[i][2] == 0)
				fprintf(stderr,"Missing value for %s\n",argv[i]);
			else if (argv[i][0] == 'T')
				threshold = atoi(argv[i]+2);
                        else if (argv[i][0] == 'L')
                                len = atoi(argv[i]+2);
			else if (argv[i][0] == 'B')
				beg = atoll(argv[i]+2);
			else if (argv[i][0] == 'E')
                                end = atoll(argv[i]+2);
			else if  (argv[i][0] == 'P')
				pat = argv[i]+2;
		}
		else   // filename
		{
			if (stsfile == NULL)
				stsfile = argv[i];
			else if (seqfile ==NULL)
				seqfile = argv[i];
			else if (outR == NULL) 
				outR = argv[i];
                        else if (outN == NULL)
                                outN = argv[i];
			else 
				fprintf(stderr,"Argument \"%s\" ignored\n",argv[i]);
		}
	}

	if (outR==NULL) {
	    return Usage();
	}
	mapping_machine_color *mp = new mapping_machine_color();
	e_PCR = mp;
	len--;
	e_PCR->probe_init(len);
	if (pat == NULL) {
	    pat = pat1;
	    pat[len] = 0;
	    for (i = 0; i < len; i++) pat[i] = '0';
	    for (i = 0; i < 7; i++) {
	    	pat[i] = pat[len-i-1] = '1';
	    }
	} 
	e_PCR->SetWordSize(pat);
	e_PCR->SetThreshold(threshold);
	e_PCR->setColorspace(colorspace);

	return !run_pcr(e_PCR, stsfile, seqfile, outR, outN);
}

int main (int argc, char **argv)
{
	mapping_machine *probem = new mapping_machine();
	
	return epcr_main("pbm", argc, argv,probem);
}
