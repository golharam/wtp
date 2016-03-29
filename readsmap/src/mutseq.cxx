#include "fasta-io.h"
#include "zutil.h"

static void output(int & cur, int pos, const char *seq, int len)
{
    int i;
    if (pos > len) pos = len; 
    for (i = cur; i < pos; i++) printf("%c", seq[i]);
    cur = pos; 
}

main(int argc, char *argv[])
{
    if (argc < 3) 
	fatal("mutseq seqfile mutfile\nMutate sequence according to operations in mutfile\n Position of the operations in mutfile must be sorted in increasing order\n");
    char *seqfile = argv[1];
    char *mfile = argv[2];
    char mfileout[1000];
    sprintf(mfileout, "%s.pos", mfile); 
    char line[100000];
    FILE *mfp = ckopen(mfile, "r");
    FILE *mfpo = ckopen(mfileout, "w");
    int x = 0;

    FastaFile fafile(SEQTYPE_NT);

    if (!fafile.Open(seqfile,"r"))
        return 0;
    FastaSeq faseq;
    if (!fafile.Read(faseq)) return 0;
    printf(">%s\n", faseq.Defline());
    int cur = 0, pos;
    char com, str[100000];
    while (fgets(line, sizeof line, mfp)) {
	sscanf(line, "%d %c %s", &pos, &com, str);
	output(cur, pos, faseq.Sequence(), faseq.Length());
	fprintf(mfpo, "%d %d\n", pos+x-5, pos+x+5);
	if (com =='I') {
	    printf("%s", str);
	    x += strlen(str);
	} else if (com == 'R') {
	    printf("%c", str[0]);
	    cur++;
	} else { //deletion
	    int dlen = atoi(str);
	    cur += dlen;
	    x -= dlen;
	}
    }
    output(cur, faseq.Length(), faseq.Sequence(), faseq.Length());
    printf("\n");
    fclose(mfp);
    fclose(mfpo);
}
