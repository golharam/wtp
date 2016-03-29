#include "fasta-io.h"
#include "zutil.h"
#include "miniassem.h"

main(int argc, char *argv[])
{
    if (argc < 2) {
	fatal("miniassem start end readsfile seq relen mini_over\n");
    }
    int relen = atoi(argv[5]);
    int start = atoi(argv[1]);
    int end = atoi(argv[2]);
    char *rfile = argv[3];
    char *seq = argv[4];
    int mini = atoi(argv[6]);
    
    FastaFile fafile(SEQTYPE_NT);

    if (!fafile.Open(seq,"r"))
        return 0;
    FastaSeq faseq;
    if (!fafile.Read(faseq)) return 0;
    fafile.Close();

    graph_assem *g = new graph_assem();

    g->init(100000, relen);
    g->inputReads(rfile);
    g->check_overlap(mini);
    g->assembly(faseq.Sequence(), start, end, mini);
}
