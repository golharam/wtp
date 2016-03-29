SRCDIR = src
CFLAGS = -g -O2 -march=nocona -pipe -ffast-math
BIN=.
TARGETS := $(BIN)/map $(BIN)/mapreads $(BIN)/remduphits $(BIN)/radsort \
	$(BIN)/panelGrouping \
	$(BIN)/pairing $(BIN)/seqstat $(BIN)/makepattern $(BIN)/rescue \
	$(BIN)/simu $(BIN)/simu_del $(BIN)/miniassem $(BIN)/mutseq \
	$(BIN)/lowcov $(BIN)/getreads $(BIN)/checkpattern $(BIN)/remap \
	$(BIN)/denovoadp $(BIN)/mpindel_summ $(BIN)/mpindel_parse \
	$(BIN)/comp_dbSNP $(BIN)/merge_map $(BIN)/add_reads $(BIN)/merge_pat \
	$(BIN)/pairing_barcode $(BIN)/prepairing $(BIN)/classify $(BIN)/buildschemaDB \
	$(BIN)/newSchema $(BIN)/schema2adj $(BIN)/schemalimit1

DBSCHEMA = schemas/DBschema

ifndef PREFIX
  PREFIX=/share/apps/corona
endif
ifndef DESTDIR
  DESTDIR=
endif

.PHONY: all install install_targets install_schemas install_lite dbschema $(DBSCHEMA)

all : $(TARGETS) dbschema

dbschema: $(DBSCHEMA) 

$(DBSCHEMA) : $(BIN)/buildschemaDB 
	cd schemas && (ls schema_[0-9]* > lst)
	cd schemas && ../buildschemaDB lst > DBschema

install: install_targets install_schemas

install_lite: install

install_targets: $(TARGETS)
	install -d -m 775 ${DESTDIR}${PREFIX}/bin
	install -m 775 $? ${DESTDIR}${PREFIX}/bin/

#install -m 664 templates ${DESTDIR}${PREFIX}/conf/templates
#install -m 664 templates ${DESTDIR}${PREFIX}/bin/templates


install_schemas: $(wildcard schemas/*schema_*) $(DBSCHEMA)
	install -d -m 775 ${DESTDIR}${PREFIX}/etc/schemas/
	install -m 444 $? ${DESTDIR}${PREFIX}/etc/schemas/
	install -m 444 $(DBSCHEMA) ${DESTDIR}${PREFIX}/etc/schemas

clean: 
	-rm -f *.o
	-rm -f $(TARGETS)


clean_lite:  
	-rm -f *.o
	-rm -f $(TARGETS)

LFS_CFLAGS    := $(shell getconf LFS_CFLAGS 2>/dev/null) -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE
LFS_LDFLAGS   := $(shell getconf LFS_LDFLAGS 2>/dev/null)
LFS_LIBS      := $(shell getconf LFS_LIBS 2>/dev/null)

# _LARGEFILE_SOURCE gets you O_LARGEFILE for open, and fseeko and ftello
# interfaces.
# _LARGEFILE64_SOURCE gets you the xxx64() versions of calls that take
# off64_t, etc., arguments.
# _FILE_OFFSET_BITS=64 makes off_t, etc., be 64 bits.

CFLAGS        += $(LFS_CFLAGS)
LDFLAGS       = $(LFS_LDFLAGS)
LDLIBS        = $(LFS_LIBS)

CXXFLAGS=$(CFLAGS)

$(BIN)/add_reads: $(SRCDIR)/add_reads.cxx zutil.o util.o
	$(CXX) $(CFLAGS) $(SRCDIR)/add_reads.cxx zutil.o util.o -o $(BIN)/add_reads

$(BIN)/merge_map: $(SRCDIR)/merge_map.cxx zutil.o util.o zcompress.o
	$(CXX) $(CFLAGS) $(SRCDIR)/merge_map.cxx zutil.o util.o zcompress.o -o $(BIN)/merge_map -lz

$(BIN)/buildschemaDB: $(SRCDIR)/buildschemaDB.cxx zutil.o util.o
	$(CXX) $(CFLAGS) $(SRCDIR)/buildschemaDB.cxx zutil.o util.o -o $(BIN)/buildschemaDB

$(BIN)/newSchema: $(SRCDIR)/newSchema.cxx zutil.o util.o
	$(CXX) $(CFLAGS) $(SRCDIR)/newSchema.cxx zutil.o util.o -o $(BIN)/newSchema

$(BIN)/merge_pat : $(SRCDIR)/merge_pat.cxx zutil.o util.o
	$(CXX) $(CFLAGS) $(SRCDIR)/merge_pat.cxx zutil.o util.o -o $(BIN)/merge_pat

$(BIN)/comp_dbSNP: $(SRCDIR)/comp_dbSNP.cxx
	$(CXX) $(CFLAGS) $(SRCDIR)/comp_dbSNP.cxx -o $(BIN)/comp_dbSNP

$(BIN)/denovoadp: $(SRCDIR)/denovoadp.cxx zutil.o util.o
	$(CXX) $(CFLAGS) $(SRCDIR)/denovoadp.cxx zutil.o util.o -o $(BIN)/denovoadp

$(BIN)/mpindel_summ: $(SRCDIR)/mpindel_summ.cxx
	$(CXX) $(CFLAGS) $(SRCDIR)/mpindel_summ.cxx -o $(BIN)/mpindel_summ

$(BIN)/mpindel_parse: $(SRCDIR)/mpindel_parse.cxx zcompress.o
	 $(CXX) $(CFLAGS) $(SRCDIR)/mpindel_parse.cxx zcompress.o -o $(BIN)/mpindel_parse -lz

$(BIN)/makepattern: $(SRCDIR)/makepattern.cxx zutil.o util.o
	$(CXX) $(CFLAGS) $(SRCDIR)/makepattern.cxx zutil.o util.o -o $(BIN)/makepattern

$(BIN)/seqstat: $(SRCDIR)/seqstat.cxx zutil.o util.o
	$(CXX) $(CFLAGS) $(SRCDIR)/seqstat.cxx zutil.o util.o -o $(BIN)/seqstat

$(BIN)/mutseq: $(SRCDIR)/mutseq.cxx zutil.o util.o fasta-io.o
	$(CXX) $(CFLAGS) $(SRCDIR)/mutseq.cxx zutil.o util.o fasta-io.o -o $(BIN)/mutseq

$(BIN)/lowcov: $(SRCDIR)/lowcov.cxx zutil.o util.o
	$(CXX) $(CFLAGS) $(SRCDIR)/lowcov.cxx zutil.o util.o -o $(BIN)/lowcov

$(BIN)/getreads: $(SRCDIR)/getreads.cxx zutil.o util.o
	$(CXX) $(CFLAGS) $(SRCDIR)/getreads.cxx zutil.o util.o -o $(BIN)/getreads

$(BIN)/panelGrouping: $(SRCDIR)/panelGrouping.cxx zcompress.o 
	$(CXX) $(CFLAGS) -Wall $(SRCDIR)/panelGrouping.cxx zcompress.o -o  $(BIN)/panelGrouping -lz

$(BIN)/pairing:  $(SRCDIR)/pairing.cxx pairTags.o zutil.o getComp.o util.o rescue.o map.o fasta-io.o zcompress.o
	$(CXX) $(CFLAGS) $(SRCDIR)/pairing.cxx  pairTags.o zutil.o getComp.o util.o rescue.o map.o fasta-io.o zcompress.o -o $(BIN)/pairing -lz

$(BIN)/pairing_barcode: $(SRCDIR)/pairing_barcode.cxx zutil.o util.o fasta-io.o
	 $(CXX) $(CFLAGS) $(SRCDIR)/pairing_barcode.cxx zutil.o util.o fasta-io.o -o $(BIN)/pairing_barcode

$(BIN)/prepairing: $(SRCDIR)/prepairing.cxx zutil.o util.o fasta-io.o
	 $(CXX) $(CFLAGS) $(SRCDIR)/prepairing.cxx zutil.o util.o fasta-io.o -o $(BIN)/prepairing

$(BIN)/rescue: $(SRCDIR)/rescuemain.cxx getComp.o rescue.o fasta-io.o zutil.o util.o map.o zcompress.o
	$(CXX) $(CFLAGS) $(SRCDIR)/rescuemain.cxx rescue.o getComp.o fasta-io.o zutil.o util.o map.o zcompress.o -lz -o $(BIN)/rescue

$(BIN)/remap : rescue.o getComp.o fasta-io.o zutil.o util.o map.o  zcompress.o $(SRCDIR)/remap.cxx
	$(CXX) $(CFLAGS) $(SRCDIR)/remap.cxx getComp.o rescue.o fasta-io.o zutil.o util.o map.o zcompress.o -lz -o $(BIN)/remap

$(BIN)/radsort: $(SRCDIR)/radsort.cxx zutil.o util.o
	$(CXX) $(CFLAGS) $(SRCDIR)/radsort.cxx zutil.o util.o -o $(BIN)/radsort

$(BIN)/mapreads: $(SRCDIR)/mapreads.cxx zutil.o util.o fasta-io.o zcompress.o
	$(CXX) $(CFLAGS) $(SRCDIR)/mapreads.cxx zutil.o util.o fasta-io.o zcompress.o -o $(BIN)/mapreads -lz

$(BIN)/mapreads-mp: $(SRCDIR)/mapreads-mp.cxx zutil.o util.o fasta-io.o zcompress.o
	$(CXX) $(CFLAGS) $(SRCDIR)/mapreads.cxx zutil.o util.o fasta-io.o zcompress.o -o $(BIN)/mapreads-mp -lz

$(BIN)/checkpattern: $(SRCDIR)/checkpattern.cxx zutil.o
	$(CXX) $(CFLAGS) $(SRCDIR)/checkpattern.cxx zutil.o -o $(BIN)/checkpattern

$(BIN)/schema2adj: $(SRCDIR)/schema2adj.cxx zutil.o
	$(CXX) $(CFLAGS) $(SRCDIR)/schema2adj.cxx zutil.o -o $(BIN)/schema2adj

$(BIN)/schemalimit1: $(SRCDIR)/schemalimit1.cxx zutil.o
	 $(CXX) $(CFLAGS) $(SRCDIR)/schemalimit1.cxx zutil.o -o $(BIN)/schemalimit1

$(BIN)/map: $(SRCDIR)/map.h map.o zutil.o util.o $(SRCDIR)/mapmain.cxx fasta-io.o zcompress.o
	$(CXX) $(CFLAGS) $(SRCDIR)/mapmain.cxx map.o -o $(BIN)/map util.o zutil.o fasta-io.o zcompress.o -lz -lpthread

$(BIN)/map-mp: $(SRCDIR)/map.h map.o zutil.o util.o $(SRCDIR)/mapmain-mp.cxx fasta-io.o zcompress.o
	$(CXX) $(CFLAGS) $(SRCDIR)/mapmain-mp.cxx map.o -o $(BIN)/map-mp util.o zutil.o fasta-io.o zcompress.o -lz

$(BIN)/classify: $(SRCDIR)/map.h map.o zutil.o util.o $(SRCDIR)/classify.cxx fasta-io.o zcompress.o
	$(CXX) $(CFLAGS) $(SRCDIR)/classify.cxx map.o -o $(BIN)/classify util.o zutil.o fasta-io.o zcompress.o -lz

$(BIN)/remduphits: $(SRCDIR)/remduphits.cxx zutil.o zcompress.o
	$(CXX) $(CFLAGS) $(SRCDIR)/remduphits.cxx -o $(BIN)/remduphits zutil.o zcompress.o -lz

$(BIN)/simu: simu.o fasta-io.o util.o zutil.o $(SRCDIR)/simu_main.cxx 
	$(CXX) $(CFLAGS) simu.o $(SRCDIR)/simu_main.cxx util.o fasta-io.o zutil.o -o $(BIN)/simu

$(BIN)/simu_del: simu.o fasta-io.o util.o zutil.o $(SRCDIR)/FNdel.cxx
	$(CXX) $(CFLAGS) simu.o $(SRCDIR)/FNdel.cxx util.o fasta-io.o zutil.o -o $(BIN)/simu_del

$(BIN)/miniassem: miniassem.o zutil.o util.o $(SRCDIR)/mytemplate.h fasta-io.o $(SRCDIR)/minimain.cxx
	$(CXX) $(CFLAGS) miniassem.o zutil.o util.o fasta-io.o $(SRCDIR)/minimain.cxx -o $(BIN)/miniassem

simu.o: $(SRCDIR)/simu.cxx $(SRCDIR)/util.h $(SRCDIR)/zutil.h
	$(CXX) $(CFLAGS) -c $(SRCDIR)/simu.cxx

map.o: $(SRCDIR)/map.h $(SRCDIR)/map.cxx $(SRCDIR)/zcompress.h $(SRCDIR)/util.h $(SRCDIR)/zutil.h
	$(CXX) $(CFLAGS) -c $(SRCDIR)/map.cxx

getComp.o: $(SRCDIR)/getComp.cxx $(SRCDIR)/getComp.h
	$(CXX) $(CFLAGS) -c $(SRCDIR)/getComp.cxx

fasta-io.o : $(SRCDIR)/fasta-io.cxx $(SRCDIR)/fasta-io.h $(SRCDIR)/util.h $(SRCDIR)/zutil.h
	$(CXX) $(CFLAGS) -c $(SRCDIR)/fasta-io.cxx

util.o : $(SRCDIR)/util.cxx $(SRCDIR)/util.h
	$(CXX) $(CFLAGS) -c $(SRCDIR)/util.cxx

zutil.o: $(SRCDIR)/zutil.cxx $(SRCDIR)/zutil.h
	$(CXX) $(CFLAGS) -c $(SRCDIR)/zutil.cxx

rescue.o: $(SRCDIR)/rescue.cxx $(SRCDIR)/rescue.h $(SRCDIR)/map.h $(SRCDIR)/util.h $(SRCDIR)/zutil.h
	$(CXX) $(CFLAGS) -c $(SRCDIR)/rescue.cxx

miniassem.o: $(SRCDIR)/miniassem.cxx $(SRCDIR)/miniassem.h $(SRCDIR)/util.h $(SRCDIR)/zutil.h
	$(CXX) $(CFLAGS) -c $(SRCDIR)/miniassem.cxx

zcompress.o: $(SRCDIR)/zcompress.cxx $(SRCDIR)/zcompress.h $(SRCDIR)/util.h $(SRCDIR)/zutil.h
	$(CXX) $(CFLAGS) -Wall -c $(SRCDIR)/zcompress.cxx

pairTags.o: $(SRCDIR)/pairTags.cxx $(SRCDIR)/pairTags.h
	$(CXX) $(CFLAGS) -Wall -c $(SRCDIR)/pairTags.cxx
