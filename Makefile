
DEST = pkg

all : 
	make -C readsmap
	cp readsmap/map $(DEST)
	cp readsmap/mapreads $(DEST)
	cp readsmap/remduphits $(DEST)
	chmod 755 $(DEST)/map
	chmod 755 $(DEST)/mapreads
	chmod 755 $(DEST)/remduphits
	chmod 755 $(DEST)/extendMappedReads.py
	chmod 755 bin/count_tags.pl
	chmod 755 bin/split_read_mapper.sh
	chmod 755 bin/ntr_finder.sh

clean :
	make -C readsmap clean
	rm $(DEST)/map
	rm $(DEST)/mapreads
	rm $(DEST)/remduphits
