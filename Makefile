#C_OPT=-O3 -funroll-loops -march=pentiumpro
#C_OPT=-O3 -march=pentiumpro
#C_DEBUG=-W -Wall -g
#C_DEF=-DUSE_OBSTACK

CFLAGS= -O -g
CFLAGS_interactiv=$(CFLAGS) -DTRACK_SEARCH
CC=gcc

BZ_SRCS=bz_main.c bz_align.c bz_extend.c bz_chain.c bz_dna.c bz_protein.c bz_print.c \
	bz_table.c bz_census.c bz_hit19.c bz_inner.c util.c seq.c args.c \
	edit.c dna.c charvec.c nib.c bz_site.c bz_iupac.c \
	recognize.c bz_weight.c det-marginal-prob.c seedlist.c

BZ_HEADERS=args.h bz_dna.h bz_protein.h bz_print.h dna.h nib.h bz_all.h bz_hit19.h \
	bz_table.h edit.h seq.h bz_census.h bz_main.h charvec.h \
	int_types.h util.h bz_site.h bz_iupac.h bz_iupac.decode recognize.h astack.c \
	bz_weight.h det-marginal-prob.h seedlist.h

SR_SRCS=strip_rpts.c util.c seq.c charvec.c nib.c
RC_SRCS=revcomp.c util.c seq.c charvec.c nib.c
RR_SRCS=restore_rpts.c util.c seq.c charvec.c nib.c dna.c
RT_SRCS=repeats_tag.c util.c

#all: blastz iblastz

blastz: $(BZ_SRCS) $(BZ_HEADERS)
	$(CC) $(CFLAGS) $(BZ_SRCS) \
		$(M) \
		-lm \
		-o $@

iblastz: $(BZ_SRCS) $(BZ_HEADERS)
	$(CC) $(CFLAGS_interactiv) $(BZ_SRCS) \
		$(M) \
		-lm \
		-o $@

iblastz2: $(BZ_SRCS) $(BZ_HEADERS)
	$(CC) $(CFLAGS_interactiv) $(BZ_SRCS) \
		$(M) \
		-lm \
		-o $@

createIUPACmatrixes: createIUPACmatrixes.c
	$(CC)  createIUPACmatrixes.c -o $@

bz_iupac.decode: createIUPACmatrixes
	createIUPACmatrixes > $@

strip_rpts : $(SR_SRCS)
	$(CC) $(CFLAGS) $(SR_SRCS) -o $@

revcomp : $(RC_SRCS)
	$(CC) $(CFLAGS) $(RC_SRCS) -o $@

restore_rpts : $(RR_SRCS)
	$(CC) $(CFLAGS) $(RR_SRCS) -o $@

repeats_tag : $(RT_SRCS)
	$(CC) $(CFLAGS) $(RT_SRCS) -o $@

# This test determines the additional regions aligned, for one very small
# dataset, using an approach to human-mouse alignments pioneered by Arian Smit.
# We start with human (softmasked) and mouse (softmasked) and the corresponding
# RepeatMasker ".out" files, human.out and mouse.out. Lineage-specific are
# removed, the sequences are aligned, and the blastz output files are adjusted
# as if run on the original sequences.

test : blastz strip_rpts revcomp restore_rpts repeats_tag
	primate.pl human.out
	./repeats_tag human.out.primspec > hum.new.rpts
	rodent.pl mouse.out
	./repeats_tag mouse.out.rodspec > mus.new.rpts
	./strip_rpts human hum.new.rpts > hum.stripped
	./strip_rpts mouse mus.new.rpts > mus.stripped
	./blastz hum.stripped mus.stripped B=0 | \
		./restore_rpts - \
		human "\"human\" 1 10001" \
		mouse "\"mouse\" 1 10000" \
		hum.new.rpts mus.new.rpts \
		> restored.blastz
	./revcomp mus.stripped > mus.stripped.rev
	./blastz hum.stripped mus.stripped.rev B=0 | \
		./restore_rpts - \
		human "\"human\" 1 10001" \
		mouse "\"mouse-\" 1 10000" \
		hum.new.rpts mus.new.rpts \
		reverse >> restored.blastz

	./blastz human mouse > regular.blastz
	wc regular.blastz restored.blastz

README : README.tex
	latex README.tex
	dvips -Ppdf -G0 -o README.ps README.dvi
	ps2pdf README.ps README.pdf


tarfile : README
	mkdir blastz-source
	cp README.pdf Makefile *.c *.h human human.out mouse mouse.out \
		primate.pl rodent.pl blastz-source
	tar -cf - blastz-source | gzip > blastz.tar.gz
	rm -rf blastz-source

dist : tarfile
	cp blastz.tar.gz /home/nog/httpd/html/dist








