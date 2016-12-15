#!/usr/bin/make
#makefile
#defines all rules and should NOT be edited
#except as necessary

include makefile.include

#########################
#   Local Directories   #
#########################

#intentionally blank
DIRLIST =

#code directories
ODIR = obj
BDIR = bin
IDIR = inc
SDIR = src

DIRLIST += $(ODIR) $(BDIR) $(IDIR) $(SDIR)

CSUFFIX=.c

BINC=$(CC)
BINSUFFIX=$(CSUFFIX)
BINFLAGS=$(CFLAGS)

OFLAG=$(OFLAG_DEBUG)
#specific flags
CFLAGS=$(CFLAGS_DEBUG)

#working directory
MAINDIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

_OBJ = $(patsubst %.o,%$(CSUFFIX).o,$(_COBJ))

DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

#build
BIN = $(patsubst %,$(BDIR)/%,$(_BIN))

##########################################################################
# Do not change anything below here (unless you know what you're doing)  #
##########################################################################

.PRECIOUS: $(ODIR)/%$(CSUFFIX).o\

.PHONY: all
all: $(BIN)

#prog
$(ODIR)/%$(CSUFFIX).o: $(SDIR)/%$(CSUFFIX) $(DEPS) $(CLASS)
	$(CC) -c -o $@ $< $(FLAGS_BASE) $(OFLAG) $(LDFLAGS) $(CFLAGS)

$(BDIR)/%: $(OBJ) $(ODIR)/%$(BINSUFFIX).o
	$(BINC) $^ $(FLAGS_BASE) $(OFLAG) $(LDFLAGS) $(BINFLAGS) -o $@

#clean up
CLEAN_CMD = rm -f *~;
RESET_CMD =
CLEAN_CMD += rm -f $(ODIR)/* $(IDIR)/*~ $(SDIR)/*~;
RESET_CMD += rm -f $(BDIR)/*;

.PHONY: clean reset dir
clean:
	$(CLEAN_CMD)

reset: clean
	$(RESET_CMD)

dir:
	mkdir -p $(DIRLIST)

#force specific flags
.PHONY: debug final install uninstall link
debug: reset
	make all CFLAGS="$(CFLAGS_DEBUG)" OFLAG="$(OFLAG_DEBUG)"

final: reset
	make all CFLAGS="$(CFLAGS_FINAL)" OFLAG="$(OFLAG_FINAL)"

#install and make links in default bin directory
link: uninstall
	ln -s $(MAINDIR)/bin/* $(INSTALLDIR)/

install: final uninstall
	cp $(MAINDIR)/bin/* $(INSTALLDIR)/

uninstall:
	cd $(INSTALLDIR)/ && rm -f $(_BIN);
	cd $(MAINDIR);
