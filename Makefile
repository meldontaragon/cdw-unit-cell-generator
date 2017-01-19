#!/usr/bin/make

#  Copyright (C) 2016-2017 David C. Miller

#  This file is part of TMD CDW Unit Cell Generator

#  TMD CDW Unit Cell Generatort is free software: you can redistribute it
#  and/or modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation, either version
#  3 of the License, or (at your option) any later version.

#  TMD CDW Unit Cell Generator is distributed in the hope that it will be
#  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.

#  You should have received a copy of the GNU Lesser General Public
#  License along with TMD CDW Unit Cell Generator.  If not, see
#  <http://www.gnu.org/licenses/>.

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
	$(CC) -c -o $@ $< $(FLAGS_BASE) $(OFLAG) $(CFLAGS)

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
debug: dir reset
	make all CFLAGS="$(CFLAGS_DEBUG)" OFLAG="$(OFLAG_DEBUG)"

final: dir reset
	make all CFLAGS="$(CFLAGS_FINAL)" OFLAG="$(OFLAG_FINAL)"

#install and make links in default bin directory
link: uninstall
	ln -s $(MAINDIR)/bin/* $(INSTALLDIR)/

install: final uninstall
	cp $(MAINDIR)/bin/* $(INSTALLDIR)/

uninstall:
	cd $(INSTALLDIR)/ && rm -f $(_BIN);
	cd $(MAINDIR);
