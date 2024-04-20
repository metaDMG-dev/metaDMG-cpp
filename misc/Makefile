#assumed HTSSRC is a absolute path
CC  ?= gcc
CXX ?= g++

LIBS = -lz -lm -lbz2 -llzma -lpthread -lcurl

CRYPTO_TRY=$(shell echo 'int main(){}'|g++ -x c++ - -lcrypto 2>/dev/null; echo $$?)
ifeq "$(CRYPTO_TRY)" "0"
$(info Crypto library is available to link; adding -lcrypto to LIBS)
LIBS += -lcrypto
else
$(info Crypto library is not available to link; will not use -lcrypto)
endif


ifdef HTSSRC
$(info misc-HTSSRC defined $(HTSSRC))
CPPFLAGS += -I"$(HTSSRC)"
LIBS := "$(HTSSRC)/libhts.a" $(LIBS)
else
$(info HTSSRC not defined, assuming systemwide installation)
LIBS += -lhts
endif

FLAGS = $(CPPFLAGS) -O3 $(LDFLAGS)

.PHONY: all clean install

prefix      = /usr/local
exec_prefix = $(prefix)
bindir      = $(exec_prefix)/bin

INSTALL = install
INSTALL_DIR = $(INSTALL) -dm0755
INSTALL_PROGRAM = $(INSTALL) -Dm0755

PROGRAMS = compressbam extract_reads
PROGRAMS_MISC = MAP dfit_optim

all: $(PROGRAMS) $(PROGRAMS_MISC)

kmath.o: bfgs.cpp bfgs.h
	$(CC) $(FLAGS) -c bfgs.cpp

Likelihood.o: Likelihood.cpp
	$(CC) $(FLAGS) -c Likelihood.cpp

M3Load.o: M3Load.cpp
	$(CC) $(FLAGS) -c M3Load.cpp

compressbam: compressbam.cpp
	$(CXX) $(FLAGS) compressbam.cpp -o compressbam $(LIBS)

extract_reads: extract_reads.cpp ../shared.o
	$(CXX) $(FLAGS) extract_reads.cpp -o extract_reads ../shared.o $(LIBS)

MAP: MAP.cpp M3Load.o bfgs.o Likelihood.o
	$(CXX) $(FLAGS) -std=c++11 MAP.cpp M3Load.o bfgs.o Likelihood.o -o MAP $(LIBS) -D__WITH_MAIN__

dfit_optim: dfit_optim.cpp bfgs.o
	$(CXX) $(FLAGS) dfit_optim.cpp bfgs.o -o dfit_optim $(LIBS) -D__WITH_MAIN__
clean:
	rm -f *.o $(PROGRAMS) $(PROGRAMS_MISC) *~

install: $(PROGRAMS)
	$(INSTALL_DIR) $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) $(PROGRAMS) $(DESTDIR)$(bindir)

install-misc: $(PROGRAMS_MISC)
	$(INSTALL_DIR) $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) $(PROGRAMS_MISC) $(DESTDIR)$(bindir)
