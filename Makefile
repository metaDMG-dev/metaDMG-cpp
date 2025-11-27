# --- Flags og kompileringsvalg ---
FLAGS     := -O3
CFLAGS    := $(FLAGS)
CXXFLAGS  := $(FLAGS)
CPPFLAGS  := $(CPPFLAGS) -Wall -Wextra
LDFLAGS   := -lgsl -lgslcblas

CC  ?= gcc
CXX ?= g++

CXXSRC := $(wildcard *.cpp)
CSRC   := $(wildcard *.c)
OBJ    := $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)



# --- Brugervalgte paths og biblioteker ---
ifdef PREFIX
CPPFLAGS += -I$(PREFIX)/include
LDFLAGS  := -L$(PREFIX)/lib $(LDFLAGS)
endif

ifdef EXTRA_INC
CPPFLAGS += $(addprefix -I,$(EXTRA_INC))
endif

ifdef EXTRA_LIB
LDFLAGS  := $(addprefix -L,$(EXTRA_LIB)) $(LDFLAGS)
endif

ifdef EXTRA_LIBS
LIBS += $(EXTRA_LIBS)
endif

# --- Crypto library detektion ---
HAVE_CRYPTO := $(shell echo 'int main(){}'|$(CXX) -x c++ - -lcrypto -o /dev/null 2>/dev/null && echo 0 || echo 1)
LIBS := -lz -lm -lbz2 -llzma -lpthread -lcurl -lhts

# --- Mål og bygning ---
PROGRAM = metaDMG-cpp
all: version.h $(PROGRAM) misc

ifeq ($(HAVE_CRYPTO),0)
  $(info Crypto library is available to link; adding -lcrypto)
  LIBS += -lcrypto
else
  $(info Crypto library is not available to link; skipping -lcrypto)
endif

# --- Versionsnummer og version.h ---
PACKAGE_VERSION := 0.4
ifneq ("$(wildcard .git)","")
  PACKAGE_VERSION := $(shell git describe --always --dirty)
endif

version.h:
	@echo '#define METADAMAGE_VERSION "$(PACKAGE_VERSION)"' > version.h.tmp
	@cmp -s version.h.tmp version.h || mv version.h.tmp version.h
	@rm -f version.h.tmp

$(PROGRAM): version.h $(OBJ)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(OBJ) $(LIBS) $(LDFLAGS)

.PHONY: misc
misc: $(OBJ)
	$(MAKE) -C misc

# --- Automatisk afhængighedshåndtering ---
-include $(OBJ:.o=.d)

%.o: %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(CC) -MM $(CPPFLAGS) $(CFLAGS) $< > $*.d

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $< > $*.d

# --- Rens og tests ---
.PHONY: clean test testclean force

clean: 
	rm -f *.o *.d $(PROGRAM) version.h *~
	$(MAKE) -C misc clean

testclean:
	rm -f test/*.bam.bin test/data/*.bam
	rm -rf test/output test/logfile version.h

test:
	@echo "Running unit tests for metaDMG"
	cd test ; ./testAll.sh

force:
