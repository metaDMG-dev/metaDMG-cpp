# --- Flags og kompileringsvalg ---
FLAGS     := -O2
CFLAGS    := $(FLAGS)
CXXFLAGS  := $(FLAGS)
CPPFLAGS  := $(CPPFLAGS) -Wall -Wextra
LDFLAGS   := -lgsl

CC  ?= gcc
CXX ?= g++

CXXSRC := $(wildcard *.cpp)
CSRC   := $(wildcard *.c)
OBJ    := $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)

# --- Brugervalgte paths og biblioteker ---
#EIGEN_INC := /opt/homebrew/Cellar/eigen/3.4.0_1/include
#CPPFLAGS  := $(CPPFLAGS) -I$(EIGEN_INC)

ifdef PREFIX
CPPFLAGS += -I$(PREFIX)/include
LDFLAGS  += -L$(PREFIX)/lib
endif

ifdef EXTRA_INC
CPPFLAGS += $(addprefix -I,$(EXTRA_INC))
endif

ifdef EXTRA_LIB
LDFLAGS  += $(addprefix -L,$(EXTRA_LIB))
endif

ifdef EXTRA_LIBS
LIBS += $(EXTRA_LIBS)
endif

# --- Crypto library detektion ---
HAVE_CRYPTO := $(shell echo 'int main(){}'|$(CXX) -x c++ - -lcrypto -o /dev/null 2>/dev/null && echo 0 || echo 1)
LIBS := -lz -lm -lbz2 -llzma -lpthread -lcurl

# --- Mål og bygning ---
PROGRAM = metaDMG-cpp
all: version.h $(PROGRAM) misc

ifeq ($(HAVE_CRYPTO),0)
  $(info Crypto library is available to link; adding -lcrypto)
  LIBS += -lcrypto
else
  $(info Crypto library is not available to link; skipping -lcrypto)
endif

# --- HTSLIB håndtering ---
ifndef HTSSRC
  $(info HTSSRC not defined; cloning htslib from GitHub)
  HTSSRC := $(CURDIR)/htslib
endif

ifeq ($(HTSSRC),systemwide)
  $(info HTSSRC set to systemwide; using systemwide installation)
  LIBS += -lhts
  LIBHTS :=
else
  # Use HTSSRC directly for include path
  CPPFLAGS += -I$(HTSSRC)
  LIBHTS := $(HTSSRC)/libhts.a
  LIBS := $(LIBHTS) $(LIBS)

  $(PROGRAM): $(LIBHTS)
  $(LIBHTS): .clone_htslib
endif

.PHONY: .clone_htslib
.clone_htslib:
	@if [ ! -d "$(HTSSRC)" ]; then \
		echo "Cloning htslib into $(HTSSRC)..."; \
		git clone --recursive https://github.com/samtools/htslib.git $(HTSSRC) || { echo "Clone failed"; exit 1; }; \
	else \
		echo "$(HTSSRC) already exists, skipping clone."; \
	fi
	@$(MAKE) -C $(HTSSRC) libhts.a || { echo "Failed to build libhts.a"; exit 1; }
	$(info CPPFLAGS after htslib setup: $(CPPFLAGS) HTSSRC=$(HTSSRC))

# --- Versionsnummer og version.h ---
PACKAGE_VERSION := 0.4
ifneq ("$(wildcard .git)","")
  PACKAGE_VERSION := $(shell git describe --always --dirty)
endif

version.h:
	@echo '#define METADAMAGE_VERSION "$(PACKAGE_VERSION)"' > version.h.tmp
	@cmp -s version.h.tmp version.h || mv version.h.tmp version.h
	@rm -f version.h.tmp

$(PROGRAM): version.h $(OBJ) $(LIBHTS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(OBJ) $(LIBS) $(LDFLAGS)

.PHONY: misc
misc:
	$(MAKE) -C misc HTSSRC="$(HTSSRC)"

# --- Automatisk afhængighedshåndtering ---
-include $(OBJ:.o=.d)

%.o: %.c .clone_htslib
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(CC) -MM $(CPPFLAGS) $(CFLAGS) $< > $*.d

%.o: %.cpp .clone_htslib
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $< > $*.d

# --- Rens og tests ---
.PHONY: clean test testclean force

clean: testclean
	rm -f *.o *.d $(PROGRAM) version.h *~
	rm -rf $(HTSSRC)
	$(MAKE) -C misc clean

testclean:
	rm -f test/*.bam.bin test/data/*.bam
	rm -rf test/output test/logfile version.h

test:
	@echo "Running unit tests for metaDMG"
	cd test ; ./testAll.sh

force:
