# --- Flags og kompileringsvalg ---
FLAGS     ?= -O2 -Wall -Wextra
CFLAGS    += $(FLAGS) -MMD -MP
CXXFLAGS  += $(FLAGS) -MMD -MP
# CPPFLAGS  +=  <- dont think we need this

LDFLAGS  +=
LDLIBS   += -lgsl -lgslcblas -lz -lm -lpthread
LDHTS    := -lbz2 -llzma -lcurl

CC  ?= gcc
CXX ?= g++

CXXSRC := $(wildcard *.cpp)
CSRC   := $(wildcard *.c)
OBJ    := $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)


#$(info PREFIX before = '$(PREFIX)')
ifneq ($(strip $(PREFIX)),)
  override PREFIX := $(abspath $(PREFIX))
endif
#$(info PREFIX after  = '$(PREFIX)')

ifneq ($(strip $(EXTRA_INC)),)
  override EXTRA_INC := $(foreach d,$(EXTRA_INC),$(abspath $(d)))
endif

ifneq ($(strip $(EXTRA_LIB)),)
  override EXTRA_LIB := $(foreach d,$(EXTRA_LIB),$(abspath $(d)))
endif

ifneq ($(strip $(PREFIX)),)
  CPPFLAGS += -I$(PREFIX)/include
  LDFLAGS  += -L$(PREFIX)/lib
endif

ifneq ($(strip $(EXTRA_INC)),)
  CPPFLAGS += $(addprefix -I,$(EXTRA_INC))
endif

ifneq ($(strip $(EXTRA_LIB)),)
  LDFLAGS  += $(addprefix -L,$(EXTRA_LIB))
endif

ifneq ($(strip $(EXTRA_LIBS)),)
  LDLIBS += $(EXTRA_LIBS)
endif



# --- Crypto library detektion ---
HAVE_CRYPTO := $(shell echo 'int main(){}'|$(CXX) -x c++ - -lcrypto -o /dev/null 2>/dev/null && echo 0 || echo 1)


# --- Mål og bygning ---
PROGRAM = metaDMG-cpp
all: version.h $(PROGRAM) misc

ifeq ($(HAVE_CRYPTO),0)
  $(info Crypto library is available to link; adding -lcrypto)
  LDLIBS += -lcrypto
else
  $(info Crypto library is not available to link; skipping -lcrypto)
endif

# --- HTSLIB håndtering ---
# HTSSRC is an absolute path (e.g., $(CURDIR)/htslib or user-specified)

ifndef HTSSRC
  $(info HTSSRC not defined; cloning htslib from GitHub)
  HTSSRC := $(CURDIR)/htslib
  ABSPATH=$(HTSSRC) #donkykong
endif

ABSPATH=$(HTSSRC) #donkykong
ifeq ($(HTSSRC),systemwide)
  $(info HTSSRC set to systemwide; using systemwide installation)
  LDLIBS += -lhts
  LIBHTS :=
else
  # Use HTSSRC directly for include path
  CPPFLAGS += -I$(HTSSRC)
  LIBHTS := $(HTSSRC)/libhts.a
  LDLIBS := $(LIBHTS) $(LIBS) $(LDHTS) $(LDLIBS)
  $(PROGRAM): $(LIBHTS)

  ifneq ($(filter /%,$(HTSSRC)),$(HTSSRC))
    ABSPATH=../$(HTSSRC)
  endif
endif

# Ensure htslib is cloned and built only if libhts.a is missing
$(LIBHTS): .clone_htslib

.clone_htslib:
	@if [ ! -d "$(HTSSRC)" ]; then \
		echo "Cloning htslib into $(HTSSRC) with submodules..."; \
		git clone --recursive https://github.com/samtools/htslib.git $(HTSSRC) || { echo "Clone failed"; exit 1; }; \
	else \
		echo "$(HTSSRC) already exists, skipping clone."; \
	fi
	@if [ ! -f "$(LIBHTS)" ]; then \
		$(MAKE) -C $(HTSSRC) libhts.a || { echo "Failed to build libhts.a"; exit 1; }; \
	else \
		echo "$(LIBHTS) already exists, skipping build."; \
	fi
	$(info CPPFLAGS=$(CPPFLAGS) HTSSRC=$(HTSSRC) (absolute path))

# --- Versionsnummer og version.h ---
PACKAGE_VERSION := unknown
ifneq ("$(wildcard .git)","")
  PACKAGE_VERSION := $(shell git describe --tags --always --dirty)
endif

version.h:
	@echo '#define METADAMAGE_VERSION "$(PACKAGE_VERSION)"' > version.h.tmp
	@cmp -s version.h.tmp version.h || mv version.h.tmp version.h
	@rm -f version.h.tmp

$(PROGRAM): version.h $(OBJ) $(LIBHTS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $(OBJ) $(LDFLAGS) $(LDLIBS)

.PHONY: misc
misc: $(LIBHTS) $(OBJ)
	$(MAKE) -C misc HTSSRC=$(ABSPATH) PREFIX="$(PREFIX)"

# --- Automatisk afhængighedshåndtering ---
-include $(OBJ:.o=.d)

bfgs.o: bfgs.cpp $(LIBHTS)
	$(CXX) -c $(CPPFLAGS) $(filter-out -Wall -Wextra,$(CXXFLAGS)) -w $< -o $@

%.o: %.c $(LIBHTS)
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@

%.o: %.cpp $(LIBHTS)
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@


# --- Rens og tests ---
.PHONY: clean test unittest regtest testclean force

.PHONY: print
print:
	@echo "---- COMPILER ----"
	@echo "CC       = $(CC)"
	@echo "CXX      = $(CXX)"
	@echo ""
	@echo "---- FLAGS ----"
	@echo "CFLAGS   = $(CFLAGS)"
	@echo "CXXFLAGS = $(CXXFLAGS)"
	@echo "CPPFLAGS = $(CPPFLAGS)"
	@echo "LDFLAGS  = $(LDFLAGS)"
	@echo "LDLIBS   = $(LDLIBS)"
	@echo ""
	@echo "---- ORIGIN ----"
	@echo "CFLAGS origin   = $(origin CFLAGS)"
	@echo "CXXFLAGS origin = $(origin CXXFLAGS)"
	@echo "CPPFLAGS origin = $(origin CPPFLAGS)"
	@echo "LDFLAGS origin  = $(origin LDFLAGS)"
	@echo "LDLIBS origin   = $(origin LDLIBS)"

clean: 
	rm -f *.o *.d $(PROGRAM) version.h *~
ifneq ($(HTSSRC),systemwide)
	rm -rf $(HTSSRC)
endif
	$(MAKE) -C misc clean

testclean:
	rm -f test/*.bam.bin test/*.log test/testAll.sh.log test/testAll2.sh.log
	rm -rf test/bin
	rm -f test/test_softclips test/test_cigar_ops test/test_dfit_bounds
	rm -f test/*.bam
	rm -rf test/output_extract_reads
	rm -f test/data/f570b1db7c.dedup.filtered.rname.bam
	rm -rf test/output test/output_data2 test/output_data2_gd test/output_compressbam

misc/sort_taxonomy_output:
	$(MAKE) -C misc sort_taxonomy_output

# --- Unit tests: small C++ programs, no test data or samtools needed ---
UNITTEST_BIN := test/bin/test_dust_score test/bin/test_pmd_model \
                test/bin/test_softclips test/bin/test_cigar_ops test/bin/test_dfit_bounds

# objects each unit test links against, besides htslib and LDLIBS
test/bin/test_dust_score:  shared.o
test/bin/test_softclips:   profile.o
test/bin/test_cigar_ops:   profile.o
test/bin/test_dfit_bounds: dfit_optim.o bfgs.o mrand.o
# test_pmd_model includes main_pmd.cpp directly; let the linker drop unused code
ifeq ($(shell uname -s),Darwin)
  UNIT_GC := -Wl,-dead_strip
else
  UNIT_GC := -Wl,--gc-sections
endif
test/bin/test_pmd_model: UNIT_EXTRA := -ffunction-sections -fdata-sections $(UNIT_GC)

test/bin/%: test/%.cpp version.h $(LIBHTS)
	@mkdir -p test/bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(UNIT_EXTRA) -o $@ $< $(filter %.o,$^) $(LDFLAGS) $(LDLIBS)

-include $(UNITTEST_BIN:=.d)

unittest: $(UNITTEST_BIN) misc/sort_taxonomy_output
	@echo "Running unit tests for metaDMG"
	@set -e; for t in $(UNITTEST_BIN); do echo "== $$t"; ./$$t; done
	./misc/test_sort_taxonomy_output.sh

# --- Regression tests: run metaDMG-cpp on test data and compare output ---
regtest: $(PROGRAM) misc
	@echo "Running regression tests for metaDMG"
	cd test ; ./testAll2.sh

test: unittest regtest

force:
