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
PACKAGE_VERSION := 0.4.5
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
.PHONY: clean test testclean force

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
	rm -f test/data/*.bam
	rm -rf test/output test/output_data2 test/output_data2_gd

test:
	@echo "Running unit tests for metaDMG"
	cd test ; ./testAll2.sh

force:
