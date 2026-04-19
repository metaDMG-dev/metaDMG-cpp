# Makefile for the misc directory

CC  ?= gcc
CXX ?= g++

FLAGS     ?= -O2 -Wall -Wextra
CFLAGS    += $(FLAGS) -MMD -MP
CXXFLAGS  += $(FLAGS) -MMD -MP

LDLIBS    += -lz -lm -lpthread
LDHTS     := -lbz2 -llzma -lcurl

# --- Brugervalgte paths og biblioteker ---

ifneq ($(strip $(PREFIX)),)
  override PREFIX := $(abspath $(PREFIX))
endif

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

# --- Crypto library detection ---
HAVE_CRYPTO := $(shell echo 'int main(){}' | $(CXX) -x c++ - -lcrypto -o /dev/null 2>/dev/null && echo 0 || echo 1)
ifeq ($(HAVE_CRYPTO),0)
  LDLIBS += -lcrypto
endif

# --- HTSLIB håndtering ---
ifeq ($(HTSSRC),systemwide)
  LDLIBS += -lhts
  LIBHTS :=
else
  ifneq ($(filter /%,$(HTSSRC)),$(HTSSRC))
    override HTSSRC := $(shell [ -d "$(HTSSRC)" ] && cd "$(HTSSRC)" && pwd || echo "__INVALID__")
    ifeq ($(HTSSRC),__INVALID__)
      $(error HTSSRC is not a valid directory in misc folder)
    endif
  endif
  CPPFLAGS += -I$(HTSSRC)
  LIBHTS := $(HTSSRC)/libhts.a
  LDLIBS := $(LIBHTS) $(LDLIBS) $(LDHTS)
endif

PROGRAMS := compressbam extract_reads
PROGRAMS_MISC := MAP dfit_optim
OBJ := kmath.o Likelihood.o M3Load.o

all: $(PROGRAMS) $(PROGRAMS_MISC)

$(PROGRAMS) $(PROGRAMS_MISC): $(LIBHTS)
$(OBJ): $(LIBHTS)

$(LIBHTS):
	@if [ ! -f "$(LIBHTS)" ]; then \
		echo "Error: $(LIBHTS) not found; ensure parent Makefile has built htslib"; \
		exit 1; \
	fi

../shared.o:
	@if [ ! -f "../shared.o" ]; then \
		echo "Error: ../shared.o not found; ensure parent Makefile has built it"; \
		exit 1; \
	fi

KMATH_CXXFLAGS := $(filter-out -Wall -Wextra,$(CXXFLAGS)) -w
kmath.o: bfgs.cpp bfgs.h
	$(CXX) -c $(CPPFLAGS) $(KMATH_CXXFLAGS) bfgs.cpp -o kmath.o

Likelihood.o: Likelihood.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) Likelihood.cpp -o Likelihood.o

M3Load.o: M3Load.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) M3Load.cpp -o M3Load.o

compressbam: compressbam.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) compressbam.cpp -o compressbam $(LDFLAGS) $(LDLIBS)

extract_reads: extract_reads.cpp ../shared.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) extract_reads.cpp ../shared.o -o extract_reads $(LDFLAGS) $(LDLIBS)

MAP: MAP.cpp M3Load.o kmath.o Likelihood.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -D__WITH_MAIN__ -std=c++11 MAP.cpp M3Load.o kmath.o Likelihood.o -o MAP $(LDFLAGS) $(LDLIBS)

dfit_optim: dfit_optim.cpp kmath.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -D__WITH_MAIN__ dfit_optim.cpp kmath.o -o dfit_optim $(LDFLAGS) $(LDLIBS)

# --- Dependencies ---
-include $(OBJ:.o=.d)

# --- Print variables ---
.PHONY: print
print:
	@echo "CC       = $(CC)"
	@echo "CXX      = $(CXX)"
	@echo "CFLAGS   = $(CFLAGS)"
	@echo "CXXFLAGS = $(CXXFLAGS)"
	@echo "CPPFLAGS = $(CPPFLAGS)"
	@echo "LDFLAGS  = $(LDFLAGS)"
	@echo "LDLIBS   = $(LDLIBS)"
	@echo "HTSSRC   = $(HTSSRC)"
	@echo "PREFIX   = $(PREFIX)"
	@echo "EXTRA_INC= $(EXTRA_INC)"
	@echo "EXTRA_LIB= $(EXTRA_LIB)"
	@echo ""
	@echo "---- ORIGIN ----"
	@echo "CC       origin = $(origin CC)"
	@echo "CXX      origin = $(origin CXX)"
	@echo "CFLAGS   origin = $(origin CFLAGS)"
	@echo "CXXFLAGS origin = $(origin CXXFLAGS)"
	@echo "CPPFLAGS origin = $(origin CPPFLAGS)"
	@echo "LDFLAGS  origin = $(origin LDFLAGS)"
	@echo "LDLIBS   origin = $(origin LDLIBS)"
	@echo "HTSSRC   origin = $(origin HTSSRC)"

.PHONY: clean
clean:
	rm -f *.o *.d $(PROGRAMS) $(PROGRAMS_MISC) *~
