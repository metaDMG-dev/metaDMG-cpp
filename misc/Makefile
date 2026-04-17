# Makefile for the misc directory
# Builds utilities for metaDMG-cpp, using htslib and shared.o from parent

CC  ?= gcc
CXX ?= g++

FLAGS     := -O2
CFLAGS    := $(FLAGS)
CXXFLAGS  := $(FLAGS)
CPPFLAGS  := $(CPPFLAGS) -Wall -Wextra
LDFLAGS   := $(LDFLAGS)
LIBS      := -lz -lm -lpthread 
LDHTS     := -lbz2 -llzma -lcurl

# --- Brugervalgte paths og biblioteker ---

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
  LIBS += $(EXTRA_LIBS)
endif


# --- Crypto library detection ---
HAVE_CRYPTO := $(shell echo 'int main(){}'|$(CXX) -x c++ - -lcrypto -o /dev/null 2>/dev/null && echo 0 || echo 1)
ifeq ($(HAVE_CRYPTO),0)
  $(info Crypto library is available to link; adding -lcrypto)
  LIBS += -lcrypto
else
  $(info Crypto library is not available to link; skipping -lcrypto)
endif

$(info [MISC DEBUG] HTSSRC='$(HTSSRC)')
$(info [MISC DEBUG] PREFIX='$(PREFIX)')
$(info [MISC DEBUG] EXTRA_INC='$(EXTRA_INC)')
$(info [MISC DEBUG] EXTRA_LIB='$(EXTRA_LIB)')

ifeq ($(HTSSRC),systemwide)
  $(info Using systemwide htslib installation)
  LIBS += -lhts
  LIBHTS :=
else
  # Convert relative HTSSRC to absolute if needed
  ifneq ($(filter /%,$(HTSSRC)),$(HTSSRC))
    override HTSSRC := $(shell [ -d "$(HTSSRC)" ] && cd "$(HTSSRC)" && pwd || echo "__INVALID__")
    ifeq ($(HTSSRC),__INVALID__)
      $(error HTSSRC is not a valid directory in misc folder)
    endif
    $(info Converted HTSSRC to absolute path: $(HTSSRC))
  endif
  CPPFLAGS += -I$(HTSSRC)
  LIBHTS := $(HTSSRC)/libhts.a
  LIBS := $(LIBHTS) $(LIBS) $(LDHTS)
endif

PROGRAMS := compressbam extract_reads
PROGRAMS_MISC := MAP dfit_optim
OBJ := kmath.o Likelihood.o M3Load.o

all: $(PROGRAMS) $(PROGRAMS_MISC)

# Ensure htslib and shared.o are available
$(PROGRAMS) $(PROGRAMS_MISC): $(LIBHTS)
$(OBJ): $(LIBHTS)

# Check if libhts.a exists
$(LIBHTS):
	@if [ ! -f "$(LIBHTS)" ]; then \
		echo "Error: $(LIBHTS) not found; ensure parent Makefile has built htslib"; \
		exit 1; \
	fi

# Check if shared.o exists for extract_reads
../shared.o:
	@if [ ! -f "../shared.o" ]; then \
		echo "Error: ../shared.o not found; ensure parent Makefile has built it"; \
		exit 1; \
	fi

kmath.o: bfgs.cpp bfgs.h
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) bfgs.cpp -o kmath.o
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) bfgs.cpp > kmath.d

Likelihood.o: Likelihood.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) Likelihood.cpp -o Likelihood.o
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) Likelihood.cpp > Likelihood.d

M3Load.o: M3Load.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) M3Load.cpp -o M3Load.o
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) M3Load.cpp > M3Load.d

compressbam: compressbam.cpp
	@echo "[MISC DEBUG] LINK compressbam"
	@echo "[MISC DEBUG] CXX=$(CXX)"
	@echo "[MISC DEBUG] CPPFLAGS=$(CPPFLAGS)"
	@echo "[MISC DEBUG] CXXFLAGS=$(CXXFLAGS)"
	@echo "[MISC DEBUG] LIBS=$(LIBS)"
	@echo "[MISC DEBUG] LDFLAGS=$(LDFLAGS)"
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) compressbam.cpp -o compressbam $(LIBS) $(LDFLAGS)

extract_reads: extract_reads.cpp ../shared.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) extract_reads.cpp ../shared.o -o extract_reads $(LIBS) $(LDFLAGS)

MAP: MAP.cpp M3Load.o kmath.o Likelihood.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -std=c++11 MAP.cpp M3Load.o kmath.o Likelihood.o -o MAP $(LIBS) $(LDFLAGS) -D__WITH_MAIN__

dfit_optim: dfit_optim.cpp kmath.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) dfit_optim.cpp kmath.o -o dfit_optim $(LIBS) $(LDFLAGS) -D__WITH_MAIN__

# --- Automatisk afhængighedshåndtering ---
-include $(OBJ:.o=.d)

.PHONY: clean
clean:
	rm -f *.o *.d $(PROGRAMS) $(PROGRAMS_MISC) *~
