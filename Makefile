#modied from htslib makefile
FLAGS2=-O3 -std=c++11 

CFLAGS += $(FLAGS2)
CXXFLAGS += $(FLAGS2)

CSRC = $(wildcard *.c) 
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)

all: metaDMG-cpp

PACKAGE_VERSION  = 0.2

ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git describe --always --dirty)
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif

version.h:
	echo '#define METADAMAGE_VERSION "$(PACKAGE_VERSION)"' > $@


# Adjust $(HTSSRC) to point to your top-level htslib directory
ifdef HTSSRC
$(info HTSSRC defined)
HTS_INCDIR=$(realpath $(HTSSRC))
HTS_LIBDIR=$(realpath $(HTSSRC))/libhts.a
else
$(info HTSSRC not defined, assuming systemwide installation -lhts)
endif


-include $(OBJ:.o=.d)

ifdef LDFLAGS
FLAGS += $(LDFLAGS)
else
FLAGS += $(FLAGS2)
endif

ifdef HTSSRC
%.o: %.c
	$(CC) -c  $(CFLAGS) -I$(HTS_INCDIR) $*.c
	$(CC) -MM $(CFLAGS)  -I$(HTS_INCDIR) $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  -I$(HTS_INCDIR) $*.cpp
	$(CXX) -MM $(CXXFLAGS)  -I$(HTS_INCDIR) $*.cpp >$*.d

metaDMG-cpp: version.h $(OBJ)
	$(CXX) $(FLAGS) -o metaDMG-cpp *.o $(HTS_LIBDIR) -lz -llzma -lbz2 -lpthread -lcurl -lgsl
else
%.o: %.c
	$(CC) -c  $(CFLAGS)  $*.c
	$(CC) -MM $(CFLAGS)  $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  $*.cpp
	$(CXX) -MM $(CXXFLAGS)  $*.cpp >$*.d

metaDMG-cpp: version.h $(OBJ)
	$(CXX) $(FLAGS)  -o metaDMG-cpp *.o -lz -llzma -lbz2 -lpthread -lcurl -lhts -lgsl
endif

clean:	
	rm  -f metaDMG-cpp *.o *.d version.h

force:

test:
	
