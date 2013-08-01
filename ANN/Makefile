CXX           = g++ -m64
CXXFLAGS      = -g -std=gnu++11 -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2 -Xlinker -zmuldefs -Wall -Wno-error=unused-variable -Wno-error=sign-compare -Wno-error=unused-value -Wno-error=unused-but-set-variable -pthread -lm -ldl -rdynamic
LD            = g++ -m64
LDFLAGS       = -g
SOFLAGS       = -shared


INCLUDEDIR       = ./include/
SRCDIR           = ./src/
OUTLIB           = ./lib/

.SUFFIXES: .cc,.C, .hh
.PREFIXES: ./lib/

all:  lib main

lib: 	$(OUTLIB)SimpleMatrix.o

# analysis functions
main: 	$(SRCDIR)main.cpp \
	$(OUTLIB)SimpleMatrix.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $(OUTLIB)*.o $ $<	

$(OUTLIB)SimpleMatrix.o: $(SRCDIR)SimpleMatrix.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<


clean:
	rm -f $(OUTLIB)*.o
