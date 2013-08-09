##
## 
##

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTGLIBS      = $(shell $(ROOTSYS)/bin/root-config --glibs) -lMLP -lXMLIO -lTMVA -lTreePlayer

CXX           = g++ -m64
CXXFLAGS      = -g -std=c++11 -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2 -Xlinker -zmuldefs -Wall -Wno-error=unused-variable -Wno-error=sign-compare -Wno-error=unused-value -Wno-error=unused-but-set-variable
CXXFLAGS     += $(ROOTCFLAGS) 
LD            = g++ -m64

GLIBS         = $(ROOTGLIBS)

INCLUDEDIR    = ./
SRCDIR        = ./src/
BINDIR        = ./bin/
LIBDIR        = ./lib/

HGGAPPLIB     = ../lib/

all: $(BINDIR)MakeAllVarPlots
.PHONY: all

$(LIBDIR)weightManager.o: $(SRCDIR)weightManager.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<

$(LIBDIR)varCorrector.o: $(SRCDIR)varCorrector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<

$(LIBDIR)plotManager.o: $(SRCDIR)plotManager.cc \
			$(LIBDIR)weightManager.o \
			$(LIBDIR)varCorrector.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<

$(LIBDIR)HistogramStack.o: $(SRCDIR)HistogramStack.cc 
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<

$(BINDIR)MakeAllVarPlots:  $(SRCDIR)MakeAllVarPlots.C \
			$(LIBDIR)plotManager.o \
			$(SRCDIR)getCategories.C \
			$(LIBDIR)HistogramStack.o \
			$(HGGAPPLIB)ArgParser.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $(LIBDIR)/*.o $(HGGAPPLIB)ArgParser.o $(GLIBS) $<	



