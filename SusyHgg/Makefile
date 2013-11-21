ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs) -lMLP -lXMLIO -lTMVA
#ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs) -L TMVA/lib
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs) -lTMVA -lMLP -lXMLIO

ROOFITLIBS    = -L$(ROOFITSYS)/lib -lRooFitCore -lRooFit -lRooStats -lMinuit -lFoam

CXX           = g++ -m64
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2 -Xlinker -zmuldefs -Wall -Wno-error=unused-variable -Wno-error=sign-compare -Wno-error=unused-value -Wno-error=unused-but-set-variable -fopenmp -std=gnu++0x
LD            = g++ -m64
LDFLAGS       = -g
SOFLAGS       = -shared

#PG da qui per macosx
#PG -----------------

ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit
gGLIBS          = $(filter-out -lNew, $(NGLIBS))

CXXFLAGS      += $(ROOTCFLAGS)
CXXFLAGS      += $(FASTJETFLAGS)
LIBS           = $(ROOTLIBS)

NGLIBS         = $(ROOTGLIBS) 
#NGLIBS        += -lMinuit -lTMVA.1 -lMLP -lTreePlayer
NGLIBS        += -lMinuit -lMLP -lTreePlayer
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./include/
INCLUDEDIRCOMMON = ./
INCLUDEDIRROOFIT = $(ROOFITSYS)/include/

HGGAPPDIR        = ../
HGGAPPLIB        = $(HGGAPPDIR)lib/
HGGAPPINC        = $(HGGAPPDIR)include/
SRCDIR           = ./src/
EXESRCDIR        = ./exesrc/
SCRIPTS          = ./scripts/
CXX	         += -I$(INCLUDEDIR) -I$(INCLUDEDIRCOMMON) -I. -I$(INCLUDEDIRROOFIT)
OUTLIB	         = ./lib/

.SUFFIXES: .cc,.C, .hh
.PREFIXES: ./lib/

.PHONY: all
all: CreateWorkspace CreateFits CountEntries

.PHONY: clean
clean: 
	rm -f $(OUTLIB)*.o
	rm -f CreateWorkspace CreateFits CountEntries


CreateWorkspace:	$(SRCDIR)CreateWorkspace.C \
			$(HGGAPPLIB)ArgParser.o
	$(CXX) $(CXXFLAGS) $(ROOFITLIBS) -I $(HGGAPPINC) -I $(INCLUDEDIR) -o $@ $(GLIBS) $^ $<	

CreateFits:		$(SRCDIR)CreateFits.C \
			$(OUTLIB)Fitter.o \
			$(HGGAPPLIB)ArgParser.o
	$(CXX) $(CXXFLAGS) $(ROOFITLIBS) -I $(HGGAPPINC) -I $(INCLUDEDIR) -o $@ $(GLIBS) $(HGGAPPLIB)ArgParser.o $(OUTLIB)*.o $<	

CountEntries:		$(SRCDIR)CountEntries.C \
			$(HGGAPPLIB)ArgParser.o
	$(CXX) $(CXXFLAGS) $(ROOFITLIBS) -I $(HGGAPPINC) -I $(INCLUDEDIR) -o $@ $(GLIBS) $(HGGAPPLIB)ArgParser.o $<	



$(OUTLIB)Fitter.o:		$(SRCDIR)Fitter.cpp
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<
