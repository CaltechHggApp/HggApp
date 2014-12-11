ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs) -lMLP -lXMLIO -lTMVA
#ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs) -L TMVA/lib
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs) -lTMVA -lMLP -lXMLIO

ROOFITLIBS    = -L$(ROOFITSYS)/lib -lRooFitCore -lRooFit -lRooStats

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
all: CreateWorkspace CreateFits CreatePlots CountEntries MakeSignalInjection runSusyHgg runSignalInjectionNew runSMSValidation passMetFilters

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

CreatePlots:		$(SRCDIR)CreatePlots.C \
			$(OUTLIB)Plotter.o \
			$(HGGAPPLIB)ArgParser.o
	$(CXX) $(CXXFLAGS) $(ROOFITLIBS) -I $(HGGAPPINC) -I $(INCLUDEDIR) -o $@ $(GLIBS) $(HGGAPPLIB)ArgParser.o $(OUTLIB)*.o $<	

CountEntries:		$(SRCDIR)CountEntries.C \
			$(HGGAPPLIB)ArgParser.o
	$(CXX) $(CXXFLAGS) $(ROOFITLIBS) -I $(HGGAPPINC) -I $(INCLUDEDIR) -o $@ $(GLIBS) $(HGGAPPLIB)ArgParser.o $<	

MakeSignalInjection:	$(SRCDIR)MakeSignalInjection.cpp \
			$(HGGAPPLIB)ArgParser.o
	$(CXX) $(CXXFLAGS) $(ROOFITLIBS) -I $(HGGAPPINC) -I $(INCLUDEDIR) -o $@ $(GLIBS) $(HGGAPPLIB)ArgParser.o $<	

runSusyHgg:		$(SRCDIR)runSusyHgg.cpp\
			$(OUTLIB)DataFitterNew.o \
			$(OUTLIB)MCBackgroundFitterNew.o \
			$(OUTLIB)SMSFitterNew.o \
			$(OUTLIB)FitterNew.o \
			$(OUTLIB)CombinePrep.o \
			$(OUTLIB)SusyHggTreeBase.o
	$(CXX) $(CXXFLAGS) $(ROOFITLIBS) -I $(HGGAPPINC) -I $(INCLUDEDIR) -o $@ $(GLIBS) $(HGGAPPLIB)ArgParser.o $^	

passMetFilters:		$(SRCDIR)passMetFilters.C
	$(CXX) $(CXXFLAGS) $(ROOFITLIBS) -I $(HGGAPPINC) -I $(INCLUDEDIR) -o $@ $(GLIBS) $(HGGAPPLIB)ArgParser.o $^	

runSignalInjectionNew:	$(SRCDIR)runSignalInjectionNew.cpp \
			$(OUTLIB)SignalInjectionNew.o \
			$(OUTLIB)FitterNew.o \
			$(OUTLIB)SusyHggTreeBase.o
	$(CXX) $(CXXFLAGS) $(ROOFITLIBS) -I $(HGGAPPINC) -I $(INCLUDEDIR) -o $@ $(GLIBS) $(HGGAPPLIB)ArgParser.o $^	

runSMSValidation:	$(SRCDIR)runSMSValidation.cpp \
			$(OUTLIB)SMSValidation.o \
			$(OUTLIB)SMSFitterNew.o \
			$(OUTLIB)FitterNew.o \
			$(OUTLIB)SusyHggTreeBase.o
	$(CXX) $(CXXFLAGS) $(ROOFITLIBS) -I $(HGGAPPINC) -I $(INCLUDEDIR) -o $@ $(GLIBS) $(HGGAPPLIB)ArgParser.o $^	

$(OUTLIB)Fitter.o:		$(SRCDIR)Fitter.cpp
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<

$(OUTLIB)Plotter.o:		$(SRCDIR)Plotter.cpp
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<

$(OUTLIB)SusyHggTreeBase.o:		$(SRCDIR)SusyHggTreeBase.cpp
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<

$(OUTLIB)FitterNew.o:		$(SRCDIR)FitterNew.cpp \
				$(OUTLIB)SusyHggTreeBase.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<

$(OUTLIB)DataFitterNew.o:	$(SRCDIR)DataFitterNew.cpp \
				$(OUTLIB)FitterNew.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<

$(OUTLIB)MCBackgroundFitterNew.o:	$(SRCDIR)MCBackgroundFitterNew.cpp \
				$(OUTLIB)FitterNew.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<

$(OUTLIB)SMSFitterNew.o:	$(SRCDIR)SMSFitterNew.cpp \
				$(OUTLIB)FitterNew.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<

$(OUTLIB)CombinePrep.o:		$(SRCDIR)CombinePrep.cpp
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<

$(OUTLIB)SignalInjectionNew.o:	$(SRCDIR)SignalInjectionNew.cpp \
				$(OUTLIB)FitterNew.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<

$(OUTLIB)SMSValidation.o:	$(SRCDIR)SMSValidation.cpp \
				$(OUTLIB)FitterNew.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<

MyDict.cxx: LinkDef.h
	rootcint -f $@ -c -I$(INCLUDEDIR) -I$(INCLUDEDIRROOFIT) -p $(INCLUDEDIR)RooParameterizedPdf.hh LinkDef.h

#$(OUTLIB)libRooParameterizedPdf.so:	MyDict.cxx
#	g++ -shared -o$@ `root-config --ldflags` $(CXXFLAGS) -I$(ROOTSYS)/include $^

$(OUTLIB)MyDict.o: MyDict.cxx
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<

$(OUTLIB)RooParameterizedPdf.so:	$(SRCDIR)RooParameterizedPdf.cc \
				$(OUTLIB)MyDict.o
	$(CXX) $(ROOFITLIBS) -shared $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $(GLIBS) $^
