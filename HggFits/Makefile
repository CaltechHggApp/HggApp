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
all: runWorkspaces runFits runPlots runToyWorkspaces runSimpleHypothesisTest runBiasStudy

.PHONY: clean
clean: 
	rm -f $(OUTLIB)*.o
	rm -f runWorkspaces
	rm -f runFits
	rm -f runPlots
	rm -f runToyWorkspaces
	rm -f MyDict.cxx

threadTest: $(EXESRCDIR)threadTest.C
	$(CXX) $(CXXFLAGS) -o threadTest $(GLIBS) $ $<	


runWorkspaces:		$(EXESRCDIR)runWorkspaces.C \
			$(OUTLIB)MakeSpinWorkspace.o \
			$(OUTLIB)MakeSusyHggWorkspace.o \
			$(OUTLIB)MakeSpinFits.o \
			$(OUTLIB)MakeSpinToy.o \
			$(OUTLIB)ArgParser.o \
			$(HGGAPPINC)ReadConfig.hh
	$(CXX) $(CXXFLAGS) $(ROOFITLIBS) -I $(HGGAPPINC) -I $(INCLUDEDIR) -o runWorkspaces $(OUTLIB)/*.o  $(GLIBS) $ $<	

runFits:		$(EXESRCDIR)runFits.C \
			$(OUTLIB)MakeSpinWorkspace.o \
			$(OUTLIB)MakeSpinFits.o \
			$(OUTLIB)ArgParser.o \
			$(HGGAPPINC)ReadConfig.hh
	$(CXX) $(CXXFLAGS) $(ROOFITLIBS) -I $(HGGAPPINC) -I $(INCLUDEDIR) -o runFits $(OUTLIB)/*.o  $(GLIBS) $ $<	

runPlots:		$(EXESRCDIR)runPlots.C \
			$(OUTLIB)MakeSpinPlots.o \
			$(OUTLIB)ArgParser.o \
			$(HGGAPPINC)ReadConfig.hh
	$(CXX) $(CXXFLAGS) $(ROOFITLIBS) -I $(HGGAPPINC) -I $(INCLUDEDIR) -o runPlots $(OUTLIB)/*.o  $(GLIBS) $ $<	

runToyWorkspaces:		$(EXESRCDIR)runToyWorkspaces.C \
				$(OUTLIB)MakeSpinToy.o\
				$(OUTLIB)MakeSpinFits.o\
				$(OUTLIB)ArgParser.o \
				$(HGGAPPINC)ReadConfig.hh
	$(CXX) $(CXXFLAGS) $(ROOFITLIBS) -I $(HGGAPPINC) -I $(INCLUDEDIR) -o runToyWorkspaces $(OUTLIB)/*.o  $(GLIBS) $ $<	

runSimpleHypothesisTest:		$(EXESRCDIR)runSimpleHypothesisTest.C \
				$(OUTLIB)MakeSimpleHypothesisTest.o\
				$(OUTLIB)MakeSpinFits.o\
				$(OUTLIB)ArgParser.o \
				$(HGGAPPINC)ReadConfig.hh
	$(CXX) $(CXXFLAGS) $(ROOFITLIBS) -I $(HGGAPPINC) -I $(INCLUDEDIR) -o $@ $(OUTLIB)/*.o  $(GLIBS) $ $<	

runBiasStudy:		$(EXESRCDIR)runBiasStudy.C \
				$(OUTLIB)MakeBiasStudy.o\
				$(OUTLIB)MakeSpinFits.o\
				$(OUTLIB)ArgParser.o \
				$(HGGAPPINC)ReadConfig.hh
	$(CXX) $(CXXFLAGS) $(ROOFITLIBS) -I $(HGGAPPINC) -I $(INCLUDEDIR) -o $@ $(OUTLIB)/*.o  $(GLIBS) $ $<	

optimizeSigmaE: 		$(SCRIPTS)optimizeSigmaE.C \
				$(OUTLIB)MakeSpinFits.o
	$(CXX) $(CXXFLAGS) $(ROOFITLIBS) -I $(HGGAPPINC) -I $(INCLUDEDIR) -o optimizeSigmaE $(OUTLIB)/*.o  $(GLIBS) $ $<	

$(OUTLIB)ArgParser.o:		$(SRCDIR)ArgParser.cc
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)ArgParser.o $<

$(OUTLIB)HggOutputReader2.o:	$(SRCDIR)HggOutputReader2.C
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HggOutputReader2.o $<

$(OUTLIB)GlobeReader.o:	$(SRCDIR)GlobeReader.C
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)GlobeReader.o $<

$(OUTLIB)MakeSpinWorkspace.o:	$(SRCDIR)MakeSpinWorkspace.C \
				$(OUTLIB)MixSpinDatasets.o \
				$(OUTLIB)HggOutputReader2.o \
				$(SRCDIR)selectionMaps.C \
				$(OUTLIB)GlobeReader.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)MakeSpinWorkspace.o $<

$(OUTLIB)MakeSusyHggWorkspace.o:	$(SRCDIR)MakeSusyHggWorkspace.C \
				$(OUTLIB)MakeSpinWorkspace.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<

$(OUTLIB)MakeSpinFits.o:	$(SRCDIR)MakeSpinFits.C \
				$(OUTLIB)HggOutputReader2.o \
				$(OUTLIB)MakeSpinSPlot.o \
				$(OUTLIB)RooGaussianCorr.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)MakeSpinFits.o $<

$(OUTLIB)MakeSpinPlots.o:	$(SRCDIR)MakeSpinPlots.C \
				$(OUTLIB)MakeSpinFits.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)MakeSpinPlots.o $<

$(OUTLIB)MakeSpinSeparationTest.o:	$(SRCDIR)MakeSpinSeparationTest.C \
					$(OUTLIB)MakeSpinWorkspace.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)MakeSpinSeparationTest.o $<

$(OUTLIB)MakeSpinToyWorkspace.o:	$(SRCDIR)MakeSpinToyWorkspace.C \
					$(OUTLIB)MakeSpinWorkspace.o \
					$(OUTLIB)MakeSpinSeparationTest.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)MakeSpinToyWorkspace.o $<


$(OUTLIB)MakeSpinToy.o:		$(SRCDIR)MakeSpinToy.C \
				$(OUTLIB)MakeSpinFits.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)MakeSpinToy.o $<

$(OUTLIB)MakeSpinSPlot.o:	$(SRCDIR)MakeSpinSPlot.C
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)MakeSpinSPlot.o $<

$(OUTLIB)MixSpinDatasets.o:	$(SRCDIR)MixSpinDatasets.C \
				$(OUTLIB)MakeSpinFits.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)MixSpinDatasets.o $<

$(OUTLIB)MakeSimpleHypothesisTest.o: $(SRCDIR)MakeSimpleHypothesisTest.cc \
									$(OUTLIB)MakeSpinFits.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<

$(OUTLIB)MakeBiasStudy.o: $(SRCDIR)MakeBiasStudy.cc \
									$(OUTLIB)MakeSpinFits.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<

$(OUTLIB)RooGaussianCorr.o:	$(SRCDIR)RooGaussianCorr.cc \
				$(OUTLIB)MyDict.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $@ $<

## rootcint -f test.cxx -c -I./include -I/cvmfs/cms.cern.ch/slc5_amd64_gcc472/lcg/roofit/5.34.07-cms/include -p include/RooGaussianCorr.hh LinkDef.h
MyDict.cxx: LinkDef.h
	rootcint -f $@ -c -I$(INCLUDEDIR) -I$(INCLUDEDIRROOFIT) -p $(INCLUDEDIR)RooGaussianCorr.hh LinkDef.h

$(OUTLIB)MyDict.o: MyDict.cxx
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)MyDict.o $<
