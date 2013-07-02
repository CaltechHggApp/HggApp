ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs) -lMLP -lXMLIO -lTMVA
#ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs) -L TMVA/lib
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs) -lTMVA -lMLP -lXMLIO

FASTJETFLAGS = $(shell FASTJET/bin/fastjet-config --cxxflags)
FASTJETLIBS  = $(shell FASTJET/bin/fastjet-config --libs --plugins)

CXX           = g++ -m64
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2 -Xlinker -zmuldefs -Wall -Wno-error=unused-variable -Wno-error=sign-compare -Wno-error=unused-value -Wno-error=unused-but-set-variable
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
NGLIBS        += $(FASTJETLIBS)
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./include
INCLUDEDIRCOMMON = ./
MITINCLUDE       = MitPhysics/Utils/interface/
#INCLUDEDIRTMVA   = ./TMVA/include
SRCDIR           = ./src/
MITSRC           = ./MitPhysics/Utils/src/
#CXX	         += -I$(INCLUDEDIR) -I$(INCLUDEDIRCOMMON) -I$(INCLUDEDIRTMVA) -I.
CXX	         += -I$(INCLUDEDIR) -I$(INCLUDEDIRCOMMON) -I.
OUTLIB	         = ./lib/
H2GLIB           = ./h2glib/
OUTLIBCOMMON     = $(INCLUDEDIRCOMMON)/CommonTools/lib/
OUTLIBEGAMMA	 = $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/lib/

.SUFFIXES: .cc,.C, .hh
.PREFIXES: ./lib/

all:  lib HggApp HggSelectorApp

lib: 	$(OUTLIBCOMMON)Conditions.o \
	$(OUTLIBCOMMON)Selection.o \
	$(OUTLIBCOMMON)Counters.o \
	$(OUTLIBCOMMON)TriggerMask.o \
	$(OUTLIBCOMMON)Utils.o \
	$(OUTLIBCOMMON)Skimmer.o \
	$(OUTLIBCOMMON)EfficiencyEvaluator.o \
	$(OUTLIBCOMMON)CutBasedEleIDSelector.o \
	$(OUTLIBCOMMON)EcalCleaner.o \
	$(OUTLIBEGAMMA)ElectronTrackerIsolation.o \
	$(OUTLIBEGAMMA)ElectronCaloIsolation.o \
	$(OUTLIBEGAMMA)ElectronBestCandidateSelector.o \
	$(OUTLIBEGAMMA)LikelihoodPdf.o \
	$(OUTLIBEGAMMA)LikelihoodSpecies.o \
	$(OUTLIBEGAMMA)LikelihoodPdfProduct.o \
	$(OUTLIBEGAMMA)ElectronLikelihood.o \
	$(OUTLIB)HggReducer.o \
	$(OUTLIB)HggEnergyScale.o \
	$(OUTLIB)HggVertexing.o \
	$(OUTLIB)HggEGEnergyCorrector.o \
	$(OUTLIB)VecbosEGObject.o \
	$(OUTLIB)HggMCWeight.o \
	$(OUTLIB)HggSelector.o \
	$(OUTLIB)HggScaling.o \
	$(OUTLIB)ArgParser.o 


# analysis functions
HggApp: $(SRCDIR)HggApp.C \
	$(OUTLIBCOMMON)Conditions.o \
	$(OUTLIBCOMMON)Selection.o \
	$(OUTLIBCOMMON)Counters.o \
	$(OUTLIBCOMMON)TriggerMask.o \
	$(OUTLIBCOMMON)Utils.o \
	$(OUTLIBCOMMON)Skimmer.o \
	$(OUTLIBEGAMMA)ElectronTrackerIsolation.o \
	$(OUTLIBEGAMMA)ElectronCaloIsolation.o \
	$(OUTLIBEGAMMA)ElectronBestCandidateSelector.o \
	$(OUTLIBEGAMMA)LikelihoodPdf.o \
	$(OUTLIBEGAMMA)LikelihoodSpecies.o \
	$(OUTLIBEGAMMA)LikelihoodPdfProduct.o \
	$(OUTLIBEGAMMA)ElectronLikelihood.o \
	$(OUTLIBCOMMON)EfficiencyEvaluator.o \
	$(OUTLIB)VecbosEGObject.o \
	$(OUTLIB)HggEnergyScale.o \
	$(OUTLIB)HggReducer.o \
	$(OUTLIB)HggMakePhotonTree.o \
	$(OUTLIB)GBRTree.o \
	$(OUTLIB)GBRForest.o \
	$(OUTLIB)ArgParser.o
	$(CXX) $(CXXFLAGS) -o HggApp $(OUTLIB)/*.o $(OUTLIBCOMMON)/*o $(OUTLIBEGAMMA)/*o $(GLIBS) $ $<	

HggSelectorApp: $(SRCDIR)HggSelectorApp.C \
	$(OUTLIBCOMMON)Conditions.o \
	$(OUTLIBCOMMON)Selection.o \
	$(OUTLIBCOMMON)Counters.o \
	$(OUTLIBCOMMON)TriggerMask.o \
	$(OUTLIBCOMMON)Utils.o \
	$(OUTLIBCOMMON)Skimmer.o \
	$(OUTLIBEGAMMA)ElectronTrackerIsolation.o \
	$(OUTLIBEGAMMA)ElectronCaloIsolation.o \
	$(OUTLIBEGAMMA)ElectronBestCandidateSelector.o \
	$(OUTLIBEGAMMA)LikelihoodPdf.o \
	$(OUTLIBEGAMMA)LikelihoodSpecies.o \
	$(OUTLIBEGAMMA)LikelihoodPdfProduct.o \
	$(OUTLIBEGAMMA)ElectronLikelihood.o \
	$(OUTLIBCOMMON)EfficiencyEvaluator.o \
	$(OUTLIB)VecbosEGObject.o \
	$(OUTLIB)HggMCWeight.o \
	$(OUTLIB)HggSelector.o \
	$(OUTLIB)ArgParser.o
	$(CXX) $(CXXFLAGS) -o HggSelectorApp $(OUTLIB)/*.o $(OUTLIBCOMMON)/*o $(OUTLIBEGAMMA)/*o $(GLIBS) $ $<	

HggEfficiencyMapApp: $(SRCDIR)HggEfficiencyMapApp.C \
	$(OUTLIBCOMMON)Conditions.o \
	$(OUTLIBCOMMON)Selection.o \
	$(OUTLIBCOMMON)Counters.o \
	$(OUTLIBCOMMON)TriggerMask.o \
	$(OUTLIBCOMMON)Utils.o \
	$(OUTLIBCOMMON)Skimmer.o \
	$(OUTLIBEGAMMA)ElectronTrackerIsolation.o \
	$(OUTLIBEGAMMA)ElectronCaloIsolation.o \
	$(OUTLIBEGAMMA)ElectronBestCandidateSelector.o \
	$(OUTLIBEGAMMA)LikelihoodPdf.o \
	$(OUTLIBEGAMMA)LikelihoodSpecies.o \
	$(OUTLIBEGAMMA)LikelihoodPdfProduct.o \
	$(OUTLIBEGAMMA)ElectronLikelihood.o \
	$(OUTLIBCOMMON)EfficiencyEvaluator.o \
	$(OUTLIB)VecbosEGObject.o \
	$(OUTLIB)HggMCWeight.o \
	$(OUTLIB)HggEfficiencyMap.o
	$(CXX) $(CXXFLAGS) -o HggEfficiencyMapApp $(OUTLIB)/*.o $(OUTLIBCOMMON)/*o $(OUTLIBEGAMMA)/*o $(GLIBS) $ $<	

ZeeSelectorApp: $(SRCDIR)ZeeSelectorApp.C \
	$(OUTLIBCOMMON)Conditions.o \
	$(OUTLIBCOMMON)Selection.o \
	$(OUTLIBCOMMON)Counters.o \
	$(OUTLIBCOMMON)TriggerMask.o \
	$(OUTLIBCOMMON)Utils.o \
	$(OUTLIBCOMMON)Skimmer.o \
	$(OUTLIBEGAMMA)ElectronTrackerIsolation.o \
	$(OUTLIBEGAMMA)ElectronCaloIsolation.o \
	$(OUTLIBEGAMMA)ElectronBestCandidateSelector.o \
	$(OUTLIBEGAMMA)LikelihoodPdf.o \
	$(OUTLIBEGAMMA)LikelihoodSpecies.o \
	$(OUTLIBEGAMMA)LikelihoodPdfProduct.o \
	$(OUTLIBEGAMMA)ElectronLikelihood.o \
	$(OUTLIBCOMMON)EfficiencyEvaluator.o \
	$(OUTLIB)HggEGEnergyCorrector.o \
	$(OUTLIB)VecbosEGObject.o \
	$(OUTLIB)ZeeSelector.o
	$(CXX) $(CXXFLAGS) -o ZeeSelectorApp $(OUTLIB)/*.o $(OUTLIBCOMMON)/*o $(OUTLIBEGAMMA)/*o $(GLIBS) $ $<	

MuMuGammaSelectorApp: $(SRCDIR)MuMuGammaSelectorApp.C \
	$(OUTLIBCOMMON)Conditions.o \
	$(OUTLIBCOMMON)Selection.o \
	$(OUTLIBCOMMON)Counters.o \
	$(OUTLIBCOMMON)TriggerMask.o \
	$(OUTLIBCOMMON)Utils.o \
	$(OUTLIBCOMMON)Skimmer.o \
	$(OUTLIBEGAMMA)ElectronTrackerIsolation.o \
	$(OUTLIBEGAMMA)ElectronCaloIsolation.o \
	$(OUTLIBEGAMMA)ElectronBestCandidateSelector.o \
	$(OUTLIBEGAMMA)LikelihoodPdf.o \
	$(OUTLIBEGAMMA)LikelihoodSpecies.o \
	$(OUTLIBEGAMMA)LikelihoodPdfProduct.o \
	$(OUTLIBEGAMMA)ElectronLikelihood.o \
	$(OUTLIBCOMMON)EfficiencyEvaluator.o \
	$(OUTLIB)VecbosEGObject.o \
	$(OUTLIB)MuMuGammaSelector.o
	$(CXX) $(CXXFLAGS) -o MuMuGammaSelectorApp $(OUTLIB)/*.o $(OUTLIBCOMMON)/*o $(OUTLIBEGAMMA)/*o $(GLIBS) $ $<	

$(OUTLIB)HggSelector.o: $(SRCDIR)HggSelector.cc \
			$(OUTLIB)HggMassResolution.o \
			$(OUTLIB)HggPhotonID.o \
			$(OUTLIB)Vecbos.o \
			$(OUTLIB)VecbosEGObject.o \
			$(OUTLIB)HggDict.o \
			$(OUTLIB)GBRTree.o \
			$(OUTLIB)HggMCWeight.o \
			$(OUTLIB)HggPhotonID.o \
			$(OUTLIB)HggMuonID.o \
			$(OUTLIB)HggMassResolution.o \
			$(OUTLIB)HggEnergyScale.o \
			$(OUTLIB)MitDict.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HggSelector.o $<

$(OUTLIB)HggEfficiencyMap.o: $(SRCDIR)HggEfficiencyMap.cc \
			$(OUTLIB)HggPhotonID.o \
			$(OUTLIB)HggMuonID.o \
			$(OUTLIB)Vecbos.o \
			$(OUTLIB)VecbosEGObject.o \
			$(OUTLIB)HggDict.o \
			$(OUTLIB)MitDict.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HggEfficiencyMap.o $<

$(OUTLIB)ZeeSelector.o: $(SRCDIR)ZeeSelector.cc \
			$(OUTLIB)Vecbos.o \
			$(OUTLIB)VecbosEGObject.o \
			$(OUTLIB)HggDict.o \
			$(OUTLIB)GBRTree.o \
			$(OUTLIB)MitDict.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)ZeeSelector.o $<

$(OUTLIB)MuMuGammaSelector.o: $(SRCDIR)MuMuGammaSelector.cc \
			$(OUTLIB)Vecbos.o \
			$(OUTLIB)VecbosEGObject.o \
			$(OUTLIB)HggDict.o \
			$(OUTLIB)GBRTree.o \
			$(OUTLIB)HggPhotonID.o \
			$(OUTLIB)MitDict.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)MuMuGammaSelector.o $<

$(OUTLIB)HggPhotonID.o: $(SRCDIR)HggPhotonID.cc \
			$(OUTLIB)VecbosEGObject.o \
			$(OUTLIB)HggDict.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HggPhotonID.o $<

$(OUTLIB)HggMuonID.o: $(SRCDIR)HggMuonID.cc \
			$(OUTLIB)VecbosEGObject.o \
			$(OUTLIB)HggDict.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HggMuonID.o $<

$(OUTLIB)HggMCWeight.o: $(SRCDIR)HggMCWeight.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HggMCWeight.o $<


$(OUTLIB)HggMassResolution.o: $(SRCDIR)HggMassResolution.cc \
				$(OUTLIB)VecbosEGObject.o \
				$(OUTLIB)HggDict.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HggMassResolution.o $<

$(OUTLIB)HggEGEnergyCorrector.o: $(SRCDIR)HggEGEnergyCorrector.cc \
				$(OUTLIB)Vecbos.o \
				$(OUTLIB)GBRTree.o \
				$(OUTLIB)GBRForest.o\
				$(OUTLIB)MitDict.o \
				$(OUTLIB)VecbosEGObject.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HggEGEnergyCorrector.o $<

$(OUTLIB)HggScaling.o: $(SRCDIR)HggScaling.cc \
			$(OUTLIB)VecbosEGObject.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HggScaling.o $<

$(OUTLIB)HggVertexing.o: $(SRCDIR)HggVertexingNew.cc \
			$(OUTLIB)Vecbos.o \
			$(OUTLIB)VecbosEGObject.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HggVertexing.o $<

$(OUTLIB)HggMakePhotonTree.o: $(SRCDIR)HggMakePhotonTree.cc \
			      $(OUTLIB)Vecbos.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HggMakePhotonTree.o $<

$(OUTLIB)HggEnergyScale.o: $(SRCDIR)HggEnergyScale.cc \
			$(OUTLIB)Vecbos.o \
			$(OUTLIB)MitDict.o \
			$(OUTLIB)VecbosEGObject.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HggEnergyScale.o $<

$(OUTLIB)HggReducer.o: $(SRCDIR)HggReducer.cc \
		$(OUTLIB)HggEGEnergyCorrector.o \
		$(OUTLIB)HggVertexing.o \
		$(OUTLIB)HggEnergyScale.o \
		$(OUTLIB)HggDict.o \
		$(OUTLIB)GBRTree.o \
		$(OUTLIB)GBRForest.o \
		$(OUTLIB)MitDict.o \
		$(OUTLIB)VecbosEGObject.o \
		$(OUTLIB)HggScaling.o \
		$(OUTLIB)VecbosJetCorrector.o
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HggReducer.o $(OUTLIB)HggDict.o $<

$(OUTLIB)VecbosEGObject.o: $(SRCDIR)VecbosEGObject.cc
		$(CXX) $(CXXFLAGS) -c -fPIC -I$(INCLUDEDIR) -o $(OUTLIB)VecbosEGObject.o $<

HggDict.cc: LinkDef.h $(OUTLIB)VecbosEGObject.o
	rootcint -l -f HggDict.cc -c -I$(INCLUDEDIR) -p $(INCLUDEDIR)/VecbosEGObject.hh LinkDef.h

$(OUTLIB)HggDict.o: HggDict.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HggDict.o $<

$(OUTLIB)GBRForest.o: $(SRCDIR)GBRForest.cxx $(OUTLIB)GBRTree.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)GBRForest.o $<

$(OUTLIB)GBRTree.o: $(SRCDIR)GBRTree.cxx
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)GBRTree.o $<

MitPhysicsDict.cc: MitPhysicsUtilsLinkDef.h $(OUTLIB)GBRForest.o $(OUTLIB)GBRTree.o
	rootcint -l -f MitPhysicsDict.cc -c -I$(INCLUDEDIR) -p $(INCLUDEDIR)/GBRForest.h $(INCLUDEDIR)/GBRTree.h MitPhysicsUtilsLinkDef.h

$(OUTLIB)MitDict.o: MitPhysicsDict.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)MitDict.o $<


$(OUTLIB)JetCorrectorParameters.o: $(SRCDIR)JetCorrectorParameters.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)JetCorrectorParameters.o $<

$(OUTLIB)SimpleJetCorrectionUncertainty.o: $(SRCDIR)SimpleJetCorrectionUncertainty.cc \
	$(OUTLIB)JetCorrectorParameters.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)SimpleJetCorrectionUncertainty.o $<

$(OUTLIB)JetCorrectionUncertainty.o: $(SRCDIR)JetCorrectionUncertainty.cc \
	$(OUTLIB)SimpleJetCorrectionUncertainty.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)JetCorrectionUncertainty.o $<

$(OUTLIB)SimpleJetCorrector.o: $(SRCDIR)SimpleJetCorrector.cc \
	$(OUTLIB)JetCorrectionUncertainty.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)SimpleJetCorrector.o $<

$(OUTLIB)FactorizedJetCorrector.o: $(SRCDIR)FactorizedJetCorrector.cc \
	$(OUTLIB)JetCorrectionUncertainty.o \
	$(OUTLIB)SimpleJetCorrector.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)FactorizedJetCorrector.o $<

$(OUTLIB)VecbosJetCorrector.o: $(SRCDIR)VecbosJetCorrector.cc \
	$(OUTLIB)JetCorrectionUncertainty.o \
	$(OUTLIB)FactorizedJetCorrector.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)VecbosJetCorrector.o $<

$(OUTLIB)Jet.o: $(SRCDIR)Jet.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)Jet.o $<

$(OUTLIB)MET.o: $(SRCDIR)MET.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)MET.o $<

$(OUTLIB)CaloTower.o: $(SRCDIR)CaloTower.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)CaloTower.o $<

# auxiliary functions to compute selections/efficiencies
$(OUTLIBCOMMON)Conditions.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Conditions.C
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIBCOMMON)Conditions.o $<
$(OUTLIBCOMMON)Utils.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Utils.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIBCOMMON)Utils.o $<
$(OUTLIBCOMMON)Skimmer.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Skimmer.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIBCOMMON)Skimmer.o $<
$(OUTLIBCOMMON)Counters.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Counters.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIBCOMMON)Counters.o $<
$(OUTLIBCOMMON)Selection.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Selection.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIBCOMMON)Selection.o $<
$(OUTLIBCOMMON)TriggerMask.o: $(INCLUDEDIRCOMMON)/CommonTools/src/TriggerMask.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIBCOMMON)TriggerMask.o $<
$(OUTLIBCOMMON)EfficiencyEvaluator.o: $(INCLUDEDIRCOMMON)/CommonTools/src/EfficiencyEvaluator.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)EfficiencyEvaluator.o $<
$(OUTLIBCOMMON)CutBasedEleIDSelector.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/CutBasedEleIDSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)CutBasedEleIDSelector.o $<
$(OUTLIBCOMMON)EcalCleaner.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/EcalCleaner.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)EcalCleaner.o $<
$(OUTLIBEGAMMA)ElectronTrackerIsolation.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/ElectronTrackerIsolation.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBEGAMMA)ElectronTrackerIsolation.o $<
$(OUTLIBEGAMMA)ElectronCaloIsolation.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/ElectronCaloIsolation.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBEGAMMA)ElectronCaloIsolation.o $<
$(OUTLIBEGAMMA)ElectronBestCandidateSelector.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/ElectronBestCandidateSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBEGAMMA)ElectronBestCandidateSelector.o $<
$(OUTLIBEGAMMA)LikelihoodPdf.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/LikelihoodPdf.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBEGAMMA)LikelihoodPdf.o $<
$(OUTLIBEGAMMA)LikelihoodSpecies.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/LikelihoodSpecies.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBEGAMMA)LikelihoodSpecies.o $<
$(OUTLIBEGAMMA)LikelihoodPdfProduct.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/LikelihoodPdfProduct.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBEGAMMA)LikelihoodPdfProduct.o $<
$(OUTLIBEGAMMA)ElectronLikelihood.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/ElectronLikelihood.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBEGAMMA)ElectronLikelihood.o $<

$(OUTLIB)VecbosBase.o: $(SRCDIR)VecbosBase.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)VecbosBase.o $<

$(OUTLIB)ArgParser.o: $(SRCDIR)ArgParser.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)ArgParser.o $<


$(OUTLIB)Vecbos.o: $(SRCDIR)Vecbos.cc \
		$(OUTLIB)VecbosBase.o \
		$(OUTLIB)CaloTower.o \
		$(OUTLIB)Jet.o \
		$(OUTLIB)JetCorrectionUncertainty.o \
		$(OUTLIB)MET.o 
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)Vecbos.o $<

# GEN VECBOS STUSY#
#$(OUTLIB)GenVecbos.o: $(SRCDIR)GenVecbos.C
#	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)GenVecbos.o $<

#$(OUTLIB)GenWjets.o: $(SRCDIR)GenWjets.C $(OUTLIB)GenVecbos.o $(OUTLIB)Vecbos.o
#	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)GenWjets.o $<

#$(OUTLIB)GenZjets.o: $(SRCDIR)GenZjets.C $(OUTLIB)GenVecbos.o $(OUTLIB)Vecbos.o
#	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)GenZjets.o $<

VecbosApp.clean:
	rm -f VecbosApp

clean:
	rm -f $(OUTLIB)*.o $(OUTLIBCOMMON)*.o $(OUTLIBEGAMMA)*.o
	rm -f HggDict.cc
	rm -f HggDict.h
	rm -f MitPhysicsDict.cc
	rm -f MitPhysicsDict.h
	rm -f VecbosApp
	rm -f HggApp
	rm -f HggSelectorApp
