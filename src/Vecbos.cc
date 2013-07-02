#include "Vecbos.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"

#include "cajun/json/reader.h"
#include "cajun/json/elements.h"

#include <CommonTools/include/TriggerMask.hh>


using namespace std;
using namespace stdcomb;
using namespace bits;

Vecbos::Vecbos(TTree *tree) : VecbosBase(tree)
{
  // By default, we do NOT want to be verbose.
  verbose = false;
  jsonFile = "";
  lastFile = "";

  //initialize the JES objects for calo and PF
  jecUnc_calo = 
    (JetCorrectionUncertainty*) new JetCorrectionUncertainty("data/JES/START38_V13_AK5Calo_Uncertainty.txt");
  jecUnc_PF = 
    (JetCorrectionUncertainty*) new JetCorrectionUncertainty("data/JES/START38_V13_AK5PF_Uncertainty.txt");

}

Vecbos::~Vecbos()
{
  // By this time, the destructor of VecbosBase has not yet been called.
  // This means that the tree has not yet been deleted.
  // So, we do nothing here.
}

struct Vecbos::JetConfig{
  fastjet::JetDefinition theJetDef;
  fastjet::ActiveAreaSpec theAreaSpec;
};

void Vecbos::setJsonGoodRunList(const string& jsonFilePath)
{
  jsonFile=jsonFilePath;
}

void Vecbos::fillRunLSMap()
{
  
  if (jsonFile == "")
    {
      std::cout << "Cannot fill RunLSMap. json file not configured" << std::endl;
      return;
    }

  std::ifstream jsonFileStream;
  jsonFileStream.open(jsonFile.c_str());
  if (!jsonFileStream.is_open())
    {
      std::cout << "Unable to open file " << jsonFile << std::endl;
      return;
    }

  json::Object elemRootFile;
  json::Reader::Read(elemRootFile, jsonFileStream);

  for (json::Object::const_iterator itRun=elemRootFile.Begin();itRun!=elemRootFile.End();++itRun)
    {
      const json::Array& lsSegment = (*itRun).element;
      LSSegments thisRunSegments; 
      for (json::Array::const_iterator lsIterator=lsSegment.Begin();lsIterator!=lsSegment.End();++lsIterator)
	{
	  json::Array lsSegment=(*lsIterator);
	  json::Number lsStart=lsSegment[0];	   
	  json::Number lsEnd=lsSegment[1];
	  aLSSegment thisSegment;
	  thisSegment.first=lsStart.Value();
	  thisSegment.second=lsEnd.Value();
	  thisRunSegments.push_back(thisSegment);
	  //	   std::pair<int, int> lsSegment=std::pair<int, int>(atoi(,lsIterator[1]); 
	}
      goodRunLS.insert(aRunsLSSegmentsMapElement(atoi((*itRun).name.c_str()),thisRunSegments));
    }


  std::cout << "[GoodRunLSMap]::Good Run LS map filled with " << goodRunLS.size() << " runs" << std::endl;
  for (runsLSSegmentsMap::const_iterator itR=goodRunLS.begin(); itR!=goodRunLS.end(); ++itR)
    {
      std::cout << "[GoodRunLSMap]::Run " << (*itR).first <<  " LS ranges are: ";
      for (LSSegments::const_iterator iSeg=(*itR).second.begin();iSeg!=(*itR).second.end();++iSeg)
	std::cout << "[" << (*iSeg).first << "," << (*iSeg).second << "] "; 
      std::cout << std::endl;
    }
}

bool Vecbos::isGoodRunLS()
{
  //std::cout << "GoodRunLS" << std::endl;
  runsLSSegmentsMap::const_iterator thisRun=goodRunLS.find(runNumber);
  if (thisRun == goodRunLS.end())
    return false;
  //std::cout << runNumber << " found in the good run map" << std::endl;
  for (LSSegments::const_iterator iSeg=goodRunLS[runNumber].begin();iSeg!=goodRunLS[runNumber].end();++iSeg)
    {
      //std::cout << "Range is [" << (*iSeg).first << "," << (*iSeg).second << "]" << std::endl;
      if ( static_cast<unsigned int>(lumiBlock) >= (*iSeg).first && static_cast<unsigned int>(lumiBlock) <= (*iSeg).second)
	return true;
    }
  return false;
}

std::string Vecbos::getHLTPathForRun(int runN, std::string fullname) {
  TString fullName = TString(fullname.c_str());
  TObjArray* selectionTokens = fullName.Tokenize(":");
  if (selectionTokens->GetEntries()!=2) {
    std::cout << "Wrong trigger strings " << selectionTokens->GetEntries() << std::endl;
    return std::string("NOPATH");
  }
  TString RunRange =((TObjString*)(*selectionTokens)[0])->GetString();
  TString HLTPathName =((TObjString*)(*selectionTokens)[1])->GetString();
  
  TObjArray* runs = RunRange.Tokenize("-");
  if (runs->GetEntries()!=2) {
    std::cout << "Wrong trigger run range strings " << runs->GetEntries() << std::endl;
    return std::string("NOPATH");    
  }
  
  const char *minStr = (((TObjString*)(*runs)[0])->GetString()).Data();
  const char *maxStr = (((TObjString*)(*runs)[1])->GetString()).Data();

  int min = atoi(minStr);
  int max = atoi(maxStr);

  if(runN>=min && runN<=max) return std::string(HLTPathName.Data());
  else return std::string("NOPATH");
}


bool Vecbos::reloadTriggerMask(bool newVersion)
{
  //  std::cout << "[ReloadTriggerMask]::Reloading trigger mask" << std::endl;
  if(newVersion) {
    std::vector<int> triggerMask;
    for (std::vector< std::string >::const_iterator fIter=requiredTriggers.begin();fIter!=requiredTriggers.end();++fIter)
      {   
        for(unsigned int i=0; i<nameHLT->size(); i++)
          {
	    //std::cout << i << " " << nameHLT->at(i).c_str() << std::endl;
	    //            if( !strcmp ((*fIter).c_str(), nameHLT->at(i).c_str() ) )
	    if( nameHLT->at(i).find(*fIter) != string::npos)
              {
                triggerMask.push_back( indexHLT[i] ) ;
		//std::cout << "[ReloadTriggerMaskNew]::Requiring bit " << indexHLT[i]  << " " << nameHLT->at(i).c_str() << std::endl;
                break;
              }
          }
      }
    m_requiredTriggers = triggerMask;
    //std::cout << m_requiredTriggers.size() <<" " << std::endl;
  } else {
    TString fileName=((TChain*)fChain)->GetFile()->GetName();
    if ( TString(lastFile) != fileName )
      {

        std::cout << "[ReloadTriggerMask]::File has changed reloading trigger mask" << std::endl;
        lastFile = fileName;
        TTree *treeCond;
        std::cout << "[ReloadTriggerMask]::Opening " << fileName << std::endl;
        treeCond = (TTree*)((TChain*)fChain)->GetFile()->Get("Conditions");
        int           nHLT_;
        std::vector<std::string>  *nameHLT_;
        std::vector<unsigned int> *indexHLT_;

        //To get the pointers for the vectors
        nameHLT_=0;
        indexHLT_=0;

        treeCond->SetBranchAddress("nHLT", &nHLT_);
        treeCond->SetBranchAddress("nameHLT", &nameHLT_);
        treeCond->SetBranchAddress("indexHLT", &indexHLT_);
        treeCond->GetEntry(0);

        std::vector<int> triggerMask;
        for (std::vector< std::string >::const_iterator fIter=requiredTriggers.begin();fIter!=requiredTriggers.end();++fIter)
          {
            for(unsigned int i=0; i<nameHLT_->size(); i++) 
              {
                if( !strcmp ((*fIter).c_str(), nameHLT_->at(i).c_str() ) ) 
                  {
                    triggerMask.push_back( indexHLT_->at(i) ) ;
                    break;
                  }
              }
          }
        m_requiredTriggers = triggerMask;
        for (unsigned int i=0;i<m_requiredTriggers.size();++i)
          std::cout << "[ReloadTriggerMask]::Requiring bit " << m_requiredTriggers[i] << " " << requiredTriggers[i] << std::endl;
      }
  }
  return true;
}

bool Vecbos::reloadTriggerMask(int runN)
{
  std::vector<int> triggerMask;

  // load the triggers required for 
  for (std::vector< std::string >::const_iterator fIter=requiredTriggers.begin();fIter!=requiredTriggers.end();++fIter)
    {   
      std::string pathName = getHLTPathForRun(runN,*fIter);
      for(unsigned int i=0; i<nameHLT->size(); i++)
        {
           if(nameHLT->at(i).find(pathName) != string::npos)
            {
              triggerMask.push_back( indexHLT[i] ) ;
              break;
            }
        }
    }
  m_requiredTriggers = triggerMask;
  return true;
}

void Vecbos::WriteHistos(vector<TH1D*> histos, TFile* file, string dirname){
  file->cd();
  if (!file->GetDirectory(dirname.c_str()))
      file->mkdir(dirname.c_str());
  file->cd(dirname.c_str());
  for(int i=0; i< int(histos.size()); i++) 
    histos[i]->Write();
  file->cd();
}
void Vecbos::WriteHistos(vector<TProfile*> histos, TFile* file, string dirname){
  file->cd();
  if (!file->GetDirectory(dirname.c_str()))
      file->mkdir(dirname.c_str());
  file->cd(dirname.c_str());
  for(int i=0; i< int(histos.size()); i++) 
    histos[i]->Write();
  file->cd();
}
void Vecbos::WriteHistos(vector<TH2D*> histos, TFile* file, string dirname){
  file->cd();
  if (!file->GetDirectory(dirname.c_str()))
      file->mkdir(dirname.c_str());
  file->cd(dirname.c_str());
  for(int i=0; i< int(histos.size()); i++) 
    histos[i]->Write();
  file->cd();
}
void Vecbos::ReadParameters(const char* filename) {
  char buffer[256];  
  ifstream reader (filename);
  map<string, double> ndata; 
  map<string, string> sdata;
  string label, equal, svalue;
  char cvalue[40];
  double value;

  if ( ! reader.is_open())
    { cout << "Error opening file"; exit (1); }
  
  while ( reader.getline (buffer,256) )
    {
      istringstream i(buffer);
      i >> label
	>> equal
	>> cvalue;
      if (equal == "=") {
	value = atof(cvalue);
	ndata[label] = value;
      } else if (equal == "==") {
	svalue = (string) cvalue;
	sdata[label] = svalue;
      }
      equal = "";
    }
  
  map<string, double>::const_iterator iter;

  cout << "---------------------------------------------------" << endl;
  cout << "Application is running with the following setting: " << endl;
  for (iter=ndata.begin(); iter != ndata.end(); iter++) {
    cout << iter->first << " " << iter->second << endl;
  }

  map<string, string>::const_iterator stiter;
  for (stiter=sdata.begin(); stiter != sdata.end(); stiter++) {
    cout << stiter->first << " " << stiter->second << endl;
  }
  cout << "---------------------------------------------------" << endl;
  
  AssignParameters(ndata);
  AssignParameters(sdata);
  
  return;
}

void Vecbos::AssignParameters(map<string, string> sdata) {
  process = string (sdata["process"]);
}

void Vecbos::AssignParameters(map<string, double> ndata) {
  
  // Int Parameters         
  // NExtractions = int (data["NExtractions"]);
  
  // Double Parameters
  barrellimit = ndata["barrellimit"];
  endcaplimit = ndata["endcaplimit"];
}

void Vecbos::InitParameters() {
  process     = "Zjets";
  barrellimit = 1.3;
  endcaplimit = 3.0;
}

vector<Jet> Vecbos::SortJet(vector<Jet> v){
  vector<Jet> sorted;
  vector<pair<double,int> > pT;
  int N = v.size();

  for(int i=0; i<N; i++) {
    double pt = sqrt(pow(v[i].px(),2.)+pow(v[i].py(),2.));
    pT.push_back(std::make_pair(pt,i));
  }

  // sort from smallest to larger, then reverse.
  std::sort(pT.begin(), pT.end());
  std::reverse(pT.begin(), pT.end());
   
  for(unsigned int i=0; i<pT.size(); i++) 
    sorted.push_back(v[pT.at(i).second]);
  
  return sorted;
}

vector<Jet> Vecbos::SortJetByEt(vector<Jet> v){
  vector<Jet> sorted;
  vector<pair<double,int> > pT;
  int N = v.size();

  for(int i=0; i<N; i++) {
    double pt = v[i].et(); // No need to change variable names.
    pT.push_back(std::make_pair(pt,i));
  }

  // sort from smallest to larger, then reverse.
  std::sort(pT.begin(), pT.end());
  std::reverse(pT.begin(), pT.end());
   
  for(unsigned int i=0; i<pT.size(); i++) 
    sorted.push_back(v[pT.at(i).second]);
  
  return sorted;
}

vector<Jet> Vecbos::GetUncorrJets() {

  vector<Jet> jets;
  jets.clear();
  for(int j=0; j<nAK5Jet; j++) {
    float uncorrEt = uncorrEnergyAK5Jet[j]*fabs(sin(thetaAK5Jet[j]));
    TLorentzVector p4;
    p4.SetPtEtaPhiE(uncorrEt,etaAK5Jet[j],phiAK5Jet[j],uncorrEnergyAK5Jet[j]);
    float efrac = emFracAK5Jet[j];
    float hfrac = hadFracAK5Jet[j];
    Jet ak5j(p4,efrac,hfrac);
    jets.push_back(ak5j);
  }

  return jets;

}

vector<Jet> Vecbos::GetCorrJets(int scaleEnergy) {

  vector<Jet> jets;
  jets.clear();
  for(int j=0; j<nAK5Jet; j++) {
    TLorentzVector p4;
    if(scaleEnergy==0) p4.SetPxPyPzE(pxAK5Jet[j],pyAK5Jet[j],pzAK5Jet[j],energyAK5Jet[j]);
    else {
      float mass = sqrt(pow(energyAK5Jet[j],2)-pow(pxAK5Jet[j],2)-pow(pyAK5Jet[j],2)-pow(pzAK5Jet[j],2));
      float ptUnscaled = sqrt(pow(pxAK5Jet[j],2)+pow(pyAK5Jet[j],2));

      // estimate the uncertainty
      jecUnc_PF->setJetEta(etaAK5Jet[j]);
      jecUnc_PF->setJetPt(ptUnscaled);

      // apply the uncertainty
      float pt = ptUnscaled + scaleEnergy*jecUnc_PF->getUncertainty(true)*ptUnscaled;
      float p = pt/fabs(sin(thetaAK5Jet[j]));
      float energy = sqrt(p*p+mass*mass);

      p4.SetPtEtaPhiE(pt,etaAK5Jet[j],phiAK5Jet[j],energy);
    }
    float efrac = emFracAK5Jet[j];
    float hfrac = hadFracAK5Jet[j];
    Jet ak5j(p4,efrac,hfrac);
    jets.push_back(ak5j);
  }

  return jets;

}

vector<Jet> Vecbos::GetGenJets() {

  vector<Jet> jets;
  /*
  jets.clear();
  for(int j=0; j<nAK5GenJet; j++) {
    TLorentzVector p4(pxAK5GenJet[j],pyAK5GenJet[j],pzAK5GenJet[j],energyAK5GenJet[j]);
    float efrac = 0.;
    float hfrac = 0.;
    Jet ak5j(p4,efrac,hfrac);
    ak5j.setPFJetID(0,0,0,0,0,0);
    jets.push_back(ak5j);
  }
  */
  return jets;

}

vector<BTagJet> Vecbos::GetBTagCorrJets() {
  
  vector<BTagJet> btags;
  btags.clear();
  for(int j=0; j<nAK5Jet; j++) {
    BTagJet btag;
    //    btag.combinedSecondaryVertexBJetTags = combinedSecondaryVertexBJetTagsAK5Jet[j];
    //btag.combinedSecondaryVertexMVABJetTags = combinedSecondaryVertexMVABJetTagsAK5Jet[j];
    //btag.jetBProbabilityBJetTags = jetBProbabilityBJetTagsAK5Jet[j];
    //btag.jetProbabilityBJetTags = jetProbabilityBJetTagsAK5Jet[j];
    btag.simpleSecondaryVertexBJetTags = simpleSecondaryVertexHighEffBJetTagsAK5Jet[j];
    //btag.softMuonBJetTags = softMuonBJetTagsAK5Jet[j];
    btag.trackCountingHighPurBJetTags = trackCountingHighPurBJetTagsAK5Jet[j];
    btag.trackCountingHighEffBJetTags = trackCountingHighEffBJetTagsAK5Jet[j];
    btags.push_back(btag);
  }

  return btags;

}

vector<Jet> Vecbos::GetUncorrPFJets() {

  vector<Jet> jets;
  jets.clear();
  /*
  for(int j=0; j<nAK5PFNoPUJet; j++) {
    float uncorrEt = uncorrEnergyAK5PFNoPUJet[j]*fabs(sin(thetaAK5PFNoPUJet[j]));
    TLorentzVector p4;
    p4.SetPtEtaPhiE(uncorrEt,etaAK5PFNoPUJet[j],phiAK5PFNoPUJet[j],uncorrEnergyAK5PFNoPUJet[j]);
    float totenergy = (chargedHadronEnergyAK5PFNoPUJet[j]+neutralHadronEnergyAK5PFNoPUJet[j]+
                       chargedEmEnergyAK5PFNoPUJet[j]+neutralEmEnergyAK5PFNoPUJet[j]);
    float efrac = (chargedEmEnergyAK5PFNoPUJet[j]+neutralEmEnergyAK5PFNoPUJet[j]) / totenergy;
    float hfrac = (chargedHadronEnergyAK5PFNoPUJet[j]+neutralHadronEnergyAK5PFNoPUJet[j]) / totenergy;
    Jet sisconej(p4,efrac,hfrac);

    // PF jet ID variables
    float neutralHadFrac = neutralHadronEnergyAK5PFNoPUJet[j]/energyAK5PFNoPUJet[j];
    float neutralEmFraction = neutralEmEnergyAK5PFNoPUJet[j]/energyAK5PFNoPUJet[j];
    float nConstituents = chargedHadronMultiplicityAK5PFNoPUJet[j] + neutralHadronMultiplicityAK5PFNoPUJet[j] +
      photonMultiplicityAK5PFNoPUJet[j] + electronMultiplicityAK5PFNoPUJet[j] + muonMultiplicityAK5PFNoPUJet[j] +
      HFHadronMultiplicityAK5PFNoPUJet[j] + HFEMMultiplicityAK5PFNoPUJet[j];
    float chargedHadFraction = chargedHadronEnergyAK5PFNoPUJet[j]/energyAK5PFNoPUJet[j];
    float chargedMultiplicity = chargedHadronMultiplicityAK5PFNoPUJet[j] + electronMultiplicityAK5PFNoPUJet[j] + muonMultiplicityAK5PFNoPUJet[j];
    float chargedEmFraction = chargedEmEnergyAK5PFNoPUJet[j]/energyAK5PFNoPUJet[j];
    //sisconej.setPFNoPUJetID(neutralHadFrac,neutralEmFraction,nConstituents,
    //chargedHadFraction,chargedMultiplicity,chargedEmFraction);
    jets.push_back(sisconej);
  }
  */
  return jets;

}

vector<Jet> Vecbos::GetCorrPFJets(int scaleEnergy) {

  vector<Jet> jets;
  jets.clear();
  for(int j=0; j<nAK5PFNoPUJet; j++) {
    TLorentzVector p4;
    if(scaleEnergy==0) p4.SetPxPyPzE(pxAK5PFNoPUJet[j],pyAK5PFNoPUJet[j],pzAK5PFNoPUJet[j],energyAK5PFNoPUJet[j]);
    else {
      float mass = sqrt(pow(energyAK5PFNoPUJet[j],2)-pow(pxAK5PFNoPUJet[j],2)-pow(pyAK5PFNoPUJet[j],2)-pow(pzAK5PFNoPUJet[j],2));
      float ptUnscaled = sqrt(pow(pxAK5PFNoPUJet[j],2)+pow(pyAK5PFNoPUJet[j],2));

      // estimate the uncertainty
      jecUnc_PF->setJetEta(etaAK5PFNoPUJet[j]);
      jecUnc_PF->setJetPt(ptUnscaled);

      // apply the uncertainty
      float pt = ptUnscaled + scaleEnergy*jecUnc_PF->getUncertainty(true)*ptUnscaled;
      float p = pt/fabs(sin(thetaAK5PFNoPUJet[j]));
      float energy = sqrt(p*p+mass*mass);

      p4.SetPtEtaPhiE(pt,etaAK5PFNoPUJet[j],phiAK5PFNoPUJet[j],energy);
    }
    float totenergy = (chargedHadronEnergyAK5PFNoPUJet[j]+neutralHadronEnergyAK5PFNoPUJet[j]+
                       chargedEmEnergyAK5PFNoPUJet[j]+neutralEmEnergyAK5PFNoPUJet[j]);
    float efrac = (chargedEmEnergyAK5PFNoPUJet[j]+neutralEmEnergyAK5PFNoPUJet[j]) / totenergy;
    float hfrac = (chargedHadronEnergyAK5PFNoPUJet[j]+neutralHadronEnergyAK5PFNoPUJet[j]) / totenergy;
    Jet sisconej(p4,efrac,hfrac);

    // PF jet ID variables
    float neutralHadFrac = neutralHadronEnergyAK5PFNoPUJet[j]/energyAK5PFNoPUJet[j];
    float neutralEmFraction = neutralEmEnergyAK5PFNoPUJet[j]/energyAK5PFNoPUJet[j];
    float nConstituents = chargedHadronMultiplicityAK5PFNoPUJet[j] + neutralHadronMultiplicityAK5PFNoPUJet[j] +
      photonMultiplicityAK5PFNoPUJet[j] + electronMultiplicityAK5PFNoPUJet[j] + muonMultiplicityAK5PFNoPUJet[j] +
      HFHadronMultiplicityAK5PFNoPUJet[j] + HFEMMultiplicityAK5PFNoPUJet[j];
    float chargedHadFraction = chargedHadronEnergyAK5PFNoPUJet[j]/energyAK5PFNoPUJet[j];
    float chargedMultiplicity = chargedHadronMultiplicityAK5PFNoPUJet[j] + electronMultiplicityAK5PFNoPUJet[j] + muonMultiplicityAK5PFNoPUJet[j];
    float chargedEmFraction = chargedEmEnergyAK5PFNoPUJet[j]/energyAK5PFNoPUJet[j];
    //sisconej.setPFNoPUJetID(neutralHadFrac,neutralEmFraction,nConstituents,
    //                    chargedHadFraction,chargedMultiplicity,chargedEmFraction);
    jets.push_back(sisconej);
  }

  return jets;

}

vector<BTagJet> Vecbos::GetBTagCorrPFJets()
{
  vector<BTagJet> btags;
  btags.clear();
  for(int j=0; j<nAK5PFNoPUJet; j++) {
    BTagJet btag;
    //btag.combinedSecondaryVertexBJetTags = combinedSecondaryVertexBJetTagsAK5PFNoPUJet[j];
    //btag.combinedSecondaryVertexMVABJetTags = combinedSecondaryVertexMVABJetTagsAK5PFNoPUJet[j];
    //btag.jetBProbabilityBJetTags = jetBProbabilityBJetTagsAK5PFNoPUJet[j];
    //btag.jetProbabilityBJetTags = jetProbabilityBJetTagsAK5PFNoPUJet[j];
    btag.simpleSecondaryVertexBJetTags = simpleSecondaryVertexHighEffBJetTagsAK5PFNoPUJet[j];
    //    btag.softMuonBJetTags = softMuonBJetTagsAK5PFNoPUJet[j];
    btag.trackCountingHighPurBJetTags = trackCountingHighPurBJetTagsAK5PFNoPUJet[j];
    btag.trackCountingHighEffBJetTags = trackCountingHighEffBJetTagsAK5PFNoPUJet[j];
    btags.push_back(btag);
  }

  return btags;

}

vector<Jet> Vecbos::GetUncorrPUPFJets() {

  vector<Jet> jets;
  jets.clear();
  /*
  for(int j=0; j<nAK5PFPUcorrJet; j++) {
    float uncorrEt = uncorrEnergyAK5PFPUcorrJet[j]*fabs(sin(thetaAK5PFPUcorrJet[j]));
    TLorentzVector p4;
    p4.SetPtEtaPhiE(uncorrEt,etaAK5PFPUcorrJet[j],phiAK5PFPUcorrJet[j],uncorrEnergyAK5PFPUcorrJet[j]);
    float totenergy = (chargedHadronEnergyAK5PFPUcorrJet[j]+neutralHadronEnergyAK5PFPUcorrJet[j]+
                       chargedEmEnergyAK5PFPUcorrJet[j]+neutralEmEnergyAK5PFPUcorrJet[j]);
    float efrac = (chargedEmEnergyAK5PFPUcorrJet[j]+neutralEmEnergyAK5PFPUcorrJet[j]) / totenergy;
    float hfrac = (chargedHadronEnergyAK5PFPUcorrJet[j]+neutralHadronEnergyAK5PFPUcorrJet[j]) / totenergy;
    Jet sisconej(p4,efrac,hfrac);

    // PF jet ID variables
    float neutralHadFrac = neutralHadronEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    float neutralEmFraction = neutralEmEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    float nConstituents = chargedHadronMultiplicityAK5PFPUcorrJet[j] + neutralHadronMultiplicityAK5PFPUcorrJet[j] +
      photonMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j] +
      HFHadronMultiplicityAK5PFPUcorrJet[j] + HFEMMultiplicityAK5PFPUcorrJet[j];
    float chargedHadFraction = chargedHadronEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    float chargedMultiplicity = chargedHadronMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j];
    float chargedEmFraction = chargedEmEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    //sisconej.setPFNoPUJetID(neutralHadFrac,neutralEmFraction,nConstituents,
    //                    chargedHadFraction,chargedMultiplicity,chargedEmFraction);
    jets.push_back(sisconej);
  }
  */
  return jets;

}

vector<Jet> Vecbos::GetCorrPUPFJets(int scaleEnergy) {

  vector<Jet> jets;
  jets.clear();
  for(int j=0; j<nAK5PFPUcorrJet; j++) {
    TLorentzVector p4;
    if(scaleEnergy==0) p4.SetPxPyPzE(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j],pzAK5PFPUcorrJet[j],energyAK5PFPUcorrJet[j]);
    else {
      float mass = sqrt(pow(energyAK5PFPUcorrJet[j],2)-pow(pxAK5PFPUcorrJet[j],2)-pow(pyAK5PFPUcorrJet[j],2)-pow(pzAK5PFPUcorrJet[j],2));
      float ptUnscaled = sqrt(pow(pxAK5PFPUcorrJet[j],2)+pow(pyAK5PFPUcorrJet[j],2));

      // estimate the uncertainty
      jecUnc_PF->setJetEta(etaAK5PFPUcorrJet[j]);
      jecUnc_PF->setJetPt(ptUnscaled);

      // apply the uncertainty
      float pt = ptUnscaled + scaleEnergy*jecUnc_PF->getUncertainty(true)*ptUnscaled;
      float p = pt/fabs(sin(thetaAK5PFPUcorrJet[j]));
      float energy = sqrt(p*p+mass*mass);

      p4.SetPtEtaPhiE(pt,etaAK5PFPUcorrJet[j],phiAK5PFPUcorrJet[j],energy);
    }
    // offset-subtracted jets can have energy<0. In this case pT is fixed at 0. Do not consider these jets (to avoid warnings)
    if(energyAK5PFPUcorrJet[j] == 0) continue;

    float totenergy = (chargedHadronEnergyAK5PFPUcorrJet[j]+neutralHadronEnergyAK5PFPUcorrJet[j]+
                       chargedEmEnergyAK5PFPUcorrJet[j]+neutralEmEnergyAK5PFPUcorrJet[j]);
    float efrac = (chargedEmEnergyAK5PFPUcorrJet[j]+neutralEmEnergyAK5PFPUcorrJet[j]) / totenergy;
    float hfrac = (chargedHadronEnergyAK5PFPUcorrJet[j]+neutralHadronEnergyAK5PFPUcorrJet[j]) / totenergy;
    Jet sisconej(p4,efrac,hfrac);

    // PF jet ID variables
    float neutralHadFrac = neutralHadronEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    float neutralEmFraction = neutralEmEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    float nConstituents = chargedHadronMultiplicityAK5PFPUcorrJet[j] + neutralHadronMultiplicityAK5PFPUcorrJet[j] +
      photonMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j] +
      HFHadronMultiplicityAK5PFPUcorrJet[j] + HFEMMultiplicityAK5PFPUcorrJet[j];
    float chargedHadFraction = chargedHadronEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    float chargedMultiplicity = chargedHadronMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j];
    float chargedEmFraction = chargedEmEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    //sisconej.setPFNoPUJetID(neutralHadFrac,neutralEmFraction,nConstituents,
    //                    chargedHadFraction,chargedMultiplicity,chargedEmFraction);
    jets.push_back(sisconej);
  }

  return jets;

}

vector<BTagJet> Vecbos::GetBTagCorrPUPFJets() {
  
  vector<BTagJet> btags;
  btags.clear();
  for(int j=0; j<nAK5PFPUcorrJet; j++) {
    // offset-subtracted jets can have energy<0. In this case pT is fixed at 0. Do not consider these jets (to avoid warnings)
    if(energyAK5PFPUcorrJet[j] == 0) continue;
    BTagJet btag;
    btag.combinedSecondaryVertexBJetTags = combinedSecondaryVertexBJetTagsAK5PFPUcorrJet[j];
    //btag.combinedSecondaryVertexMVABJetTags = combinedSecondaryVertexMVABJetTagsAK5PFPUcorrJet[j];
    //btag.jetBProbabilityBJetTags = jetBProbabilityBJetTagsAK5PFPUcorrJet[j];
    //btag.jetProbabilityBJetTags = jetProbabilityBJetTagsAK5PFPUcorrJet[j];
    btag.simpleSecondaryVertexBJetTags = simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet[j];
    //btag.softMuonBJetTags = softMuonBJetTagsAK5PFPUcorrJet[j];
    btag.trackCountingHighPurBJetTags = trackCountingHighPurBJetTagsAK5PFPUcorrJet[j];
    btag.trackCountingHighEffBJetTags = trackCountingHighEffBJetTagsAK5PFPUcorrJet[j];
    btags.push_back(btag);
  }

  return btags;

}

#ifdef USECALOTOWERS

vector<float> Vecbos::DefaultCaloThresholds() {
    vector<float> thresh;
    thresh.reserve(10);
    thresh.push_back(.9);
    thresh.push_back(1.1);
    thresh.push_back(1.4);
    thresh.push_back(1.4);
    thresh.push_back(1.2);
    thresh.push_back(1.8);
    thresh.push_back(.09);
    thresh.push_back(.45);
    thresh.push_back(.2);
    thresh.push_back(.45);
    
    return thresh;
}

vector<CaloTower> Vecbos::CreateCaloTowers(vector<float> thresh, float zpv, int type){
  float HBthresh  = thresh[0];
  float HOthresh  = thresh[1];
  float HESthresh = thresh[2];
  float HEDthresh = thresh[3];
  float HF1thresh = thresh[4];
  float HF2thresh = thresh[5];
  float EBthresh  = thresh[6];
  float EEthresh  = thresh[7];
  float EBSumthresh = thresh[8];
  float EESumthresh = thresh[9];

  vector<CaloTower> CaloTowerCollection;
  float energy;
  float E = 0.0;
  float E_em = 0.0;
  float E_had = 0.0;
  int iECAL = 0;
  int iHCAL = 0;
  float energy_ecal[25];
  float eta_ecal[25];
  float phi_ecal[25];
  TVector3 ecal;
  TVector3 hcal;
  TVector3 calo;
  ecal.SetXYZ(0.0,0.0,0.0);
  hcal.SetXYZ(0.0,0.0,0.0);
  if(type == 0 || type == 4){
    calo.SetXYZ(xCaloTowers[0],
		yCaloTowers[0],
		zCaloTowers[0]-zpv);
  }
  if(type == 1 || type == 3){
    TVector3 dum(xCaloTowers[0],
		 yCaloTowers[0],
		 zCaloTowers[0]);
    CaloPointMy ref(0.0, dum.Eta(), CaloF);
    calo.SetXYZ(ref.r()*cos(dum.Phi()),
		ref.r()*sin(dum.Phi()),
		ref.z()-zpv);
  }
    
  float real_energy;
  bool HF = false;
  bool EB = false;
  bool EE = false;
  
  int index = 0;
  for(int i = 0; i < nCaloTowers; i++){
    energy = energyCaloTowers[i];
    if(CaloIndexCaloTowers[i] != index){
      if(EB){
	if(E_em >= EBSumthresh){
	  E += E_em;
	  if(type == 3 || type == 4){
	    double wEta = 0.0;
	    double wPhi = 0.0;
	    double wSum = 0.0;
	    for(int j = 0; j < iECAL; j++){
	      double w = 4.2 + log(energy_ecal[j]/E_em);
	      if(w < 0.0) w = 0.0;
	      wSum += w;
	      wEta += w*eta_ecal[j];
	      wPhi += w*phi_ecal[j];
	    }
	    wEta /= wSum;
	    wPhi /= wSum;
	    if(type == 3){
	      CaloPointMy re(0.0, wEta, CaloF);
	      calo.SetXYZ(re.r()*cos(wPhi),
			  re.r()*sin(wPhi),
			  re.z()-zpv);
	    }
	    if(type == 4){
	      calo.SetPtEtaPhi(calo.Pt(),wEta,wPhi);
	    }
	  }
	  
	} else {
	  E_em = 0.0;
	}
      }
      if(EE){
	if(E_em >= EESumthresh){
	  E += E_em;
	  if(type == 3 || type == 4){
	    double wEta = 0.0;
	    double wPhi = 0.0;
	    double wSum = 0.0;
	    for(int j = 0; j < iECAL; j++){
	      double w = 4.2 + log(energy_ecal[j]/E_em);
	      if(w < 0.0) w = 0.0;
	      wSum += w;
	      wEta += w*eta_ecal[j];
	      wPhi += w*phi_ecal[j];
	    }
	    wEta /= wSum;
	    wPhi /= wSum;
	    if(type == 3){
	      CaloPointMy re(0.0, wEta, CaloF);
	      calo.SetXYZ(re.r()*cos(wPhi),
			  re.r()*sin(wPhi),
			  re.z()-zpv);
	    }
	    if(type == 4){
	      calo.SetPtEtaPhi(calo.Pt(),wEta,wPhi);
	    }
	    
	  }
	} else {
	  E_em = 0.0;
	}
      }
      E += E_had;
      if(HF) E += E_em;
      if(E > 0.0){
	if(iECAL > 0) ecal = ecal*(1/double(iECAL));
	if(iHCAL > 0) hcal = hcal*(1/double(iHCAL));
	if(type != 2 && type != 5){
	  CaloTowerCollection.push_back(CaloTower(E_em, E_had, calo, ecal, hcal));
	}
	
      }
      E = 0.0;
      E_em = 0.0;
      E_had = 0.0;
      iECAL = 0;
      iHCAL = 0;
      ecal.SetXYZ(0.0,0.0,0.0);
      hcal.SetXYZ(0.0,0.0,0.0);
      if(type == 0 || type == 4){
	calo.SetXYZ(xCaloTowers[i],
		    yCaloTowers[i],
		    zCaloTowers[i]-zpv);
      }
      if(type == 1 || type == 3){
	TVector3 dumy(xCaloTowers[i],
		      yCaloTowers[i],
		      zCaloTowers[i]);
	CaloPointMy Ref(0.0, dumy.Eta(),CaloF);
	calo.SetXYZ(Ref.r()*cos(dumy.Phi()),
		    Ref.r()*sin(dumy.Phi()),
		    Ref.z()-zpv);
      }
      //calo.SetXYZ(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]);
    
      index = CaloIndexCaloTowers[i];
      HF = false;
      EB = false;
      EE = false;
    } else {
      if(CALOCaloTowers[i]/100 == 1){
	EB = true;
	if(energy >= EBthresh){
	  E_em += energy;
	  TVector3 p(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]);
	  ecal = ecal + p;
	  energy_ecal[iECAL] = energy;
	  eta_ecal[iECAL] = p.Eta();
	  phi_ecal[iECAL] = p.Phi();
	  iECAL++;
	  if(type == 2){
	    calo.SetXYZ(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]-zpv);
	    ecal.SetXYZ(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]-zpv);
	    hcal.SetXYZ(0.0,0.0,0.0);
	    CaloTowerCollection.push_back(CaloTower(energy, 0.0, calo, ecal, hcal));
	  }
	  
	}
      }
      if(CALOCaloTowers[i]/100 == 2){
	EE = true;
	if(energy >= EEthresh){
	  E_em += energy;
	  TVector3 p(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]);
	  ecal = ecal + p;
	  energy_ecal[iECAL] = energy;
	  eta_ecal[iECAL] = p.Eta();
	  phi_ecal[iECAL] = p.Phi();
	  iECAL++;
	  if(type == 2){
	    calo.SetXYZ(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]-zpv);
	    ecal.SetXYZ(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]-zpv);
	    hcal.SetXYZ(0.0,0.0,0.0);
	    CaloTowerCollection.push_back(CaloTower(energy, 0.0, calo, ecal, hcal));
	  }
	}
      }
      if((CALOCaloTowers[i]%100)/10 == 1){
	if(energy >= HBthresh){
	  E_had += energy;
	  iHCAL++;
	  TVector3 p(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]);
	  hcal = hcal + p;
	  if(type == 5){
	    calo.SetXYZ(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]-zpv);
	    hcal.SetXYZ(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]-zpv);
	    ecal.SetXYZ(0.0,0.0,0.0);
	    CaloTowerCollection.push_back(CaloTower(0.0, energy, calo, ecal, hcal));
	  }
	}
      }
      if((CALOCaloTowers[i]%100)/10 == 2){
	if(energy >= HESthresh){
	  E_had += energy;
	  iHCAL++;
	  TVector3 p(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]);
	  hcal = hcal + p;
	  if(type == 5){
	    calo.SetXYZ(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]-zpv);
	    hcal.SetXYZ(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]-zpv);
	    ecal.SetXYZ(0.0,0.0,0.0);
	    CaloTowerCollection.push_back(CaloTower(0.0, energy, calo, ecal, hcal));
	  }
	}
      }
      if((CALOCaloTowers[i]%100)/10 == 3){
	if(energy >= HEDthresh){
	  E_had += energy;
	  iHCAL++;
	  TVector3 p(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]);
	  hcal = hcal + p;
	  if(type == 5){
	    calo.SetXYZ(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]-zpv);
	    hcal.SetXYZ(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]-zpv);
	    ecal.SetXYZ(0.0,0.0,0.0);
	    CaloTowerCollection.push_back(CaloTower(0.0, energy, calo, ecal, hcal));
	  }
	}
      }
      if((CALOCaloTowers[i]%100)/10 == 4){
	if(energy >= HOthresh){
	  E_had += energy;
	  if(type == 5){
	    calo.SetXYZ(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]-zpv);
	    hcal.SetXYZ(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]-zpv);
	    ecal.SetXYZ(0.0,0.0,0.0);
	    CaloTowerCollection.push_back(CaloTower(0.0, energy, calo, ecal, hcal));
	  }
	}
      }
      if((CALOCaloTowers[i]%100)/10 == 5){
	HF = true;
	if(CALOCaloTowers[i]%10 == 1){
	  if(energy >= HF1thresh){
	    E_em += energy;
	    iHCAL++;
	    TVector3 p(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]);
	    hcal = hcal + p;
	  }
	} else {
	  if(energy >= HF2thresh){
	    E_em -= energy;
	    E_had += 2.0*energy;
	    iHCAL++;
	    TVector3 p(xCaloTowers[i],yCaloTowers[i],zCaloTowers[i]);
	    hcal = hcal + p;
	  }
	}
      }
    }
  }
  return CaloTowerCollection;
}

MET Vecbos::CreateMET(vector<CaloTower> InputCollection){
  double sum_et = 0.0;
  double sum_ex = 0.0;
  double sum_ey = 0.0;
  double sum_ez = 0.0;
  for(int i = 0; i < InputCollection.size(); i++){
    double phi = InputCollection[i].phi();
    double theta = InputCollection[i].theta();
    double e = InputCollection[i].e();
    double et = e*sin(theta);
    sum_ez += e*cos(theta);
    sum_et += et;
    sum_ex += et*cos(phi);
    sum_ey += et*sin(phi);
  }
  MET met(-sum_ex, -sum_ey, -sum_ez, sum_et);
  return met;
}

MET Vecbos::CreateMET(vector<Jet> InputCollection){
  double sum_et = 0.0;
  double sum_ex = 0.0;
  double sum_ey = 0.0;
  double sum_ez = 0.0;
  for(int i = 0; i < InputCollection.size(); i++){
    double phi = InputCollection[i].phi();
    double theta = InputCollection[i].theta();
    double e = InputCollection[i].e();
    double et = e*sin(theta);
    sum_ez += e*cos(theta);
    sum_et += et;
    sum_ex += et*cos(phi);
    sum_ey += et*sin(phi);
  }
  MET met(-sum_ex, -sum_ey, -sum_ez, sum_et);
  return met;
}
 
MET Vecbos::CreateMET(vector<TLorentzVector> InputCollection){
  double sum_et = 0.0;
  double sum_ex = 0.0;
  double sum_ey = 0.0;
  double sum_ez = 0.0;
  for(int i = 0; i < InputCollection.size(); i++){
    double phi = InputCollection[i].Phi();
    double theta = InputCollection[i].Theta();
    double e = InputCollection[i].E();
    double et = e*sin(theta);
    sum_ez += e*cos(theta);
    sum_et += et;
    sum_ex += et*cos(phi);
    sum_ey += et*sin(phi);
  }
  MET met(-sum_ex, -sum_ey, -sum_ez, sum_et);
  return met;
}
#endif

vector<int> Vecbos::AlgoBruteForce(int nMin, int nMax){
  vector<int> ca;
  vector<int> cb;
  vector<int> bestCB;
  float totalDeltaR=0;
  float BestTotalDeltaR=1000.0;
  
  for(int i1 = 0; i1 < nMax; i1++) ca.push_back(i1);
  for(int i1 = 0; i1 < nMin; i1++) cb.push_back(i1);

  for(int cnt = 0; cnt < TMath::Factorial(nMin); cnt++){
    totalDeltaR = lengtht(cb);
    if(totalDeltaR < BestTotalDeltaR){
      BestTotalDeltaR = totalDeltaR;
      bestCB=cb;
    }
    next_permutation(cb.begin(), cb.end());
  }
  while(next_combination(ca.begin(), ca.end(), cb.begin(), cb.end())){
    for(int cnt = 0; cnt < TMath::Factorial(nMin); cnt++){
      totalDeltaR = lengtht(cb);
      if(totalDeltaR < BestTotalDeltaR){
	BestTotalDeltaR = totalDeltaR;
	bestCB=cb;
      }
      
      next_permutation(cb.begin(), cb.end());
    }
  }
  
  return bestCB;
}

vector<int> Vecbos::AlgoSwitchMethod(int nMin, int nMax) {
  vector<int> bestCB;
  for(int i1 = 0; i1 < nMin; i1++){
    int minInd=0;
    for(int i2 = 1; i2 < nMax; i2++) {
      if(AllDist[i1][i2] < AllDist[i1][minInd])minInd = i2;
    }
    bestCB.push_back(minInd);
  }

  bool inside = true;
  while(inside){
    inside = false;
    for(int i1 = 0; i1 < nMin; i1++){
      for(int i2 = i1+1; i2 < nMin; i2++){
	if(bestCB[i1] == bestCB[i2]){
	  inside = true;
	  if(AllDist[i1][(bestCB[i1])] <= AllDist[i2][(bestCB[i2])]){
	    AllDist[i2][(bestCB[i2])] = 1000.0;
	    int minInd = 0;
	    for(int i3 = 1; i3 < nMax; i3++){
	      if(AllDist[i2][i3] < AllDist[i2][minInd]) minInd = i3;
	    }
	    bestCB[i2] = minInd;
	  } else{
	    AllDist[i1][(bestCB[i1])] = 1000.0;
	    int minInd = 0;
	    for(int i3 = 1; i3 < nMax; i3++){
	      if(AllDist[i1][i3] < AllDist[i1][minInd]) minInd = i3;
	    }
	    bestCB[i1] = minInd;
	  }
	}
      }
    }
  }
  
  return bestCB;
}
   
double Vecbos::lengtht(vector<int> best){
  double myLength = 0;
  int row = 0;
  for(vector<int>::iterator it=best.begin(); it != best.end(); it++){
    myLength+=AllDist[row][*it];
    row++;
  }
  return myLength;
}

vector<pair<Jet,Jet> > Vecbos::OneToOneMatch(vector<Jet> v1, vector<Jet> v2){
  vector<pair<Jet,Jet> > Matched;

  const int nSrc = v1.size();
  const int nMtc = v2.size();
  
  const int nMin = min(nSrc, nMtc);
  const int nMax = max(nSrc, nMtc);

  if(nMin < 1) return Matched;
  
  if( nSrc <= nMtc ){
    for(int iSr = 0; iSr < nSrc; iSr++){
      vector <float> tempAllDist;
      for(int iMt = 0; iMt < nMtc; iMt++){
	tempAllDist.push_back(DeltaR(v1[iSr].eta(), v1[iSr].phi(),
				     v2[iMt].eta(), v2[iMt].phi()));
      }
      AllDist.push_back(tempAllDist);
      tempAllDist.clear();
    }
  } else {
    for(int iMt = 0; iMt < nMtc; iMt++){
      vector <float> tempAllDist;
      for(int iSr = 0; iSr < nSrc; iSr++){
	tempAllDist.push_back(DeltaR(v1[iSr].eta(), v1[iSr].phi(),
				     v2[iMt].eta(), v2[iMt].phi()));
      }
      AllDist.push_back(tempAllDist);
      tempAllDist.clear();
    }
  }

  long nLoopToDo = (int)(TMath::Factorial(nMax) / TMath::Factorial(nMax-nMin));
  vector<int> bestCB;
  
  if(nLoopToDo < 10000 && nLoopToDo > 0) {
    bestCB = AlgoBruteForce(nMin,nMax);
  } else {
    bestCB = AlgoSwitchMethod(nMin,nMax);
  }

  // note: if v2 is sorted correctly, then the resulting collection
  // will remain sorted in _v2_
  
  if(nSrc <= nMtc){
    vector<pair<int, int> > sorter;
    for(int c = 0; c != nMin; c++){
      sorter.push_back(std::make_pair(bestCB[c],c));
    }
    std::sort(sorter.begin(), sorter.end());
    std::reverse(sorter.begin(), sorter.end());
    for(int c = 0; c != nMin; c++){
      Matched.push_back(std::make_pair(v1[sorter.at(c).second], 
				       v2[sorter.at(c).first]));
    }
  } else {
    for(int c = 0; c != nMin; c++){
      Matched.push_back(std::make_pair(v1[bestCB[c]], v2[c]));
    }
  }
  AllDist.clear();
  return Matched;

}

Jet Vecbos::CorrectJet(Jet J, double vtx){
  TLorentzVector v;
  CaloPointMy Ref(0.0, J.eta(),CaloF);
  v.SetPtEtaPhiE(J.p()/(cosh(Ref.etaReference(vtx))),
		 Ref.etaReference(vtx),
		 J.phi(),
		 J.e());
 Jet Jnew = Jet(v, J.EmFrac(), J.HadFrac());
  return Jnew;
}

vector<Jet> Vecbos::FastJetAlgorithm(vector<TLorentzVector> InputCollection, double Rparam, double thePtMin){
  string JetFinder = "kt_algorithm";
  string Strategy = "Best";
  double theDcut = -1;
  double theNjets = -1;
  string UE_subtraction = "no";
  bool theDoSubtraction = false;
  double theGhost_EtaMax = 6.0;
  int theActive_Area_Repeats = 5;
  double theGhostArea = 0.01;

  theJetConfig = new JetConfig;
  theJetConfig->theAreaSpec=fastjet::ActiveAreaSpec(theGhost_EtaMax, theActive_Area_Repeats, theGhostArea);

  fastjet::JetFinder jet_finder = fastjet::kt_algorithm;

  fastjet::Strategy strategy = fastjet::Best;

  /*
  int theMode = 0;
  
  if((theNjets!=-1)&&(theDcut==-1)){
    theMode = 3;
  } else if((theNjets==-1)&&(theDcut!=-1)){
    theMode = 2;
  } else if((theNjets!=-1)&&(theDcut!=-1)){
    theMode = 1;
  } else {
    theMode = 0;
  }
  */
  theJetConfig->theJetDef = fastjet::JetDefinition(jet_finder, Rparam, strategy);

  std::vector<fastjet::PseudoJet> input_vectors;
  int index_ = 0;
  for(unsigned int i = 0; i < InputCollection.size(); i++){
    double px = InputCollection[i].Px();
    double py = InputCollection[i].Py();
    double pz = InputCollection[i].Pz();
    double E = InputCollection[i].E();
    fastjet::PseudoJet PsJet(px,py,pz,E);
    PsJet.set_user_index(index_);
    input_vectors.push_back(PsJet);
    index_++;
  }

  vector<Jet> output;
  if(index_ == 0) return output;

  std::vector<fastjet::PseudoJet> theJets;

  // running without subtraction; need to add code for with subtraction
  fastjet::ClusterSequence clust_seq(input_vectors, theJetConfig->theJetDef);

  if((theNjets==-1)&&(theDcut==-1)){
    theJets=clust_seq.inclusive_jets(thePtMin);
  } else if((theNjets!=-1)&&(theDcut==-1)){
    theJets=clust_seq.exclusive_jets(theNjets);
  } else if((theNjets==-1)&&(theDcut!=-1)){
    theJets=clust_seq.exclusive_jets(theDcut);
  } else if((theNjets!=-1)&&(theDcut!=-1)){
    theJets=clust_seq.inclusive_jets(thePtMin);
  } else {
    theJets=clust_seq.inclusive_jets(thePtMin);
  }
  
  
  
  //here, for the reco jets, need to loop through constituents to get fractions
  for(std::vector<fastjet::PseudoJet>::const_iterator itJet=theJets.begin();
      itJet!=theJets.end();itJet++){
    
    double px = (*itJet).px();
    double py = (*itJet).py();
    double pz = (*itJet).pz();
    double E = (*itJet).E();
    TLorentzVector J(px,py,pz,E);
    output.push_back(Jet(J,0.0,0.0));
  }

  delete theJetConfig;
  return output;
}
  
vector<Jet> Vecbos::FastJetAlgorithm(vector<CaloTower> InputCollection, double Rparam, double thePtMin){
  string JetFinder = "kt_algorithm";
  string Strategy = "Best";
  double theDcut = -1;
  double theNjets = -1;
  string UE_subtraction = "no";
  bool theDoSubtraction = false;
  double theGhost_EtaMax = 6.0;
  int theActive_Area_Repeats = 5;
  double theGhostArea = 0.01;

  theJetConfig = new JetConfig;
  theJetConfig->theAreaSpec=fastjet::ActiveAreaSpec(theGhost_EtaMax, theActive_Area_Repeats, theGhostArea);

  fastjet::JetFinder jet_finder = fastjet::kt_algorithm;

  fastjet::Strategy strategy = fastjet::Best;
  /*
  int theMode = 0;
  
  if((theNjets!=-1)&&(theDcut==-1)){
    theMode = 3;
  } else if((theNjets==-1)&&(theDcut!=-1)){
    theMode = 2;
  } else if((theNjets!=-1)&&(theDcut!=-1)){
    theMode = 1;
  } else {
    theMode = 0;
  }
  */
  theJetConfig->theJetDef = fastjet::JetDefinition(jet_finder, Rparam, strategy);

  std::vector<fastjet::PseudoJet> input_vectors;
  int index_ = 0;
  for(int i = 0; i < InputCollection.size(); i++){
    if(InputCollection[i].et() > 0.5){
      double px = InputCollection[i].px();
      double py = InputCollection[i].py();
      double pz = InputCollection[i].pz();
      double E = InputCollection[i].e();
      fastjet::PseudoJet PsJet(px,py,pz,E);
      PsJet.set_user_index(index_);
      input_vectors.push_back(PsJet);
      index_++;
    }
  }

  vector<Jet> output;
  if(index_ == 0) return output;

  std::vector<fastjet::PseudoJet> theJets;

  // running without subtraction; need to add code for with subtraction
  fastjet::ClusterSequence clust_seq(input_vectors, theJetConfig->theJetDef);

  if((theNjets==-1)&&(theDcut==-1)){
    theJets=clust_seq.inclusive_jets(thePtMin);
  } else if((theNjets!=-1)&&(theDcut==-1)){
    theJets=clust_seq.exclusive_jets(theNjets);
  } else if((theNjets==-1)&&(theDcut!=-1)){
    theJets=clust_seq.exclusive_jets(theDcut);
  } else if((theNjets!=-1)&&(theDcut!=-1)){
    theJets=clust_seq.inclusive_jets(thePtMin);
  } else {
    theJets=clust_seq.inclusive_jets(thePtMin);
  }
  
  
  for(std::vector<fastjet::PseudoJet>::const_iterator itJet=theJets.begin();
      itJet!=theJets.end();itJet++){
    double EmFrac = 0.0;
    double HadFrac = 0.0;
    double total = 0.0;
    std::vector<fastjet::PseudoJet> jet_constituents = clust_seq.constituents(*itJet);
    for(std::vector<fastjet::PseudoJet>::const_iterator itConst=jet_constituents.begin();
	itConst!=jet_constituents.end(); itConst++){
      EmFrac += InputCollection[(*itConst).user_index()].EmEnergy();
      HadFrac += InputCollection[(*itConst).user_index()].HadEnergy();
    }
    double px = (*itJet).px();
    double py = (*itJet).py();
    double pz = (*itJet).pz();
    double E = (*itJet).E();
    TLorentzVector J(px,py,pz,E);
    total = EmFrac+HadFrac;
    if(total > 0.0){
      EmFrac = EmFrac/total;
      HadFrac = HadFrac/total;
    }
    output.push_back(Jet(J,EmFrac,HadFrac));
  }

  delete theJetConfig;
  return output;
}  

vector<Jet> Vecbos::SISCone(vector<CaloTower> InputCollection, double Rparam, double thePtMin){
  fastjet::SISConePlugin* mPlugin;

  fastjet::SISConePlugin::SplitMergeScale scale = fastjet::SISConePlugin::SM_pttilde;
  mPlugin = new fastjet::SISConePlugin(Rparam, 0.75, 0, thePtMin, false, scale);

  std::vector<fastjet::PseudoJet> input_vectors;

  int index_=0;
  for(int i = 0; i < InputCollection.size(); i++){
    if(InputCollection[i].et() > 0.5){
      double px = InputCollection[i].px();
      double py = InputCollection[i].py();
      double pz = InputCollection[i].pz();
      double E = InputCollection[i].e();
      fastjet::PseudoJet PsJet(px,py,pz,E);
      PsJet.set_user_index(index_);
      input_vectors.push_back(PsJet);
      index_++;
    }
  }

  vector<Jet> output;
  if(index_ == 0) return output;

  fastjet::JetDefinition jetDefinition(mPlugin);

  fastjet::ClusterSequence clusterSequence(input_vectors, jetDefinition);

  vector<fastjet::PseudoJet> inclusiveJets = clusterSequence.inclusive_jets(1.0);
  
  
  for(std::vector<fastjet::PseudoJet>::const_iterator itJet=inclusiveJets.begin();
      itJet!=inclusiveJets.end();itJet++){
    double EmFrac = 0.0;
    double HadFrac = 0.0;
    double total = 0.0;
    std::vector<fastjet::PseudoJet> jet_constituents = clusterSequence.constituents(*itJet);
    for(std::vector<fastjet::PseudoJet>::const_iterator itConst=jet_constituents.begin();
	itConst!=jet_constituents.end(); itConst++){
      EmFrac += InputCollection[(*itConst).user_index()].EmEnergy();
      HadFrac += InputCollection[(*itConst).user_index()].HadEnergy();
    }
    double px = (*itJet).px();
    double py = (*itJet).py();
    double pz = (*itJet).pz();
    double E = (*itJet).E();
    TLorentzVector J(px,py,pz,E);
    total = EmFrac+HadFrac;
    if(total > 0.0){
      EmFrac = EmFrac/total;
      HadFrac = HadFrac/total;
    }
    output.push_back(Jet(J,EmFrac,HadFrac));
  }
 

  delete mPlugin;
  return output;
}

vector<Jet> Vecbos::SISCone(vector<TLorentzVector> InputCollection, double Rparam, double thePtMin){
  fastjet::SISConePlugin* mPlugin;

  fastjet::SISConePlugin::SplitMergeScale scale = fastjet::SISConePlugin::SM_pttilde;
  mPlugin = new fastjet::SISConePlugin(Rparam, 0.75, 0, thePtMin, false, scale);

  std::vector<fastjet::PseudoJet> input_vectors;

  int index_=0;
  for(int i = 0; i < InputCollection.size(); i++){
    
    double px = InputCollection[i].Px();
    double py = InputCollection[i].Py();
    double pz = InputCollection[i].Pz();
    double E = InputCollection[i].E();
    fastjet::PseudoJet PsJet(px,py,pz,E);
    PsJet.set_user_index(index_);
    input_vectors.push_back(PsJet);
    index_++;
  }

  vector<Jet> output;
  if(index_ == 0) return output;

  fastjet::JetDefinition jetDefinition(mPlugin);

  fastjet::ClusterSequence clusterSequence(input_vectors, jetDefinition);

  vector<fastjet::PseudoJet> inclusiveJets = clusterSequence.inclusive_jets(1.0);
  
  
  for(std::vector<fastjet::PseudoJet>::const_iterator itJet=inclusiveJets.begin();
      itJet!=inclusiveJets.end();itJet++){
    
    double px = (*itJet).px();
    double py = (*itJet).py();
    double pz = (*itJet).pz();
    double E = (*itJet).E();
    TLorentzVector J(px,py,pz,E);
    
    output.push_back(Jet(J,0.0,0.0));
  }
 

  delete mPlugin;
  return output;
}

vector<Jet> Vecbos::CMSIterativeConeAlgorithm(vector<CaloTower> InputCollection, double Etmin, double seed){
  vector<CaloTower> input;
  vector<Jet> output;
  vector<pair<double,int> > Et;
  for(int i = 0; i < InputCollection.size(); i++){
    if(InputCollection[i].et() > Etmin){  
      Et.push_back(std::make_pair(InputCollection[i].et(), i));
    }
  }
  std::sort(Et.begin(), Et.end());
  std::reverse(Et.begin(), Et.end());
  for(int i = 0; i < Et.size(); i++){
    input.push_back(InputCollection[Et.at(i).second]);
  }
  
  while( !input.empty() && input[0].et() > seed){
    // for(int i = 0; i < input.size(); i++){
//       //cout << input[i].et() << endl;
//     }

    double eta0 = input[0].eta();
    double phi0 = input[0].phi();

    double eta=0;
    double phi=0;
    double et=0;

    vector<CaloTower> cone;
    vector<CaloTower> notcone;
    for(int iteration = 0; iteration<100; iteration++){
      cone.clear();
      notcone.clear();
      eta=0;
      phi=0;
      et=0;
      //vector<CaloTower>::iterator inp=input.begin();
      for(int itower = 0; itower < input.size(); itower++){
	double etat = input[itower].eta();
	double phit = input[itower].phi();
	double dphi = DeltaPhi(phi0, phit); 
	if(((etat-eta0)*(etat-eta0)+dphi*dphi) < (0.5*0.5)){
	  cone.push_back(input[itower]);
	  eta += input[itower].et()*input[itower].eta();
	  dphi = input[itower].phi()-phi0;
	  if(dphi>M_PI) dphi-=2*M_PI;
	  else if(dphi<=-M_PI) dphi+=2*M_PI;
	  phi+=input[itower].et()*dphi;
	  et+=input[itower].et();
	} else {
	  notcone.push_back(input[itower]);
	}
      }
      eta=eta/et;
      phi=phi0+phi/et;
      if(phi>M_PI)phi -= 2*M_PI;
      else if(phi<=-M_PI)phi+=2*M_PI;

      if((fabs(eta-eta0) < 0.001) && (fabs(phi-phi0) < 0.001)) break; //stable cone found
      eta0=eta;
      phi0=phi;
    }

    TLorentzVector J(0,0,0,0);
    
    double Em = 0.0;
    double Had = 0.0;
    for(int i = 0; i < cone.size(); i++){
      TLorentzVector v;
      v.SetPtEtaPhiE(cone[i].et(),
		     cone[i].eta(),
		     cone[i].phi(),
		     cone[i].e());
      J+=v;
      Em += cone[i].EmEnergy();
      Had += cone[i].HadEnergy();
    }
    if(cone.size()){ 
      //J = J * (1/weights);
      double total = Em+Had;
      if(total > 0.0){
	output.push_back(Jet(J,(Em/total),(Had/total)));
      } else {
	output.push_back(Jet(J,0.0,0.0));
      }
			 
    }
    input.clear();
    for(int i = 0; i < notcone.size(); i++){
      input.push_back(notcone[i]);
    }
  }
  return output;
}

vector<Jet> Vecbos::CMSIterativeConeAlgorithm(vector<TLorentzVector> InputCollection, double Etmin, double seed){
  vector<TLorentzVector> input;
  vector<Jet> output;
  vector<pair<double,int> > Et;
  for(int i = 0; i < InputCollection.size(); i++){
    if(InputCollection[i].Et() > Etmin){  
      Et.push_back(std::make_pair(InputCollection[i].Et(), i));
    }
  }
  std::sort(Et.begin(), Et.end());
  std::reverse(Et.begin(), Et.end());
  for(int i = 0; i < Et.size(); i++){
    input.push_back(InputCollection[Et.at(i).second]);
  }
  
  while( !input.empty() && input[0].Et() > seed){
    // for(int i = 0; i < input.size(); i++){
//       //cout << input[i].et() << endl;
//     }

    double eta0 = input[0].Eta();
    double phi0 = input[0].Phi();

    double eta=0;
    double phi=0;
    double et=0;

    vector<TLorentzVector> cone;
    vector<TLorentzVector> notcone;
    for(int iteration = 0; iteration<100; iteration++){
      cone.clear();
      notcone.clear();
      eta=0;
      phi=0;
      et=0;
      
      for(int itower = 0; itower < input.size(); itower++){
	double etat = input[itower].Eta();
	double phit = input[itower].Phi();
	double dphi = DeltaPhi(phi0, phit); 
	if(((etat-eta0)*(etat-eta0)+dphi*dphi) < (0.5*0.5)){
	  cone.push_back(input[itower]);
	  eta += input[itower].Et()*input[itower].Eta();
	  dphi = input[itower].Phi()-phi0;
	  if(dphi>M_PI) dphi-=2*M_PI;
	  else if(dphi<=-M_PI) dphi+=2*M_PI;
	  phi+=input[itower].Et()*dphi;
	  et+=input[itower].Et();
	} else {
	  notcone.push_back(input[itower]);
	}
      }
      eta=eta/et;
      phi=phi0+phi/et;
      if(phi>M_PI)phi -= 2*M_PI;
      else if(phi<=-M_PI)phi+=2*M_PI;

      if((fabs(eta-eta0) < 0.001) && (fabs(phi-phi0) < 0.001)) break; //stable cone found
      eta0=eta;
      phi0=phi;
    }

    TLorentzVector J(0,0,0,0);
    
    double Em = 0.0;
    double Had = 0.0;
    for(int i = 0; i < cone.size(); i++){
      TLorentzVector v;
      v.SetPtEtaPhiE(cone[i].Pt(),
		     cone[i].Eta(),
		     cone[i].Phi(),
		     cone[i].E());
      J+=v;
      //Em += cone[i].EmEnergy();
      //Had += cone[i].HadEnergy();
    }
    if(cone.size()){ 
      //J = J * (1/weights);
      //double total = Em+Had;
      output.push_back(Jet(J,0.0,0.0));
			 
    }
    input.clear();
    for(int i = 0; i < notcone.size(); i++){
      input.push_back(notcone[i]);
    }
  }
  return output;
}



double Vecbos::EventEMF(vector<Jet> vJ, double ptCut, double etaCut){
  int n = 0; //num of recjets satisfying cuts
  double sumWeightedJetPt = 0.;
  double sumJetPt = 0.;
  double eemf = -999.0;
  for(int i = 0; i < vJ.size(); i++){
    if(vJ[i].pt() > ptCut && fabs(vJ[i].eta()) < etaCut){
      ++n;
      sumJetPt += vJ[i].pt();
      sumWeightedJetPt += vJ[i].pt()*vJ[i].EmFrac();
    }
  }

  if(sumJetPt != 0.0) eemf = sumWeightedJetPt / sumJetPt;
  return eemf;
}

double Vecbos::JetCHF(Jet J, double dRcut, int& N){
  double SumPt = 0.;
  N = 0;
  
  
  double etaCorr = J.eta();
  TVector3 vPV(PVxPV[iPV],PVyPV[iPV],PVzPV[iPV]);
  for(int i = 0; i < nTrack; i++){
    TVector3 pT(pxTrack[i],pyTrack[i],pzTrack[i]);
    TVector3 vT(trackVxTrack[i],trackVyTrack[i],trackVzTrack[i]);
    vT = vT-vPV;
    if(fabs(vT.z()) > 0.1) continue;
    if(fabs(vT.Mag()) > 1.0) continue;
    //if(abs(chargeTrack[i]) > 1) continue;
    double pt = sqrt(pxTrack[i]*pxTrack[i]+pyTrack[i]*pyTrack[i]);
    if(pt < 1.2) continue;
    if(pt > 500.) continue;
    if(trackNormalizedChi2Track[i] > 20.0) continue;
    if(abs(transvImpactParTrack[i]/transvImpactParErrorTrack[i])>.6) continue;
    if(fabs(pT.Eta()) > 2.4) continue;
    if(trackValidHitsTrack[i] < 5) continue;
    
    double dR = DeltaR(etaCorr, J.phi(), double(pT.Eta()), double(pT.Phi()));
    if(dR <= dRcut){
      SumPt += pt;
      
      N++;
    }
  }

  return (SumPt/J.pt());
}

double Vecbos::EventCHF(vector<Jet> vJ, double ptCut, 
			    double etaCut){
  double echf = -999.0;
  int NJet = 0;
  double SumCHF = 0.;
  int NTrack;
   
  for(int i = 0; i < vJ.size(); i++){
    double etaCorr = vJ[i].eta();
    
    if(vJ[i].pt() < ptCut)continue;
    if(fabs(etaCorr) > etaCut) continue;

    
    //NJet++;
    double ThisJetCHF = JetCHF(vJ[i], 0.75, NTrack);
    if(NTrack < 4) continue;
    NJet++;
    SumCHF += ThisJetCHF;
  }
  if(NJet > 0){ 
    echf = SumCHF / double(NJet);
  }

  return echf;
}


bool Vecbos::AlpgenIdSelection(double alpgenid, string sample) {

  bool passed = false;
  
  // W+jets
  if(sample == "Wjets" && (alpgenid >=1000. && alpgenid <1999.))
    passed = true;

  // Z+jets
  if(sample == "Zjets" && (alpgenid >=2000. && alpgenid <2999.))
    passed = true;
  
  // ttbar+jets
  if(sample == "ttbarjets" && (alpgenid >=3000.&& alpgenid <3999.))
    passed = true;

  return passed;
}

vector<TLorentzVector> Vecbos::ParticlesFromMc(int status) {
    
  vector<TLorentzVector> theParticles;
  /*
    theParticles.reserve(nMc);
  
  for(int i = 0; i != nMc; ++i) 
    if(statusMc[i] == status) {
      TLorentzVector theParticle;
      theParticle.SetPtEtaPhiE(pMc[i]*sin(thetaMc[i]),
			       etaMc[i],
			       phiMc[i],
			       energyMc[i]);
      theParticles.push_back(theParticle);
    }
  */
  return theParticles;
}

vector<TLorentzVector> Vecbos::ParticlesFromMcWithId(int status, const vector<int>& allowed) {
  
  vector<TLorentzVector> theParticles;

  /*
  theParticles.reserve(nMc);
  bool isOk = false;
  
  for(int i = 0; i != nMc; ++i) {
    isOk = false;
    for(int j = 0; j != allowed.size(); ++j) 
      isOk = ( (statusMc[i] == status && idMc[i] == allowed.at(j)) || isOk );
    
    // isOk is the logical OR of: a. particle has correct status AND particle is allowed.
    // b. particle was Ok last time.
    
    if(isOk)
      {
	TLorentzVector theParticle;
	theParticle.SetPtEtaPhiE(pMc[i]*sin(thetaMc[i]),
				 etaMc[i],
				 phiMc[i],
				 energyMc[i]);
	theParticles.push_back(theParticle);
      }
  }
  */  
  return theParticles;
}

vector<TLorentzVector> Vecbos::ParticlesFromMcWithNotId(int status, const vector<int>& forbidden) {
  
  vector<TLorentzVector> theParticles;
  /*
  theParticles.reserve(nMc);
  bool isOk = true;
  
  for(int i = 0; i != nMc; ++i) {
    isOk = true;
    for(int j = 0; j != forbidden.size(); ++j)
      isOk = (statusMc[i] == status && idMc[i] != forbidden.at(j) && isOk);
    
    // isOk is the logical AND of: a. particle has correct status.
    // b. particle is not forbidden.
    // c. particle was Ok last time.
      
    if(isOk)
      {
	TLorentzVector theParticle;
	theParticle.SetPtEtaPhiE(pMc[i]*sin(thetaMc[i]),
				 etaMc[i],
				 phiMc[i],
				 energyMc[i]);
	theParticles.push_back(theParticle);
      }
  }
  */
  return theParticles;
}

vector<TLorentzVector> Vecbos::ParticlesFromMcCharged() {

  int c[10] = {11,-11,13,-13,211,-211,321,-321,2212,-2212};
  vector<int> charged(c, c + sizeof(c)/sizeof(c[0]));
  vector<TLorentzVector> chargedParticles = Vecbos::ParticlesFromMcWithId(1, charged);
  
  return chargedParticles;
}

vector<TLorentzVector> Vecbos::CloseInEtaPhi(const vector<TLorentzVector>& set, const TLorentzVector& v, double maxDistance) {
  
  double eta2 = v.Eta();
  double phi2 = v.Phi();

  vector<TLorentzVector> theCloseOnes;
  double distance = 0.;
  for(int i=0; i!=set.size(); ++i) {
    double eta1 = set.at(i).Eta();
    double phi1 = set.at(i).Phi();
    if(DeltaR(eta1, phi1, eta2, phi2) < maxDistance)
      theCloseOnes.push_back(set.at(i));
  }
  return theCloseOnes;
}

vector<TLorentzVector> Vecbos::Tracks(double thePtCut) {
  vector<TLorentzVector> theTracks;
  for(int i=0; i!=nTrack; i++) {
    double ptTrack = sqrt(pxTrack[i]*pxTrack[i] +
			  pyTrack[i]*pyTrack[i]);
    if(ptTrack >thePtCut) {
      TVector3 p3T(pxTrack[i],pyTrack[i],pzTrack[i]);
      TLorentzVector thisTrack(p3T,p3T.Mag());
      theTracks.push_back(thisTrack);
    }
  }
  return theTracks;
}

/// sigma ieta ieta of the seed cluster (ecal-driven/tracker-driven)
float Vecbos::SigmaiEiE(int electron) {
  float see;
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[electron], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[electron];
    see = sqrt(covIEtaIEtaSC[sc]);
  } else {
    int sc = PFsuperClusterIndexEle[electron];
    if(sc>-1) {
      see = sqrt(covIEtaIEtaPFSC[sc]);
    } else {
      see = 999.;
    }
  }
  return see;
}

/// sigma iphi iphi of the seed cluster (ecal-driven/tracker-driven)
float Vecbos::SigmaiPiP(int electron) {
  float see;
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[electron], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[electron];
    see = sqrt(covIPhiIPhiSC[sc]);
  } else {
    int sc = PFsuperClusterIndexEle[electron];
    if(sc>-1) {
      see = sqrt(covIPhiIPhiPFSC[sc]);
    } else {
      see = 999.;
    }
  }
  return see;
}

float Vecbos::likelihoodRatio(int eleIndex, ElectronLikelihood &lh) {
  LikelihoodMeasurements measurements;
  Utils anaUtils;
  bool inEB=anaUtils.fiducialFlagECAL(fiducialFlagsEle[eleIndex],isEB);
  measurements.pt = GetPt(pxEle[eleIndex],pyEle[eleIndex]);
  if(inEB && fabs(etaEle[eleIndex])<1.0) measurements.subdet = 0;
  else if (inEB && fabs(etaEle[eleIndex])>=1.0) measurements.subdet = 1;
  else measurements.subdet = 2;
  measurements.deltaPhi = deltaPhiAtVtxEle[eleIndex];
  measurements.deltaEta = deltaEtaAtVtxEle[eleIndex];
  measurements.eSeedClusterOverPout = eSeedOverPoutEle[eleIndex];
  measurements.eSuperClusterOverP = eSuperClusterOverPEle[eleIndex];
  measurements.hadronicOverEm = hOverEEle[eleIndex];
  measurements.sigmaIEtaIEta = SigmaiEiE(eleIndex);
  measurements.sigmaIPhiIPhi = SigmaiPiP(eleIndex);
  measurements.fBrem = fbremEle[eleIndex];
  measurements.nBremClusters = nbremsEle[eleIndex];
  int gsftrack = gsfTrackIndexEle[eleIndex];
  TVector3 pIn(pxGsfTrack[gsftrack],pyGsfTrack[gsftrack],pzGsfTrack[gsftrack]);
  measurements.OneOverEMinusOneOverP = 1./(eSuperClusterOverPEle[eleIndex]*pIn.Mag()) - 1./pIn.Mag();
  return lh.result(measurements);
}

float Vecbos::E9ESC(int electron) {
  float e9esc;
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[electron], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[electron];
    e9esc = e3x3SC[sc]/rawEnergySC[sc];
  } else {
    int sc = PFsuperClusterIndexEle[electron];
    if(sc>-1) {
      e9esc = e3x3PFSC[sc]/rawEnergyPFSC[sc];
    } else {
      e9esc = 999.;
    }
  }
  //  std::cout << "returning " <<e9esc << std::endl;
  return e9esc;


}

float Vecbos::E25ESC(int electron) {
  float e25esc;
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[electron], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[electron];
    e25esc = e5x5SC[sc]/rawEnergySC[sc];
  } else {
    int sc = PFsuperClusterIndexEle[electron];
    if(sc>-1) {
      e25esc = e5x5PFSC[sc]/rawEnergyPFSC[sc];
    } else {
      e25esc = 999.;
    }
  }
  //  std::cout << "returning " <<e25esc << std::endl;
  return e25esc;
}

float Vecbos::ESC(int electron) {
  float esc;
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[electron], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[electron];
    esc = energySC[sc];
  } else {
    int sc = PFsuperClusterIndexEle[electron];
    if(sc>-1) {
      esc = energyPFSC[sc];
    } else {
      esc = 999.;
    }
  }
  //  std::cout << "returning " <<esc << std::endl;
  return esc;

}

float Vecbos::Eseed(int electron) {
  float eseed;
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[electron], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[electron];
    eseed = eMaxSC[sc];
  } else {
    int sc = PFsuperClusterIndexEle[electron];
    if(sc>-1) {
      eseed = eMaxSC[sc];
    } else {
      eseed = 999.;
    }
  }
  //  std::cout << "returning " <<esc << std::endl;
  return eseed;

}

/// bremsstrahlung fraction
float Vecbos::FBrem(int electron) {
  /*
 int gsfTrack = gsfTrackIndexEle[electron];
  TVector3 p3OutEle(pxAtOuterGsfTrack[gsfTrack],pyAtOuterGsfTrack[gsfTrack],pzAtOuterGsfTrack[gsfTrack]);
  TVector3 p3InEle(pxGsfTrack[gsfTrack],pyGsfTrack[gsfTrack],pzGsfTrack[gsfTrack]);
  return (p3InEle.Mag() - p3OutEle.Mag())/p3InEle.Mag();
  */
  return -1;
}

// dxy parameter with respect to PV for electron tracks
double Vecbos::eleDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz) { 
  float elePt = sqrt(elePx*elePx + elePy*elePy);
  return ( - (eleVx-PVx)*elePy + (eleVy-PVy)*elePx ) / elePt; 
}

/// dsz parameter with respect to PV for electron tracks
double Vecbos::eleDszPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz) { 
  float elePt = sqrt(elePx*elePx + elePy*elePy);
  float eleP  = sqrt(elePx*elePx + elePy*elePy + elePz*elePz);
  return (eleVz-PVz)*elePt/eleP - ((eleVx-PVx)*elePx+(eleVy-PVy)*elePy)/elePt *elePz/eleP; 
}

/// returns the pt of the true photon in a photon + jet event (approx pthat)
float Vecbos::photonPt() {
  
  float pt = 90000.0;
  /*
 for(int imc=0; imc<50; imc++) {
    if(idMc[imc]==22) {
      pt = pMc[imc]*fabs(sin(thetaMc[imc]));
      break;
    }
  }
  */
  return pt;

}

void Vecbos::setRequiredTriggers(const std::vector<std::string>& reqTriggers) {
  requiredTriggers=reqTriggers;
}

bool Vecbos::hasPassedHLT() {
  Utils anaUtils;
  //std::cout << m_requiredTriggers.size() <<" " << std::endl;
  return anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);
}

bool Vecbos::triggerMatch(float eta, float phi, float Dr){
  //std::cout << "finding match for eta: "<< eta <<" phi: " << phi <<std::endl;
  bool match=false;
  for( int i=0; i<m_requiredTriggers.size(); i++ ) {//loop over require trigger paths
    int pathIndex=m_requiredTriggers[i];
 
    //std::cout << "testing trigger " << pathIndex << " with " << sizePassing[pathIndex] << " passing objects" << std::endl; 
    if( sizePassing[pathIndex]>  0 ) { //some object has passed the required trigger
      for(int np = 0; np < sizePassing[pathIndex]; np++ ){
	int iP = indexPassing[ indexPassingPerPath[pathIndex] +np];
	//std::cout << "passing object eta: " << triggerObsEta[iP] << " phi: " <<  triggerObsPhi[iP] << std::endl;
	if(DeltaR(eta, phi,triggerObsEta[iP],  triggerObsPhi[iP] ) < Dr){
	  match=true;
	  //std::cout << "MATCH!" <<std::endl;
	  break;
	}
      }            
    }
    if(match)  //it's enough if one path matches
      break;
  }
  return match;
}

vector<int> Vecbos::getHLTOutput() {
  Utils anaUtils;
  return anaUtils.getTriggers(m_requiredTriggers, firedTrg);
}

// Razor Variables

double Vecbos::CalcMR(TLorentzVector ja, TLorentzVector jb){
  double temp = (ja.P()*jb.Pz()-jb.P()*ja.Pz())*(ja.P()*jb.Pz()-jb.P()*ja.Pz());
  temp /= (ja.Pz()-jb.Pz())*(ja.Pz()-jb.Pz())-(ja.P()-jb.P())*(ja.P()-jb.P());

  temp = 2.*sqrt(temp);
  // protect against nan
  if(temp != temp) temp = -9999999.;

  return temp;
}

double Vecbos::CalcMRstar(TLorentzVector ja, TLorentzVector jb){
  double A = ja.P();
  double B = jb.P();
  double az = ja.Pz();
  double bz = jb.Pz();
  TVector3 jaT, jbT;
  jaT.SetXYZ(ja.Px(),ja.Py(),0.0);
  jbT.SetXYZ(jb.Px(),jb.Py(),0.0);

  double temp = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
                     (jbT.Dot(jbT)-jaT.Dot(jaT))*(jbT.Dot(jbT)-jaT.Dot(jaT))/(jaT+jbT).Mag2());

  return temp;
}

double Vecbos::CalcGammaMRstar(TLorentzVector ja, TLorentzVector jb){
  double A = ja.P();
  double B = jb.P();
  double az = ja.Pz();
  double bz = jb.Pz();
  TVector3 jaT, jbT;
  jaT.SetXYZ(ja.Px(),ja.Py(),0.0);
  jbT.SetXYZ(jb.Px(),jb.Py(),0.0);
  double ATBT = (jaT+jbT).Mag2();

  double temp = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
                     (jbT.Dot(jbT)-jaT.Dot(jaT))*(jbT.Dot(jbT)-jaT.Dot(jaT))/(jaT+jbT).Mag2());

  double mybeta = (jbT.Dot(jbT)-jaT.Dot(jaT))/
    sqrt(ATBT*((A+B)*(A+B)-(az+bz)*(az+bz)));

  double mygamma = 1./sqrt(1.-mybeta*mybeta);

  //gamma times MRstar                                                                                                                                                                              
  temp *= mygamma;

  return temp;
}

double Vecbos::CalcMRP(TLorentzVector ja, TLorentzVector jb, TVector3 met){

  double jaP = ja.Pt()*ja.Pt() +ja.Pz()*jb.Pz()-ja.P()*jb.P();
  double jbP = jb.Pt()*jb.Pt() +ja.Pz()*jb.Pz()-ja.P()*jb.P();
  jbP *= -1.;
  double den = sqrt((ja.P()-jb.P())*(ja.P()-jb.P())-(ja.Pz()-jb.Pz())*(ja.Pz()-jb.Pz()));

  jaP /= den;
  jbP /= den;

  double temp = jaP*met.Dot(jb.Vect())/met.Mag() + jbP*met.Dot(ja.Vect())/met.Mag();
  temp = temp*temp;

  den = (met.Dot(ja.Vect()+jb.Vect())/met.Mag())*(met.Dot(ja.Vect()+jb.Vect())/met.Mag())-(jaP-jbP)*(jaP-jbP);

  if(den <= 0.0) return -1.;

  temp /= den;
  temp = 2.*sqrt(temp);

  double bR = (jaP-jbP)/(met.Dot(ja.Vect()+jb.Vect())/met.Mag());
  double gR = 1./sqrt(1.-bR*bR);
  
  temp *= gR;

  // protect against nan
  if(temp != temp) temp = -9999999.;
  return temp;
}

double Vecbos::CalcMTR(TLorentzVector ja, TLorentzVector jb, TVector3 met){

  double temp = met.Mag()*(ja.Pt()+jb.Pt()) - met.Dot(ja.Vect()+jb.Vect());
  temp /= 2.;

  temp = sqrt(temp);

  return temp;
}

double Vecbos::GetDiJetMass(std::vector<Jet> jets) {
  if(jets.size()!=2) return -1.;
  TLorentzVector p1(jets[0].px(),jets[0].py(),jets[0].pz(),jets[0].e());
  TLorentzVector p2(jets[1].px(),jets[1].py(),jets[1].pz(),jets[1].e());
  return (p1+p2).M();
}

double Vecbos::GetDiJetPt(std::vector<Jet> jets) {
  if(jets.size()!=2) return -1.;
  TLorentzVector p1(jets[0].px(),jets[0].py(),jets[0].pz(),jets[0].e());
  TLorentzVector p2(jets[1].px(),jets[1].py(),jets[1].pz(),jets[1].e());
  return (p1+p2).Pt();
}

double Vecbos::GetDeltaEtaJet(std::vector<Jet> jets) {
  if(jets.size()!=2) return -1.;
  TVector3 p1(jets[0].px(),jets[0].py(),jets[0].pz());
  TVector3 p2(jets[1].px(),jets[1].py(),jets[1].pz());
  float eta1 = p1.Eta();
  float eta2 = p2.Eta();
  float deltaEta = eta1-eta2;
  return deltaEta;
}

double Vecbos::GetDeltaPhiJetMet(std::vector<Jet> jets, TVector3 met){ 
  if(jets.size()!=2) return -1.;
  TVector3 p1(jets[0].px(),jets[0].py(),jets[0].pz());
  return p1.DeltaPhi(met);
}

bool Vecbos::electronPassWP(int index, int WP)
{
   if(index >= nEle)   // asking for an electron that doesn't exist!?  Meow!!
      return false;
   
   // maybe we can change this to enum at some point
   if(WP != 95 && WP != 90 && WP != 85 && WP != 80 && WP != 70 && WP != 60)
      return false;

   bool IsBarrel = false;
   bool IsEndcap = false;

   int iSC = superClusterIndexEle[index];
   if(iSC < 0 || iSC >= nSC)   // something wrong!  electron don't have SC associated?
      return false;

   int iGsfTrack = gsfTrackIndexEle[index];
   if(iGsfTrack < 0 || iGsfTrack >= nGsfTrack)
      return false;
   if(expInnerLayersGsfTrack[iGsfTrack] > 0)
      return false;

   double eta = etaSC[iSC];
   if(fabs(eta) < 1.4442)
      IsBarrel = true;
   if(fabs(eta) > 1.566 && fabs(eta) < 2.5)
      IsEndcap = true;

   double TrackIsolation = dr03TkSumPtEle[index];
   double EcalIsolation = dr03EcalRecHitSumEtEle[index];
   double HcalIsolation = dr03HcalTowerSumEtEle[index];

   double PT = sqrt(pxEle[index] * pxEle[index] + pyEle[index] * pyEle[index]);

   double SigmaIEtaIEta = covIEtaIEtaSC[iSC];
   double DEta = fabs(deltaEtaAtVtxEle[index]);
   double DPhi = fabs(deltaPhiAtVtxEle[index]);
   double HOverE = hOverEEle[index];   // not used at the moment....

   double CombinedIsolation = 0;
   if(IsBarrel == true)
      CombinedIsolation = TrackIsolation + max(0.0, EcalIsolation - 1) + HcalIsolation;
   else
      CombinedIsolation = TrackIsolation + EcalIsolation + HcalIsolation;
   CombinedIsolation = CombinedIsolation - rhoFastjet * TMath::Pi() * 0.3 * 0.3;
   CombinedIsolation = CombinedIsolation / PT;

   double ConversionDistance = fabs(convDistEle[index]);
   double ConversionDCotTheta = fabs(convDcotEle[index]);

   switch(WP)
   {
   case 60:
      if(ConversionDistance < 0.02)     return false;
      if(ConversionDCotTheta < 0.02)    return false;
      if(IsBarrel == true)
      {
         if(CombinedIsolation > 0.016)  return false;
         if(SigmaIEtaIEta > 0.01)       return false;
         if(DPhi > 0.020)               return false;
         if(DEta > 0.004)               return false;
      }
      if(IsEndcap == true)
      {
         if(CombinedIsolation > 0.008)  return false;
         if(SigmaIEtaIEta > 0.031)      return false;
         if(DPhi > 0.021)               return false;
         if(DEta > 0.004)               return false;
      }
      break;
   case 70:
      if(ConversionDistance < 0.02)     return false;
      if(ConversionDCotTheta < 0.02)    return false;
      if(IsBarrel == true)
      {
         if(CombinedIsolation > 0.030)  return false;
         if(SigmaIEtaIEta > 0.01)       return false;
         if(DPhi > 0.020)               return false;
         if(DEta > 0.004)               return false;
      }
      if(IsEndcap == true)
      {
         if(CombinedIsolation > 0.016)  return false;
         if(SigmaIEtaIEta > 0.031)      return false;
         if(DPhi > 0.021)               return false;
         if(DEta > 0.005)               return false;
      }
      break;
   case 80:
      if(ConversionDistance < 0.02)     return false;
      if(ConversionDCotTheta < 0.02)    return false;
      if(IsBarrel == true)
      {
         if(CombinedIsolation > 0.040)  return false;
         if(SigmaIEtaIEta > 0.01)       return false;
         if(DPhi > 0.027)               return false;
         if(DEta > 0.005)               return false;
      }
      if(IsEndcap == true)
      {
         if(CombinedIsolation > 0.033)  return false;
         if(SigmaIEtaIEta > 0.031)      return false;
         if(DPhi > 0.021)               return false;
         if(DEta > 0.006)               return false;
      }
      break;
   case 85:
      if(ConversionDistance < 0.02)     return false;
      if(ConversionDCotTheta < 0.02)    return false;
      if(IsBarrel == true)
      {
         if(CombinedIsolation > 0.053)  return false;
         if(SigmaIEtaIEta > 0.01)       return false;
         if(DPhi > 0.039)               return false;
         if(DEta > 0.005)               return false;
      }
      if(IsEndcap == true)
      {
         if(CombinedIsolation > 0.042)  return false;
         if(SigmaIEtaIEta > 0.031)      return false;
         if(DPhi > 0.028)               return false;
         if(DEta > 0.007)               return false;
      }
      break;
   case 90:
      if(IsBarrel == true)
      {
         if(CombinedIsolation > 0.085)  return false;
         if(SigmaIEtaIEta > 0.01)       return false;
         if(DPhi > 0.071)               return false;
         if(DEta > 0.007)               return false;
      }
      if(IsEndcap == true)
      {
         if(CombinedIsolation > 0.051)  return false;
         if(SigmaIEtaIEta > 0.031)      return false;
         if(DPhi > 0.047)               return false;
         if(DEta > 0.011)               return false;
      }
      break;
   case 95:
      if(IsBarrel == true)
      {
         if(CombinedIsolation > 0.150)  return false;
         if(SigmaIEtaIEta > 0.012)      return false;
         if(DPhi > 0.800)               return false;
         if(DEta > 0.007)               return false;
      }
      if(IsEndcap == true)
      {
         if(CombinedIsolation > 0.100)  return false;
         if(SigmaIEtaIEta > 0.031)      return false;
         if(DPhi > 0.700)               return false;
         if(DEta > 0.011)               return false;
      }
      break;
   default:
      return false;
   }

   return true;
}

bool Vecbos::electronPassWP95(int index)
{
   return electronPassWP(index, 95);
}

bool Vecbos::electronPassWP90(int index)
{
   return electronPassWP(index, 90);
}

bool Vecbos::electronPassWP85(int index)
{
   return electronPassWP(index, 85);
}

bool Vecbos::electronPassWP80(int index)
{
   return electronPassWP(index, 80);
}

bool Vecbos::electronPassWP70(int index)
{
   return electronPassWP(index, 70);
}

bool Vecbos::electronPassWP60(int index)
{
   return electronPassWP(index, 60);
}

bool Vecbos::muonPassLoose(int index)
{
   if(index >= nMuon)   // imaginary muon always fail ID
      return false;

   Utils AnalysisUtilities;
   if(AnalysisUtilities.muonIdVal(muonIdMuon[index], bits::AllGlobalMuons) == false)
      return false;

   int iTrack = trackIndexMuon[index];
   if(iTrack < 0 || iTrack >= nTrack)   // something's wrong
      return false;

   if(numberOfValidStripTIBHitsTrack[iTrack]
      + numberOfValidStripTIDHitsTrack[iTrack]
      + numberOfValidStripTOBHitsTrack[iTrack]
      + numberOfValidStripTECHitsTrack[iTrack] <= 10)
      return false;

   return true;
}

bool Vecbos::muonPassTight(int index)
{
   if(index >= nMuon)
      return false;

   Utils AnalysisUtilities;
   if(AnalysisUtilities.muonIdVal(muonIdMuon[index], bits::AllGlobalMuons) == false)
      return false;
   if(AnalysisUtilities.muonIdVal(muonIdMuon[index], bits::GlobalMuonPromptTight) == false)
      return false;
   if(AnalysisUtilities.muonIdVal(muonIdMuon[index], bits::AllTrackerMuons) == false)
      return false;

   double PT = sqrt(pxMuon[index] * pxMuon[index] + pyMuon[index] * pyMuon[index]);

   double CombinedIsolation = emEt03Muon[index] + hadEt03Muon[index] + sumPt03Muon[index];
   CombinedIsolation = CombinedIsolation - rhoFastjet * TMath::Pi() * 0.3 * 0.3;
   CombinedIsolation = CombinedIsolation / PT;

   if(CombinedIsolation >= 0.15)
      return false;

   int iTrack = trackIndexMuon[index];
   if(iTrack < 0 || iTrack >= nTrack)   // something's wrong.  VERY wrong.
      return false;
   
   if(numberOfValidStripTIBHitsTrack[iTrack]
      + numberOfValidStripTIDHitsTrack[iTrack]
      + numberOfValidStripTOBHitsTrack[iTrack]
      + numberOfValidStripTECHitsTrack[iTrack] <= 10)
      return false;

   if(numberOfValidPixelBarrelHitsTrack[iTrack] + numberOfValidPixelEndcapHitsTrack[iTrack] == 0)
      return false;

   if(fabs(transvImpactParTrack[iTrack]) >= 0.2)
      return false;

   return true;
}

bool Vecbos::eventPassHcalFilter()
{
   if(fail2011Filter == true)
      return false;
   return true;
}

bool Vecbos::caloJetPassTCHEL(int index)
{
   if(index < 0 || index >= nAK5Jet)                     return false;
   if(trackCountingHighEffBJetTagsAK5Jet[index] < 1.7)   return false;
   return true;
}

bool Vecbos::caloJetPassTCHEM(int index)
{
   if(index < 0 || index >= nAK5Jet)                     return false;
   if(trackCountingHighEffBJetTagsAK5Jet[index] < 3.3)   return false;
   return true;
}

bool Vecbos::caloJetPassTCHET(int index)
{
   if(index < 0 || index >= nAK5Jet)                      return false;
   if(trackCountingHighEffBJetTagsAK5Jet[index] < 10.2)   return false;
   return true;
}

bool Vecbos::caloJetPassTCHPL(int index)
{
   if(index < 0 || index >= nAK5Jet)                      return false;
   if(trackCountingHighPurBJetTagsAK5Jet[index] < 1.19)   return false;
   return true;
}

bool Vecbos::caloJetPassTCHPM(int index)
{
   if(index < 0 || index >= nAK5Jet)                      return false;
   if(trackCountingHighPurBJetTagsAK5Jet[index] < 1.93)   return false;
   return true;
}

bool Vecbos::caloJetPassTCHPT(int index)
{
   if(index < 0 || index >= nAK5Jet)                      return false;
   if(trackCountingHighPurBJetTagsAK5Jet[index] < 3.41)   return false;
   return true;
}

bool Vecbos::caloJetPassSSVHEM(int index)
{
   if(index < 0 || index >= nAK5Jet)                              return false;
   if(simpleSecondaryVertexHighEffBJetTagsAK5Jet[index] < 1.74)   return false;
   return true;
}

bool Vecbos::caloJetPassSSVHET(int index)
{
   if(index < 0 || index >= nAK5Jet)                              return false;
   if(simpleSecondaryVertexHighEffBJetTagsAK5Jet[index] < 3.05)   return false;
   return true;
}

bool Vecbos::caloJetPassSSVHPT(int index)
{
   if(index < 0 || index >= nAK5Jet)                              return false;
   if(simpleSecondaryVertexHighPurBJetTagsAK5Jet[index] < 2.00)   return false;
   return true;
}

bool Vecbos::pfJetPassTCHEL(int index)
{
   if(index < 0 || index >= nAK5PFNoPUJet)                     return false;
   if(trackCountingHighEffBJetTagsAK5PFNoPUJet[index] < 1.7)   return false;
   return true;
}

bool Vecbos::pfJetPassTCHEM(int index)
{
   if(index < 0 || index >= nAK5PFNoPUJet)                     return false;
   if(trackCountingHighEffBJetTagsAK5PFNoPUJet[index] < 3.3)   return false;
   return true;
}

bool Vecbos::pfJetPassTCHET(int index)
{
   if(index < 0 || index >= nAK5PFNoPUJet)                      return false;
   if(trackCountingHighEffBJetTagsAK5PFNoPUJet[index] < 10.2)   return false;
   return true;
}

bool Vecbos::pfJetPassTCHPL(int index)
{
   if(index < 0 || index >= nAK5PFNoPUJet)                      return false;
   if(trackCountingHighPurBJetTagsAK5PFNoPUJet[index] < 1.19)   return false;
   return true;
}

bool Vecbos::pfJetPassTCHPM(int index)
{
   if(index < 0 || index >= nAK5PFNoPUJet)                      return false;
   if(trackCountingHighPurBJetTagsAK5PFNoPUJet[index] < 1.93)   return false;
   return true;
}

bool Vecbos::pfJetPassTCHPT(int index)
{
   if(index < 0 || index >= nAK5PFNoPUJet)                      return false;
   if(trackCountingHighPurBJetTagsAK5PFNoPUJet[index] < 3.41)   return false;
   return true;
}

bool Vecbos::pfJetPassSSVHEM(int index)
{
   if(index < 0 || index >= nAK5PFNoPUJet)                              return false;
   if(simpleSecondaryVertexHighEffBJetTagsAK5PFNoPUJet[index] < 1.74)   return false;
   return true;
}

bool Vecbos::pfJetPassSSVHET(int index)
{
   if(index < 0 || index >= nAK5PFNoPUJet)                              return false;
   if(simpleSecondaryVertexHighEffBJetTagsAK5PFNoPUJet[index] < 3.05)   return false;
   return true;
}

bool Vecbos::pfJetPassSSVHPT(int index)
{
   if(index < 0 || index >= nAK5PFNoPUJet)                              return false;
   if(simpleSecondaryVertexHighPurBJetTagsAK5PFNoPUJet[index] < 2.00)   return false;
   return true;
}




