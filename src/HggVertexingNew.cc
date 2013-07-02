// std includes
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// ROOT includes
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLorentzVector.h>

// local includes
#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CoolTools.hh"
#include "CaloTower.hh"
#include "ReadConfig.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"
#include "HggVertexingNew.hh"

#define debugVertexing 0

HggVertexing::HggVertexing(VecbosBase *b):
  base(b)
{
  useConversion = true;
  isInit=false;
  doSaveInputs=false;
  PerVtxVars = new float[5];
  PerEvtVars = new float[8];
  rescaleTrkPt = false;
}

 
HggVertexing::~HggVertexing() {
  delete PerVtxVars;
  delete PerEvtVars;
}

void HggVertexing::init(){
  //read the config file to setup TMVA
  cout << "Doing Vertexing Initialization" << endl;
  ReadConfig cfg;
  int errorCode = cfg.read(configFilePath);
  cout << "ReadConfig exited with code: " << errorCode << endl;
  if(errorCode) return;
  cfg.printAll();
  perVtxMvaWeights = cfg.getParameter("perVtxMvaWeights");
  perVtxMvaMethod  = cfg.getParameter("perVtxMvaMethod");
  perEvtMvaWeights = cfg.getParameter("perEvtMvaWeights");
  perEvtMvaMethod  = cfg.getParameter("perEvtMvaMethod");
  rescaleTrkPt   = cfg.getParameter("doRescaleTrkPt").compare("yes")==0;

  perVtxReader = new TMVA::Reader( "!Color:!Silent" );
  perEvtReader = new TMVA::Reader( "!Color:!Silent" );

  perVtxReader->AddVariable( "ptbal", &PerVtxVars[0] );
  perVtxReader->AddVariable( "ptasym", &PerVtxVars[1] );
  perVtxReader->AddVariable( "logsumpt2", &PerVtxVars[2] );
  perVtxReader->AddVariable( "limPullToConv", &PerVtxVars[3] );
  perVtxReader->AddVariable( "nConv", &PerVtxVars[4] );
  perVtxReader->BookMVA( perVtxMvaMethod, perVtxMvaWeights);

  perEvtReader->AddVariable( "diphoPt0", &PerEvtVars[0] );
  perEvtReader->AddVariable( "nVert", &PerEvtVars[1] );
  perEvtReader->AddVariable( "MVA0", &PerEvtVars[2] );
  perEvtReader->AddVariable( "MVA1", &PerEvtVars[3] );
  perEvtReader->AddVariable( "dZ1", &PerEvtVars[4] );
  perEvtReader->AddVariable( "MVA2", &PerEvtVars[5] );
  perEvtReader->AddVariable( "dZ2", &PerEvtVars[6] );
  perEvtReader->AddVariable( "nConv", &PerEvtVars[7] );
  perEvtReader->BookMVA(perEvtMvaMethod, perEvtMvaWeights);
  
  isInit=true;
}

void HggVertexing::saveInputs(TTree* outputTree){
  doSaveInputs = true;

  outputTree->Branch( "ptbal", &PerVtxVars[0] );
  outputTree->Branch( "ptasym", &PerVtxVars[1] );
  outputTree->Branch( "logsumpt2", &PerVtxVars[2] );
  outputTree->Branch( "limPullToConv", &PerVtxVars[3] );
  outputTree->Branch( "nConv", &PerVtxVars[4] );

}

float HggVertexing::evalPerVtxMVA(float ptbal, float ptasym, float logsumpt2, float limPullToConv, float nConv){
  if(debugVertexing) cout << "evalPerVtxMVA" << endl;
  PerVtxVars[0] = ptbal;
  PerVtxVars[1] = ptasym;
  PerVtxVars[2] = logsumpt2;
  PerVtxVars[3] = limPullToConv;
  PerVtxVars[4] = nConv;
  if(debugVertexing){
    cout << "\tInputs:  ";
    for(int i=0;i<5;i++) cout << PerVtxVars[i] << "  ";
    cout << endl;
  }

  return perVtxReader->EvaluateMVA(perVtxMvaMethod);
}
bool HggVertexing::isGoodVertex(int i){
  if(i<0 || i >= base->nPV) return false;
  if(fabs(base->PVzPV[i]) > 24.) return false;
  if(base->trackSizePV[i] <=0) return false;
  if(base->ndofPV[i] < 4.0)    return false;
  if(base->rhoPV[i] > 2.0) return false;
  return true;
}

std::vector<std::pair<int,float> > HggVertexing::evalPerVtxMVA(VecbosPho* pho1, VecbosPho* pho2){
  if(debugVertexing) cout << "evalPerVtxMVA" << endl;
  double *ptbal          = new double[base->nPV];
  double *ptasym         = new double[base->nPV];
  double *sumpt2         = new double[base->nPV];
  double *limPullToConv  = new double[base->nPV];

  //lists of the 4 vectors from each vertex
  std::vector<TLorentzVector> pho1_fromVtx;
  std::vector<TLorentzVector> pho2_fromVtx;
  //loop over the vertices
  for(int iVtx=0; iVtx<base->nPV; iVtx++){
    TVector3 thisVtxPos(base->PVxPV[iVtx],base->PVyPV[iVtx],base->PVzPV[iVtx]);

    TLorentzVector pho1_p4 = pho1->p4FromVtx(thisVtxPos,pho1->finalEnergy);
    TLorentzVector pho2_p4 = pho2->p4FromVtx(thisVtxPos,pho2->finalEnergy);

    pho1_fromVtx.push_back(pho1_p4);
    pho2_fromVtx.push_back(pho2_p4);
    
    //clear
    ptbal      [iVtx] = 0.0;
    ptasym     [iVtx] = 0.0;
    sumpt2     [iVtx] = 0.0;
    limPullToConv [iVtx] = -1.0; 

  }//for(int iVtx=0; iVtx<base->nPV; iVtx++)

  if(debugVertexing) std::cout << "initial loop over vertices done" << std::endl;

  std::vector<TLorentzVector> trackMomentum;
  for(int iVtx=0; iVtx<base->nPV; iVtx++) trackMomentum.push_back(TLorentzVector(0.,0.,0.,0.));
  if(debugVertexing) std::cout << "initial loop over vertices done" << std::endl;

  for(int iTrk=0; iTrk<base->nTrack; iTrk++){
    int iVtx = base->vtxIndexTrack[iTrk];

    /*    
    if(iVtx==-1){ //try to manually match the track
      float bestDR=1e6;
      int  bestIndex=-1;
      TVector3 tkVtxPos(base->trackVxTrack[iTrk], base->trackVyTrack[iTrk], base->trackVzTrack[iTrk]);
      
      for(int i=0; i<base->nPV; i++){
	TVector3 vtxPos(base->PVxPV[i],base->PVyPV[i],base->PVzPV[i]);
	if(tkVtxPos.DeltaR(vtxPos) < 0.03 && tkVtxPos.DeltaR(vtxPos) < bestDR){
	//if(fabs(tkVtxPos.Z()-vtxPos.Z()) < 0.001 && fabs(tkVtxPos.Z()-vtxPos.Z()) < bestDR){
	  bestDR = tkVtxPos.DeltaR(vtxPos);//fabs(tkVtxPos.Z()-vtxPos.Z());
	  bestIndex = i;
	}
      }
      iVtx = bestIndex;
    }
    */
    if(iVtx==-1) continue;
    if(!isGoodVertex(iVtx)) continue;
    TVector3 trackP(base->pxTrack[iTrk],base->pyTrack[iTrk],base->pzTrack[iTrk]);
    TLorentzVector thisMomentum; //use M=0
    thisMomentum.SetVectM(trackP,0);
    /*
    thisMomentum.SetPxPyPzE(base->pxTrack[iTrk],base->pyTrack[iTrk],base->pzTrack[iTrk],
			    TMath::Sqrt(TMath::Power(base->pxTrack[iTrk],2)+
					TMath::Power(base->pyTrack[iTrk],2)+
					TMath::Power(base->pzTrack[iTrk],2)));
    */
    if(thisMomentum.Pt()<1e-6) continue;

    if(thisMomentum.Pt()==0) std::cout << "WARNING: INVALID PT in HggVertexingNew" <<std::endl;

    if(rescaleTrkPt){
      float modpt = (thisMomentum.Pt() > base->ptErrorTrack[iTrk] ? thisMomentum.Pt() - base->ptErrorTrack[iTrk] : 0. );
      if(modpt > 0){
	float corr = modpt/thisMomentum.Pt();
	thisMomentum*=corr;
      }
    }

    sumpt2[iVtx]+= thisMomentum.Pt()*thisMomentum.Pt();

    //if(debugVertexing) cout << "This Track " << iTrk << "  Pt: " << thisMomentum.Pt() << endl;

    if(fabs(pho1_fromVtx[iVtx].Eta()) > 10) std::cout << "ERROR IN VERTEXING: Photon 1" <<std::endl;
    if(fabs(pho2_fromVtx[iVtx].Eta()) > 10) std::cout << "ERROR IN VERTEXING: Photon 2" <<std::endl;
    if(fabs(thisMomentum.Eta()) > 10) std::cout << "ERROR IN VERTEXING: Track" <<std::endl;
    
    if(thisMomentum.DeltaR(pho1_fromVtx[iVtx]) < 0.05 ||
    thisMomentum.DeltaR(pho2_fromVtx[iVtx]) < 0.05 ) continue;
    
    trackMomentum[iVtx] += thisMomentum;
    
  }//for(int iTrk=0; iTrk<base->nTrack; iTrk++)

  std::pair<float,float> convZ = getZConv(pho1,pho2); //z,dz
  if(debugVertexing) cout << "conv z/dz: " << convZ.first << "  " << convZ.second << endl;
  int maxVertices = ( (pho1_fromVtx[0] + pho2_fromVtx[0]).Pt() > 30 ? 3 : 5);
  double minDz = 999;

  int maxMVAIndex=-1;
  float maxMVA=-99;
  std::vector<std::pair<int,float> > vertexMVAs;
  for(int iVtx=0; iVtx<base->nPV; iVtx++){ //loop over the vertices AGAIN to compute
    //this is inelegant, but far quicker, rather than looping over the track collection for
    //every vertex
    if(!isGoodVertex(iVtx)){
    vertexMVAs.push_back(std::pair<int,float>(iVtx,-2));
      continue;
    }
    TLorentzVector thisTrackMom = trackMomentum[iVtx];
    TLorentzVector thisHiggsMom = pho1_fromVtx[iVtx] + pho2_fromVtx[iVtx];

    ptbal[iVtx] =  thisTrackMom.Px()*thisHiggsMom.Px() + thisTrackMom.Py()*thisHiggsMom.Py();
    ptbal[iVtx] = -ptbal[iVtx]/thisHiggsMom.Pt();

    ptasym[iVtx] = (thisTrackMom.Pt() - thisHiggsMom.Pt())/(thisTrackMom.Pt() + thisHiggsMom.Pt());

    //if there is a conversion, compute the pull
    if(pho1->conversion.index != -1 || pho2->conversion.index != -1){
      if(debugVertexing) cout << "Vtx " << iVtx << " Z: " << base->PVzPV[iVtx] << "\t conv Z:"
			      << convZ.first << " +- " << convZ.second << endl;
      limPullToConv[iVtx] = TMath::Abs(base->PVzPV[iVtx]-convZ.first)/convZ.second;
    }
    
    if(debugVertexing) cout << "Vtx " << iVtx << ":  " << endl;
    double mva = evalPerVtxMVA(ptbal[iVtx],ptasym[iVtx],log(sumpt2[iVtx]),
			       limPullToConv[iVtx],getNConv(pho1,pho2));
    if(debugVertexing) cout << "Output: " << mva << endl;
    vertexMVAs.push_back(std::pair<int,float>(iVtx,mva));
    if(mva>maxMVA){
      maxMVA=mva;
      maxMVAIndex=iVtx;
    }
  }

  //hmmm std::sort doesn't appear to work
  //std::sort(vertexMVAs.begin(),vertexMVAs.end(),sort_pred()); //sort by MVA value

  //leave the PerVtxVariables initialized to the values of the best vertex (for filling)
  double mva = evalPerVtxMVA(ptbal[maxMVAIndex],ptasym[maxMVAIndex],log(sumpt2[maxMVAIndex]),
			     limPullToConv[maxMVAIndex],getNConv(pho1,pho2));

  delete ptbal;
  delete ptasym;
  delete sumpt2;
  delete limPullToConv;

  if(debugVertexing) cout << "done with perVtxMVA" << endl;
  return vertexMVAs;
}


float HggVertexing::evalPerEvtMVA(VecbosPho* pho1,VecbosPho* pho2, 
				  std::vector<std::pair<int,float> > *perVertexRank){
  if(debugVertexing) cout << "PerEvtMVA" << endl;
  int bestVtx = perVertexRank->at(0).first;
  if(bestVtx <0 || bestVtx >=base->nPV) return -1e6;
  TVector3 bestVtxPos(base->PVxPV[bestVtx],base->PVyPV[bestVtx],base->PVzPV[bestVtx]);
    
  TLorentzVector pho1_p4 = pho1->p4FromVtx(bestVtxPos,pho1->finalEnergy);
  TLorentzVector pho2_p4 = pho2->p4FromVtx(bestVtxPos,pho2->finalEnergy);
  TLorentzVector higgs_p4 = pho1_p4+pho2_p4;

  if(perVertexRank->size()==0) return -1e6;
  while(perVertexRank->size()<3) perVertexRank->push_back(std::pair<int,float>(-1,-1e6));


      
  PerEvtVars[0] = higgs_p4.Pt();
  PerEvtVars[1] = base->nPV;
  PerEvtVars[2] = perVertexRank->at(0).second;
  PerEvtVars[3] = perVertexRank->at(1).second;
  PerEvtVars[4] = (perVertexRank->at(1).first!=-1 ? base->PVzPV[perVertexRank->at(1).first] - bestVtxPos.Z() : 0);
  PerEvtVars[5] = perVertexRank->at(2).second;
  PerEvtVars[6] = (perVertexRank->at(2).first!=-1 ? base->PVzPV[perVertexRank->at(2).first] - bestVtxPos.Z() : 0);
  PerEvtVars[7] = getNConv(pho1,pho2);


  if(debugVertexing){
    cout << "\tInputs:  ";
    for(int i=0;i<8;i++) cout << PerEvtVars[i] << "  ";
    cout << endl;
  }
  if(debugVertexing) cout << "PerEvtMVA --Done" << endl;
  
  return perEvtReader->EvaluateMVA( perEvtMvaMethod );
}



//method to do the vertexing
vector<pair<int,float> > HggVertexing::vertex_tmva(VecbosPho *pho1, VecbosPho *pho2,float& evtMVA){ 
  if(debugVertexing) cout << "vertexTMVA  " <<base->runNumber << "  " 
			  << base->lumiBlock << "  " << base->eventNumber << endl;
 
  //if(base->nPV == 1) return pair<int,float>(0,1); 
  if(base->nPV == 0) return std::vector<pair<int,float> >();

  std::vector<std::pair<int,float> > vertexRanks = evalPerVtxMVA(pho1,pho2);
  evtMVA = evalPerEvtMVA(pho1,pho2,&vertexRanks);
  if(debugVertexing) cout << "Done --VertexTMVA" << endl;

  return vertexRanks;
}

std::pair<float,float> HggVertexing::getZConv(VecbosPho* pho1,VecbosPho* pho2){

  std::pair<float,float> pho1_z = getZConv(pho1);
  std::pair<float,float> pho2_z = getZConv(pho2);

  float z1 = pho1_z.first;
  float dz1 = pho1_z.second;
  float z2 = pho2_z.first;
  float dz2 = pho2_z.second;

  float zconv=0;
  float dzconv=0;
  if(pho1->conversion.index!=-1){
    zconv = z1; dzconv = dz1;
    if(debugVertexing) cout << "Photon 1 has conversion   z1: " << z1 << endl;
  }
  if(pho2->conversion.index!=-1){
    zconv = z2; dzconv = dz2;
    if(debugVertexing) cout << "Photon 2 has conversion   z2: " << z2 << endl;
  }
  if(pho1->conversion.index!=-1 && pho2->conversion.index!=-1){
    float zconv = (z1/dz1/dz1 + z2/dz2/dz2)/(1./dz1/dz1 + 1./dz2/dz2 );
    float dzconv = sqrt( 1./(1./dz1/dz1 + 1./dz2/dz2)) ;
  }
  return std::pair<float,float>(zconv,dzconv);
}

//copied from MIT framework
 std::pair<float,float> HggVertexing::getZConv(VecbosPho* pho){
  if(debugVertexing) std::cout << "getZConv" << std::endl;
   const float dzpxb = 0.016;
   const float dztib = 0.331;
   const float dztob = 1.564;
   const float dzpxf = 0.082;
   const float dztid = 0.321;
   const float dztec = 0.815;

   const float dzpxbsingle = 0.036;
   const float dztibsingle = 0.456;
   const float dztobsingle = 0.362;
   const float dzpxfsingle = 0.130;
   const float dztidsingle = 0.465;
   const float dztecsingle = 1.018;

   float zconv  = -99;
   float dzconv = -99;

   VecbosConversion c = pho->conversion;
   if(c.index == -1) return std::pair<float,float>(0.,0.);


  TVector3 beamSpot(base->beamSpotX,base->beamSpotY,base->beamSpotZ);
  float zconvtrk = convCorrectedDz(&c,beamSpot) + beamSpot.Z();
  float zconvsc  = Z0EcalVtxCiC(&c,beamSpot,pho->SC.CaloPos);

  if(debugVertexing) cout << "zconv SC: " << zconvsc << endl;
  if(debugVertexing) cout << "zconv TRK: " << zconvtrk << endl;

   if( c.vtxNTracks == 2){
     if(pho->SC.CaloPos.Eta() < 1.5){
       float rho = c.CaloPos.Mag();
       if(rho<15)           { dzconv = dzpxb; zconv = zconvtrk; }
       else if( rho < 60. ) { dzconv = dztib; zconv = zconvsc; }
       else                 { dzconv = dztob; zconv = zconvsc; }
     }else{ 
       float z = c.CaloPos.Z();
       if     ( TMath::Abs(z) < 50. ) { dzconv = dzpxf; zconv = zconvtrk; }
       else if( TMath::Abs(z) < 100.) { dzconv = dztid; zconv = zconvtrk; }
       else                           { dzconv = dztec; zconv = zconvsc; }
     }
   }//if( c.vtxNTrack == 2)
   else if( c.vtxNTracks == 1){
     if(pho->SC.CaloPos.Eta() < 1.5){
       float rho = c.CaloPos.Mag();
       if(rho<15)           { dzconv = dzpxbsingle; zconv = zconvsc; }
       else if( rho < 60. ) { dzconv = dztibsingle; zconv = zconvsc; }
       else                 { dzconv = dztobsingle; zconv = zconvsc; }
     }else{
       float z = c.CaloPos.Z();
       if     ( TMath::Abs(z) < 50. ) { dzconv = dzpxfsingle; zconv = zconvsc; }
       else if( TMath::Abs(z) < 100.) { dzconv = dztidsingle; zconv = zconvsc; }
       else                           { dzconv = dztecsingle; zconv = zconvsc; }
     }     
   }

   return std::pair<float,float>(zconv,dzconv);

}

/*  
//MIT-like version
 float HggVertexing::convCorrectedDz(VecbosConversion* c,TVector3 basePos){
   TVector3 momPerp(c->pRefittedPair.Px(),c->pRefittedPair.Py(),0);
   TVector3 posPerp(c->CaloPos.X()-basePos.X(),c->CaloPos.Y()-basePos.Y(),0);

   return c->CaloPos.Z() - basePos.Z() - posPerp.Dot(momPerp)/momPerp.Pt() * (c->pRefittedPair.Pz()/momPerp.Pt());
0 }

float HggVertexing::Z0EcalVtxCiC(VecbosConversion* c, TVector3 basePos,TVector3 caloPos){
   TVector3 dirscvtx = caloPos - c->CaloPos;
   TVector3 momPerp(c->pRefittedPair.Px(),c->pRefittedPair.Py(),0);
   TVector3 posPerp(c->CaloPos.X()-basePos.X(),c->CaloPos.Y()-basePos.Y(),0);

   return c->CaloPos.Z() - posPerp.Mag() * (dirscvtx.Z()/dirscvtx.Mag());
   
 }
*/

//Globe-like version
float HggVertexing::convCorrectedDz(VecbosConversion* c,TVector3 basePos){
  if(debugVertexing) std::cout << "convCorrectedDz" << std::endl;
  return (c->vtx.Z()-basePos.Z()) - ((c->vtx.X()-basePos.X())*c->pRefittedPair.Px() + (c->vtx.Y()-basePos.Y())*c->pRefittedPair.Py())/c->pRefittedPair.Pt()*c->pRefittedPair.Pz()/c->pRefittedPair.Pt();

}

float HggVertexing::Z0EcalVtxCiC(VecbosConversion* c, TVector3 basePos,TVector3 caloPos){
  if(debugVertexing) std::cout << "Z0EcalVtxCiC" << std::endl;
  float dx1 = caloPos.X() - c->vtx.X();
  float dy1 = caloPos.Y() - c->vtx.Y();
  float dz1 = caloPos.Z() - c->vtx.Z();
  float r1 = sqrt(dx1*dx1+dy1*dy1);
  float tantheta = r1/dz1;

  float dx2 = c->vtx.X()-basePos.X();
  float dy2 = c->vtx.Y()-basePos.Y();
  float r2  = sqrt(dx2*dx2+dy2*dy2);
  float dz2 = r2/tantheta;

  return caloPos.Z() - dz1 - dz2;
}

int HggVertexing::getNConv(VecbosPho* p1, VecbosPho* p2){
  int nconv=0;
  if(p1->conversion.index != -1) nconv++;
  if(p2->conversion.index != -1) nconv++;
  return nconv;
}
