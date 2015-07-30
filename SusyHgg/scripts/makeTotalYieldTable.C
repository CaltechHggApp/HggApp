#include "TFile.h"
#include "TString.h"
#include "TH2F.h"
#include "TMath.h"

#include <cstdio>
#include <iostream>

void makeTotalYieldTable(TString dir) {
  const int nDataSysts=3;
  TString dataSysts[nDataSysts] = {"fit","sfSyst","statistics"};

  const int nMCSysts=5;
  TString mcSysts[nMCSysts] = {"btag","jec","phoE","sigE","statistics"};


  const int nSMH=4;
  TString smh[nSMH] = {"ggH","vbfH","wzH","ttH"};

  float QCDScaleUP[nSMH] = {1.072,1.002,1.021,1.038};
  float QCDScaleDOWN[nSMH] = {0.928,0.998,0.982,0.915};

  float pdfUP[nSMH] = {1.075,1.026,1.024,0.925};
  float pdfDOWN[nSMH] = {0.935,0.973,0.977,1.081};

  float lumi=1.025;
  float trigger=1.05;
  //float trigger=1.01;

  const int nSMS=2;
  TString sms[nSMS]={"sms_ChiHH","sms_ChiWH"};
  const int nPts=2;
  TString pts[nPts]={"0_125","0_175"};

  float xsec[nSMS*nPts] = {3.248,1.027,4.12,1.32};
  float xsecErr[nSMS*nPts] = {0.1624,0.0514,0.255,0.0672};


  const int nCats=5;
  TString cats[nCats] = {"HighPt","Hbb","Zbb","HighRes","LowRes"};

  float dataYield[nCats],dataUp[nCats],dataDown[nCats];
  float smhYield[nCats],smhUp[nCats],smhDown[nCats];
  float smsYield[nSMS*nPts][nCats],smsUp[nSMS*nPts][nCats],smsDown[nSMS*nPts][nCats];
  for(int i=0; i<nCats; i++) {
    dataYield[i]=dataUp[i]=dataDown[i]=0;
    smhYield[i]=smhUp[i]=smhDown[i]=0;
    for(int j=0;j<nSMS*nPts;j++) {
      smsYield[j][i]=smsUp[j][i]=smsDown[j][i]=0;
    }
  }



  //data
  TFile dataFile(dir+"/data.root");
  for(int iCat=0; iCat<nCats; iCat++) {
    dataYield[iCat] = ((TH2F*)dataFile.Get("data_"+cats[iCat]+"_SidebandRegion"))->Integral();
    for(int iSyst=0; iSyst<nDataSysts; iSyst++) {
      float thisSystUP = ((TH2F*)dataFile.Get("data_"+cats[iCat]+"_SidebandRegion_"+dataSysts[iSyst]+"Up"))->Integral();
      float thisSystDOWN = ((TH2F*)dataFile.Get("data_"+cats[iCat]+"_SidebandRegion_"+dataSysts[iSyst]+"Down"))->Integral();
      dataUp[iCat]+= TMath::Power(thisSystUP-dataYield[iCat],2);
      dataDown[iCat]+= TMath::Power(thisSystDOWN-dataYield[iCat],2);      
      std::cout << cats[iCat] << " " <<  dataSysts[iSyst] << " " << dataYield[iCat] << " " << thisSystUP << " " << thisSystDOWN << std::endl;
    }

    dataUp[iCat] = sqrt(dataUp[iCat]);
    dataDown[iCat] = sqrt(dataDown[iCat]);
  }


  for(int iSMH=0; iSMH<nSMH; iSMH++) {
    TFile smhFile(dir+"/"+smh[iSMH]+".root");
    for(int iCat=0; iCat<nCats; iCat++) {
      float thisSMHYield=((TH2F*)smhFile.Get("data_"+cats[iCat]+"_SignalRegion"))->Integral();
      smhYield[iCat]+=thisSMHYield;
      for(int iSyst=0; iSyst<nMCSysts; iSyst++) {
	//std::cout << iSMH << " " << iCat << " " <<  iSyst << std::endl;
	float thisSystUP = ((TH2F*)smhFile.Get("data_"+cats[iCat]+"_SignalRegion_"+mcSysts[iSyst]+"_Up"))->Integral();
	float thisSystDOWN = ((TH2F*)smhFile.Get("data_"+cats[iCat]+"_SignalRegion_"+mcSysts[iSyst]+"_Down"))->Integral();
	smhUp[iCat]+= TMath::Power(thisSystUP-thisSMHYield,2);
	smhDown[iCat]+= TMath::Power(thisSystDOWN-thisSMHYield,2);      
      }//for(int iSyst=0
      smhUp[iCat]+=TMath::Power((QCDScaleUP[iSMH]-1)*thisSMHYield,2);
      smhDown[iCat]+=TMath::Power((QCDScaleDOWN[iSMH]-1)*thisSMHYield,2);

      smhUp[iCat]+=TMath::Power((pdfUP[iSMH]-1)*thisSMHYield,2);
      smhDown[iCat]+=TMath::Power((pdfDOWN[iSMH]-1)*thisSMHYield,2);      
    }//for(int iCat=0
  }//for(int iSMH=0

  for(int iSMS=0; iSMS<nSMS; iSMS++) {
    TFile smsFile(dir+"/"+sms[iSMS]+".root");
    for(int iPt=0; iPt<nPts; iPt++) {
      for(int iCat=0; iCat<nCats; iCat++) {
	float thisSMSYield=((TH2F*)smsFile.Get("data_"+pts[iPt]+"_"+cats[iCat]+"_SignalRegion"))->Integral()*xsec[iSMS*2+iPt];
	smsYield[iSMS*2+iPt][iCat]+=thisSMSYield;
	for(int iSyst=0; iSyst<nMCSysts; iSyst++) {
	  //std::cout << iSMH << " " << iCat << " " <<  iSyst << std::endl;
	  float thisSystUP = ((TH2F*)smsFile.Get("data_"+pts[iPt]+"_"+cats[iCat]+"_SignalRegion_"+mcSysts[iSyst]+"_Up"))->Integral()*xsec[iSMS*2+iPt];
	  float thisSystDOWN = ((TH2F*)smsFile.Get("data_"+pts[iPt]+"_"+cats[iCat]+"_SignalRegion_"+mcSysts[iSyst]+"_Down"))->Integral()*xsec[iSMS*2+iPt];
	  smsUp[iSMS*2+iPt][iCat]+= TMath::Power(thisSystUP-thisSMSYield,2);
	  smsDown[iSMS*2+iPt][iCat]+= TMath::Power(thisSystDOWN-thisSMSYield,2);      
	}//for(int iSyst=0
	std:: cout << TMath::Power(thisSMSYield*((xsecErr[iSMS*2+iPt]+xsec[iSMS*2+iPt])/xsec[iSMS*2+iPt]-1),2) << std::endl;
	smsUp[iSMS*2+iPt][iCat]+=TMath::Power(thisSMSYield*((xsecErr[iSMS*2+iPt]+xsec[iSMS*2+iPt])/xsec[iSMS*2+iPt]-1),2);
	smsDown[iSMS*2+iPt][iCat]+=TMath::Power(thisSMSYield*((xsecErr[iSMS*2+iPt]+xsec[iSMS*2+iPt])/xsec[iSMS*2+iPt]-1),2);
      }//for(int iCat=0
    }// for(int iPt=0
  }//for(int iSMS=0
  
  for(int iCat=0; iCat<nCats; iCat++) {
    smhUp[iCat] = sqrt(smhUp[iCat]);
    smhDown[iCat] = sqrt(smhDown[iCat]);
    printf(" %7s & $ % 4.1f \\pm {+% 3.2f} $ & $ % 4.1f^{+% 3.2f}_{-% 3.2f} $ & $ % 4.1f^{+% 3.2f}_{-% 3.2f} $ \\\\\n", 
	   cats[iCat].Data(), 
	   dataYield[iCat], dataUp[iCat],
	   smhYield[iCat], smhUp[iCat], smhDown[iCat],
	   dataYield[iCat]+smhYield[iCat], sqrt(dataUp[iCat]*dataUp[iCat]+smhUp[iCat]*smhUp[iCat]), 
	   sqrt(dataDown[iCat]*dataDown[iCat]+smhDown[iCat]*smhDown[iCat])
	   );
  }

  float tot=19780*2.28E-3;

  for(int iCat=0; iCat<nCats; iCat++) {
    for(int j=0; j < nSMS*nPts; j++) {
      smsUp[j][iCat] = sqrt(smsUp[j][iCat]);
      smsDown[j][iCat] = sqrt(smsDown[j][iCat]);
    }
    printf(" %7s & $ % 2.1f^{+% 2.1f}_{-% 2.1f} \\% $ & $ % 4.1f^{+% 3.2f}_{-% 3.2f} $ & $ % 2.1f^{+% 2.1f}_{-% 2.1f} \\% $ & $ % 4.1f^{+% 3.2f}_{-% 3.2f} $ & $ % 2.1f^{+% 2.1f}_{-% 2.1f} \\% $ & $ % 4.1f^{+% 3.2f}_{-% 3.2f} $ & $ % 2.1f^{+% 2.1f}_{-% 2.1f} \\% $ & $ % 4.1f^{+% 3.2f}_{-% 3.2f} $ \\\\\n",
	   cats[iCat].Data(), 
	   smsYield[0][iCat]/(tot*xsec[0])*100., smsUp[0][iCat]/(tot*xsec[0])*100., smsDown[0][iCat]/(tot*xsec[0])*100.,
	   smsYield[0][iCat], smsUp[0][iCat], smsDown[0][iCat],
	   smsYield[1][iCat]/(tot*xsec[1])*100., smsUp[1][iCat]/(tot*xsec[1])*100., smsDown[1][iCat]/(tot*xsec[1])*100.,
	   smsYield[1][iCat], smsUp[1][iCat], smsDown[1][iCat],
	   smsYield[2][iCat]/(tot*xsec[2])*100., smsUp[2][iCat]/(tot*xsec[2])*100., smsDown[2][iCat]/(tot*xsec[2])*100.,
	   smsYield[2][iCat], smsUp[2][iCat], smsDown[2][iCat],
	   smsYield[3][iCat]/(tot*xsec[3])*100., smsUp[3][iCat]/(tot*xsec[3])*100., smsDown[3][iCat]/(tot*xsec[3])*100.,
	   smsYield[3][iCat], smsUp[3][iCat], smsDown[3][iCat]
	   );
  }


}

