#ifndef selectionMaps_C
#define selectionMaps_C


#include "TF1.h"
#include "TH2F.h"

TH2F* getSelectionMap0(){
  float etaBins[5]={0.,1.,1.48,2.,3.};
  float ptBins[5]={0.,40.,60.,160.,2000.};
  TH2F* map = new TH2F("map","",4,etaBins,4,ptBins);
  map->SetBinContent(1,1,0.0076); 
  map->SetBinContent(1,2,0.0076);
  map->SetBinContent(1,3,0.0076);
  map->SetBinContent(1,4,0.0076);

  map->SetBinContent(2,1,0.0076);
  map->SetBinContent(2,2,0.0076);
  map->SetBinContent(2,3,0.0076);
  map->SetBinContent(2,4,0.0076);

  map->SetBinContent(3,1,0.0118);
  map->SetBinContent(3,2,0.0118);
  map->SetBinContent(3,3,0.0118);
  map->SetBinContent(3,4,0.0118);

  map->SetBinContent(4,1,0.0118);
  map->SetBinContent(4,2,0.0118);
  map->SetBinContent(4,3,0.0118);
  map->SetBinContent(4,4,0.0118);

  return map;
}

TH2F* getSelectionMap1(){

  float etaBins[5]={0.,1.,1.48,2.,3.};
  float ptBins[5]={0.,40.,60.,160.,2000.};
  TH2F* map = new TH2F("map","",4,etaBins,4,ptBins);
  map->SetBinContent(1,1,0.0086); 
  map->SetBinContent(1,2,0.0081);
  map->SetBinContent(1,3,0.0076);
  map->SetBinContent(1,4,0.0076);

  map->SetBinContent(2,1,0.0086);
  map->SetBinContent(2,2,0.0081);
  map->SetBinContent(2,3,0.0076);
  map->SetBinContent(2,4,0.0076);

  map->SetBinContent(3,1,0.0138);
  map->SetBinContent(3,2,0.0128);
  map->SetBinContent(3,3,0.0118);
  map->SetBinContent(3,4,0.0118);

  map->SetBinContent(4,1,0.0138);
  map->SetBinContent(4,2,0.0128);
  map->SetBinContent(4,3,0.0118);
  map->SetBinContent(4,4,0.0118);

  return map;
  

}


TH2F* getSelectionMap2(){
  TF1 EBfunc("EBFunc","[0]+expo(1)",0,2000);
  EBfunc.SetParameter(0,0.0065);
  EBfunc.SetParameter(1,-5.23);
  EBfunc.SetParameter(2,-0.0267);

  TF1 EEfunc("EEFunc","[0]+expo(1)",0,2000);
  EEfunc.SetParameter(0,0.0117);
  EEfunc.SetParameter(1,-3.71);
  EEfunc.SetParameter(2,-0.0550);
  
  float etaBins[5]={0.,1.,1.48,2.,3.};
  float ptBins[2001];
  for(float i=0;i<2001;i++) ptBins[int(i)]=i;
  TH2F* map = new TH2F("map","",4,etaBins,2000,ptBins);
  for(int i=1;i<2001;i++){
    map->SetBinContent(1,i,EBfunc.Eval(i));
    map->SetBinContent(2,i,EBfunc.Eval(i));
    map->SetBinContent(3,i,EEfunc.Eval(i));
    map->SetBinContent(4,i,EEfunc.Eval(i));
  }
  return map;
}

TH2F* getSelectionMap3(){
  TF1 EBlfunc("EBlFunc","[0]+expo(1)",0,2000);
  EBlfunc.SetParameter(0,0.0065);
  EBlfunc.SetParameter(1,-5.23);
  EBlfunc.SetParameter(2,-0.0267);

  TF1 EBhfunc("EBhFunc","[0]+expo(1)",0,2000);
  EBhfunc.SetParameter(0,0.0071);
  EBhfunc.SetParameter(1,-6.28);
  EBhfunc.SetParameter(2,-0.0099);

  TF1 EElfunc("EElFunc","[0]+expo(1)",0,2000);
  EElfunc.SetParameter(0,0.0117);
  EElfunc.SetParameter(1,-3.71);
  EElfunc.SetParameter(2,-0.0550);

  TF1 EEhfunc("EEhFunc","[0]+expo(1)",0,2000);
  EEhfunc.SetParameter(0,0.0117);
  EEhfunc.SetParameter(1,-4.27);
  EEhfunc.SetParameter(2,-0.0541);
  
  float etaBins[5]={0.,1.,1.48,2.,3.};
  float ptBins[8001];
  for(float i=0;i<8001;i++) ptBins[int(i)]=i/4;
  TH2F* map = new TH2F("map","",4,etaBins,8000,ptBins);
  for(int i=1;i<8001;i++){
    map->SetBinContent(1,i,EBlfunc.Eval(i/4-0.125));
    map->SetBinContent(2,i,EBhfunc.Eval(i/4-0.125));
    map->SetBinContent(3,i,EElfunc.Eval(i/4-0.125));
    map->SetBinContent(4,i,EEhfunc.Eval(i/4-0.125));
  }
  return map;
}

TH2F* getSelectionMap4(bool isData){
  TF1 EBlfunc("EBlFunc","[0]+expo(1)",0,2000);
  TF1 EBhfunc("EBhFunc","[0]+expo(1)",0,2000);
  TF1 EElfunc("EElFunc","[0]+expo(1)",0,2000);
  TF1 EEhfunc("EEhFunc","[0]+expo(1)",0,2000);

  if(isData){
    EBlfunc.SetParameter(0,0.00644);
    EBlfunc.SetParameter(1,-5.346);
    EBlfunc.SetParameter(2,-0.02506);
    
    EBhfunc.SetParameter(0,0.007766);
    EBhfunc.SetParameter(1,-5.47);
    EBhfunc.SetParameter(2,-0.0433);

    EElfunc.SetParameter(0,0.01176);
    EElfunc.SetParameter(1,-4.313);
    EElfunc.SetParameter(2,-0.04477);

    EEhfunc.SetParameter(0,0.01151);
    EEhfunc.SetParameter(1,-4.139);
    EEhfunc.SetParameter(2,-0.0511);
  }else{
    EBlfunc.SetParameter(0,0.00703);
    EBlfunc.SetParameter(1,-6.585);
    EBlfunc.SetParameter(2,-0.01328);
    
    EBhfunc.SetParameter(0,0.0075);
    EBhfunc.SetParameter(1,-6.625);
    EBhfunc.SetParameter(2,-0.0155);

    EElfunc.SetParameter(0,0.01184);
    EElfunc.SetParameter(1,-4.984);
    EElfunc.SetParameter(2,-0.0344);

    EEhfunc.SetParameter(0,0.0111);
    EEhfunc.SetParameter(1,-4.492);
    EEhfunc.SetParameter(2,-0.0407);
  }
  
  float etaBins[5]={0.,1.,1.48,2.,3.};
  float ptBins[8001];
  for(float i=0;i<8001;i++) ptBins[int(i)]=i/4;
  TH2F* map = new TH2F("map","",4,etaBins,8000,ptBins);
  for(int i=1;i<8001;i++){
    map->SetBinContent(1,i,EBlfunc.Eval(i/4+0.125));
    map->SetBinContent(2,i,EBhfunc.Eval(i/4+0.125));
    map->SetBinContent(3,i,EElfunc.Eval(i/4+0.125));
    map->SetBinContent(4,i,EEhfunc.Eval(i/4+0.125));
  }
  return map;
}

TH2F* getSelectionMap5(bool isData){
  TF1 EBlfunc("EBlFunc","[0]+expo(1)",0,2000);
  TF1 EBhfunc("EBhFunc","[0]+expo(1)",0,2000);
  TF1 EElfunc("EElFunc","[0]+expo(1)",0,2000);
  TF1 EEhfunc("EEhFunc","[0]+expo(1)",0,2000);

  EBlfunc.SetParameter(0,0.00644);
  EBlfunc.SetParameter(1,-5.346);
  EBlfunc.SetParameter(2,-0.02506);
  
  EBhfunc.SetParameter(0,0.007766);
  EBhfunc.SetParameter(1,-5.47);
  EBhfunc.SetParameter(2,-0.0433);
  
  EElfunc.SetParameter(0,0.01176);
  EElfunc.SetParameter(1,-4.313);
  EElfunc.SetParameter(2,-0.04477);
  
  EEhfunc.SetParameter(0,0.01151);
  EEhfunc.SetParameter(1,-4.139);
  EEhfunc.SetParameter(2,-0.0511);
  if(!isData){
    EBlfunc.SetParameter(0,0.00794);
    
    EBhfunc.SetParameter(0,0.009266);

    EElfunc.SetParameter(0,0.01176);

    EEhfunc.SetParameter(0,0.01151);
  }
  
  float etaBins[5]={0.,1.,1.48,2.,3.};
  float ptBins[8001];
  for(float i=0;i<8001;i++) ptBins[int(i)]=i/4;
  TH2F* map = new TH2F("map","",4,etaBins,8000,ptBins);
  for(int i=1;i<8001;i++){
    map->SetBinContent(1,i,EBlfunc.Eval(i/4+0.125));
    map->SetBinContent(2,i,EBhfunc.Eval(i/4+0.125));
    map->SetBinContent(3,i,EElfunc.Eval(i/4+0.125));
    map->SetBinContent(4,i,EEhfunc.Eval(i/4+0.125));
  }
  return map;
}


TH2F* getSelectionMap6(bool isData){
  TF1 EBlfunc("EBlFunc","[0]+expo(1)",0,2000);
  TF1 EBhfunc("EBhFunc","[0]+expo(1)",0,2000);
  TF1 EElfunc("EElFunc","[0]+expo(1)",0,2000);
  TF1 EEhfunc("EEhFunc","[0]+expo(1)",0,2000);

  EBlfunc.SetParameter(0,0.00644+0.0014);
  EBlfunc.SetParameter(1,-5.346);
  EBlfunc.SetParameter(2,-0.02506);
  
  EBhfunc.SetParameter(0,0.007766+0.0014);
  EBhfunc.SetParameter(1,-5.47);
  EBhfunc.SetParameter(2,-0.0433);
  
  EElfunc.SetParameter(0,0.01176+0.004);
  EElfunc.SetParameter(1,-4.313);
  EElfunc.SetParameter(2,-0.04477);
  
  EEhfunc.SetParameter(0,0.01151+0.004);
  EEhfunc.SetParameter(1,-4.139);
  EEhfunc.SetParameter(2,-0.0511);
  
  float etaBins[5]={0.,1.,1.48,2.,3.};
  float ptBins[8001];
  for(float i=0;i<8001;i++) ptBins[int(i)]=i/4;
  TH2F* map = new TH2F("map","",4,etaBins,8000,ptBins);
  for(int i=1;i<8001;i++){
    map->SetBinContent(1,i,EBlfunc.Eval(i/4+0.125));
    map->SetBinContent(2,i,EBhfunc.Eval(i/4+0.125));
    map->SetBinContent(3,i,EElfunc.Eval(i/4+0.125));
    map->SetBinContent(4,i,EEhfunc.Eval(i/4+0.125));
  }
  return map;
}

TH2F* getSelectionMap7(bool isData){
  TF1 EBlfunc("EBlFunc","[0]+expo(1)",0,2000);
  TF1 EBhfunc("EBhFunc","[0]+expo(1)",0,2000);
  TF1 EElfunc("EElFunc","[0]+expo(1)",0,2000);
  TF1 EEhfunc("EEhFunc","[0]+expo(1)",0,2000);

  EBlfunc.SetParameter(0,0.00644+0.001);
  EBlfunc.SetParameter(1,-5.346);
  EBlfunc.SetParameter(2,-0.02506);
  
  EBhfunc.SetParameter(0,0.007766+0.001);
  EBhfunc.SetParameter(1,-5.47);
  EBhfunc.SetParameter(2,-0.0433);
  
  EElfunc.SetParameter(0,0.01176+0.012);
  EElfunc.SetParameter(1,-4.313);
  EElfunc.SetParameter(2,-0.04477);
  
  EEhfunc.SetParameter(0,0.01151+0.012);
  EEhfunc.SetParameter(1,-4.139);
  EEhfunc.SetParameter(2,-0.0511);
  
  float etaBins[5]={0.,1.,1.48,2.,3.};
  float ptBins[8001];
  for(float i=0;i<8001;i++) ptBins[int(i)]=i/4;
  TH2F* map = new TH2F("map","",4,etaBins,8000,ptBins);
  for(int i=1;i<8001;i++){
    map->SetBinContent(1,i,EBlfunc.Eval(i/4+0.125));
    map->SetBinContent(2,i,EBhfunc.Eval(i/4+0.125));
    map->SetBinContent(3,i,EElfunc.Eval(i/4+0.125));
    map->SetBinContent(4,i,EEhfunc.Eval(i/4+0.125));
  }
  return map;
}

std::vector<TH2F*> getSelectionMap8(){
  float etaBins[5]={0.,1.,1.48,2.,3.};
  float ptBins[2]={-100.,100000.};
  std::vector<TH2F*> ret;  
  ret.push_back(new TH2F("map1","",4,etaBins,1,ptBins));
  ret.back()->SetBinContent(1,1,0.020); 
  ret.back()->SetBinContent(2,1,0.025);
  ret.back()->SetBinContent(3,1,0.020);
  ret.back()->SetBinContent(4,1,0.025);

  ret.push_back(new TH2F("map2","",4,etaBins,1,ptBins));
  ret.back()->SetBinContent(1,1,0.030); 
  ret.back()->SetBinContent(2,1,0.035);
  ret.back()->SetBinContent(3,1,0.030);
  ret.back()->SetBinContent(4,1,0.035);

  return ret;
}

std::vector<TH2F*> getSelectionMap9(){
  float etaBins[5]={0.,1.,1.48,2.,3.};
  float ptBins[2]={-100.,100000.};
  std::vector<TH2F*> ret;  
  ret.push_back(new TH2F("map0","",4,etaBins,1,ptBins));
  ret.back()->SetBinContent(1,1,0.010); 
  ret.back()->SetBinContent(2,1,0.010);
  ret.back()->SetBinContent(3,1,0.020);
  ret.back()->SetBinContent(4,1,0.020);

  ret.push_back(new TH2F("map1","",4,etaBins,1,ptBins));
  ret.back()->SetBinContent(1,1,0.015); 
  ret.back()->SetBinContent(2,1,0.015);
  ret.back()->SetBinContent(3,1,0.025);
  ret.back()->SetBinContent(4,1,0.025);

  ret.push_back(new TH2F("map2","",4,etaBins,1,ptBins));
  ret.back()->SetBinContent(1,1,0.020); 
  ret.back()->SetBinContent(2,1,0.020);
  ret.back()->SetBinContent(3,1,0.030);
  ret.back()->SetBinContent(4,1,0.030);

  ret.push_back(new TH2F("map3","",4,etaBins,1,ptBins));
  ret.back()->SetBinContent(1,1,0.030); 
  ret.back()->SetBinContent(2,1,0.030);
  ret.back()->SetBinContent(3,1,0.040);
  ret.back()->SetBinContent(4,1,0.040);

  ret.push_back(new TH2F("map4","",4,etaBins,1,ptBins));
  ret.back()->SetBinContent(1,1,0.040); 
  ret.back()->SetBinContent(2,1,0.040);
  ret.back()->SetBinContent(3,1,0.050);
  ret.back()->SetBinContent(4,1,0.050);

  return ret;
}

std::vector<TH2F*> getSelectionMap10(){ //optimized on ROC curves with CiC Isolation
  float etaBins[5]={0.,1.,1.48,2.,3.};
  float ptBins[2]={-100.,100000.};
  std::vector<TH2F*> ret;  
  ret.push_back(new TH2F("map1","",4,etaBins,1,ptBins));
  ret.back()->SetBinContent(1,1,0.018); // best photon Eff + fake Rejection
  ret.back()->SetBinContent(2,1,0.03); 
  ret.back()->SetBinContent(3,1,0.027);
  ret.back()->SetBinContent(4,1,0.021);

  ret.push_back(new TH2F("map2","",4,etaBins,1,ptBins));
  ret.back()->SetBinContent(1,1,0.0285); // 90% real photon Efficiency
  ret.back()->SetBinContent(2,1,0.0435);
  ret.back()->SetBinContent(3,1,0.0465);
  ret.back()->SetBinContent(4,1,0.0285);

  return ret;
}

std::vector<TH2F*> getSelectionMap11(){ //optimized on ROC curves with CiC Isolation
  float etaBins[5]={0.,1.,1.48,2.,3.};  
  float ptBins[2]={-100.,100000.};
  std::vector<TH2F*> ret;  
  ret.push_back(new TH2F("map1","",4,etaBins,1,ptBins));
  ret.back()->SetBinContent(1,1,0.0105); // 90% fake rejection
  ret.back()->SetBinContent(2,1,0.015); 
  ret.back()->SetBinContent(3,1,0.0195);
  ret.back()->SetBinContent(4,1,0.018);

  ret.push_back(new TH2F("map2","",4,etaBins,1,ptBins));
  ret.back()->SetBinContent(1,1,0.0285); // 90% real photon Efficiency
  ret.back()->SetBinContent(2,1,0.0435);
  ret.back()->SetBinContent(3,1,0.0465);
  ret.back()->SetBinContent(4,1,0.0285);

  return ret;
}

std::vector<TH2F*> getSelectionMap12(){ //optimized on ROC curves with CiC Isolation in pT bins
  float etaBins[5]={0.,1.,1.48,2.,3.};  
  float ptBins[6]={0.,30.,40.,60.,100.,100000.};
  std::vector<TH2F*> ret;  
  ret.push_back(new TH2F("map1","",4,etaBins,5,ptBins));
  //0-30
  ret.back()->SetBinContent(1,1,0.012); // 90% fake rejection
  ret.back()->SetBinContent(2,1,0.018); 
  ret.back()->SetBinContent(3,1,0.021);
  ret.back()->SetBinContent(4,1,0.021);
  //30-40
  ret.back()->SetBinContent(1,2,0.0105);
  ret.back()->SetBinContent(2,2,0.018); 
  ret.back()->SetBinContent(3,2,0.0195);
  ret.back()->SetBinContent(4,2,0.0195);
  //40-60
  ret.back()->SetBinContent(1,3,0.0105);
  ret.back()->SetBinContent(2,3,0.0165); 
  ret.back()->SetBinContent(3,3,0.0195);
  ret.back()->SetBinContent(4,3,0.018);
  //60-100
  ret.back()->SetBinContent(1,4,0.0105);
  ret.back()->SetBinContent(2,4,0.015); 
  ret.back()->SetBinContent(3,4,0.018);
  ret.back()->SetBinContent(4,4,0.0165);
  //100-
  ret.back()->SetBinContent(1,5,0.009);
  ret.back()->SetBinContent(2,5,0.012); 
  ret.back()->SetBinContent(3,5,0.0165);
  ret.back()->SetBinContent(4,5,0.015);

  ret.push_back(new TH2F("map2","",4,etaBins,1,ptBins));
  //0-30
  ret.back()->SetBinContent(1,1,0.03); // 90% real photon Efficiency
  ret.back()->SetBinContent(2,1,0.054);
  ret.back()->SetBinContent(3,1,0.0675);
  ret.back()->SetBinContent(4,1,0.039);
  //30-40
  ret.back()->SetBinContent(1,2,0.033);
  ret.back()->SetBinContent(2,2,0.0495);
  ret.back()->SetBinContent(3,2,0.055);
  ret.back()->SetBinContent(4,2,0.033);
  //40-60
  ret.back()->SetBinContent(1,3,0.03);
  ret.back()->SetBinContent(2,3,0.0435);
  ret.back()->SetBinContent(3,3,0.0435);
  ret.back()->SetBinContent(4,3,0.027);
  //60-100
  ret.back()->SetBinContent(1,4,0.0255);
  ret.back()->SetBinContent(2,4,0.0345);
  ret.back()->SetBinContent(3,4,0.0315);
  ret.back()->SetBinContent(4,4,0.0225);
  //100-
  ret.back()->SetBinContent(1,5,0.018);
  ret.back()->SetBinContent(2,5,0.0225);
  ret.back()->SetBinContent(3,5,0.0225);
  ret.back()->SetBinContent(4,5,0.021);

  return ret;
}

std::vector<TH2F*> getSelectionMap13(){ //optimized on ROC curves with CiC Isolation in pT bins
  float etaBins[5]={0.,0.8,1.48,2.,3.};  
  float ptBins[6]={0.,30.,40.,60.,100.,100000.};
  std::vector<TH2F*> ret;  
  ret.push_back(new TH2F("map1","",4,etaBins,5,ptBins));
  //0-30
  ret.back()->SetBinContent(1,1,0.012); // 90% fake rejection
  ret.back()->SetBinContent(2,1,0.018); 
  ret.back()->SetBinContent(3,1,0.021);
  ret.back()->SetBinContent(4,1,0.021);
  //30-40
  ret.back()->SetBinContent(1,2,0.0105);
  ret.back()->SetBinContent(2,2,0.018); 
  ret.back()->SetBinContent(3,2,0.0195);
  ret.back()->SetBinContent(4,2,0.0195);
  //40-60
  ret.back()->SetBinContent(1,3,0.0105);
  ret.back()->SetBinContent(2,3,0.0165); 
  ret.back()->SetBinContent(3,3,0.0195);
  ret.back()->SetBinContent(4,3,0.018);
  //60-100
  ret.back()->SetBinContent(1,4,0.0105);
  ret.back()->SetBinContent(2,4,0.015); 
  ret.back()->SetBinContent(3,4,0.018);
  ret.back()->SetBinContent(4,4,0.0165);
  //100-
  ret.back()->SetBinContent(1,5,0.009);
  ret.back()->SetBinContent(2,5,0.012); 
  ret.back()->SetBinContent(3,5,0.0165);
  ret.back()->SetBinContent(4,5,0.015);

  ret.push_back(new TH2F("map2","",4,etaBins,1,ptBins));
  //0-30
  ret.back()->SetBinContent(1,1,0.03); // 90% real photon Efficiency
  ret.back()->SetBinContent(2,1,0.054);
  ret.back()->SetBinContent(3,1,0.0675);
  ret.back()->SetBinContent(4,1,0.039);
  //30-40
  ret.back()->SetBinContent(1,2,0.033);
  ret.back()->SetBinContent(2,2,0.0495);
  ret.back()->SetBinContent(3,2,0.055);
  ret.back()->SetBinContent(4,2,0.033);
  //40-60
  ret.back()->SetBinContent(1,3,0.03);
  ret.back()->SetBinContent(2,3,0.0435);
  ret.back()->SetBinContent(3,3,0.0435);
  ret.back()->SetBinContent(4,3,0.027);
  //60-100
  ret.back()->SetBinContent(1,4,0.0255);
  ret.back()->SetBinContent(2,4,0.0345);
  ret.back()->SetBinContent(3,4,0.0315);
  ret.back()->SetBinContent(4,4,0.0225);
  //100-
  ret.back()->SetBinContent(1,5,0.018);
  ret.back()->SetBinContent(2,5,0.0225);
  ret.back()->SetBinContent(3,5,0.0225);
  ret.back()->SetBinContent(4,5,0.021);

  return ret;
}

std::vector<TH2F*> getSelectionMap14(){ //optimized on ROC curves with CiC Isolation in pT bins
  float etaBins[5]={0.,0.8,1.48,2.,3.};  
  float ptBins[6]={0.,30.,40.,60.,100.,100000.};
  std::vector<TH2F*> ret;  
  ret.push_back(new TH2F("map1","",4,etaBins,5,ptBins));
  //0-30
  ret.back()->SetBinContent(1,1,0.0086); // 90% fake rejection
  ret.back()->SetBinContent(2,1,0.0091); 
  ret.back()->SetBinContent(3,1,0.0140);
  ret.back()->SetBinContent(4,1,0.0130);
  //30-40
  ret.back()->SetBinContent(1,2,0.0086);
  ret.back()->SetBinContent(2,2,0.0091); 
  ret.back()->SetBinContent(3,2,0.0140);
  ret.back()->SetBinContent(4,2,0.0130);
  //40-60
  ret.back()->SetBinContent(1,3,0.0086);
  ret.back()->SetBinContent(2,3,0.0091); 
  ret.back()->SetBinContent(3,3,0.0140);
  ret.back()->SetBinContent(4,3,0.0130);
  //60-100
  ret.back()->SetBinContent(1,4,0.0086);
  ret.back()->SetBinContent(2,4,0.0091); 
  ret.back()->SetBinContent(3,4,0.0140);
  ret.back()->SetBinContent(4,4,0.0130);
  //100-
  ret.back()->SetBinContent(1,5,0.0086);
  ret.back()->SetBinContent(2,5,0.0091); 
  ret.back()->SetBinContent(3,5,0.0140);
  ret.back()->SetBinContent(4,5,0.0130);

  /*
  ret.push_back(new TH2F("map2","",4,etaBins,1,ptBins));
  //0-30
  ret.back()->SetBinContent(1,1,0.03); // 90% real photon Efficiency
  ret.back()->SetBinContent(2,1,0.054);
  ret.back()->SetBinContent(3,1,0.0675);
  ret.back()->SetBinContent(4,1,0.039);
  //30-40
  ret.back()->SetBinContent(1,2,0.033);
  ret.back()->SetBinContent(2,2,0.0495);
  ret.back()->SetBinContent(3,2,0.055);
  ret.back()->SetBinContent(4,2,0.033);
  //40-60
  ret.back()->SetBinContent(1,3,0.03);
  ret.back()->SetBinContent(2,3,0.0435);
  ret.back()->SetBinContent(3,3,0.0435);
  ret.back()->SetBinContent(4,3,0.027);
  //60-100
  ret.back()->SetBinContent(1,4,0.0255);
  ret.back()->SetBinContent(2,4,0.0345);
  ret.back()->SetBinContent(3,4,0.0315);
  ret.back()->SetBinContent(4,4,0.0225);
  //100-
  ret.back()->SetBinContent(1,5,0.018);
  ret.back()->SetBinContent(2,5,0.0225);
  ret.back()->SetBinContent(3,5,0.0225);
  ret.back()->SetBinContent(4,5,0.021);
  */
  return ret;
}


#endif
