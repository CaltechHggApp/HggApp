#ifndef PassTrigger_hh
#define PassTrigger_hh

#include <string>
#include <vector>
#include <unordered_set>

#include "stdint.h"
#include "Vecbos.hh"

using namespace std;

struct evtInfo { 
  const static uint32_t NOLUMI=8191;
  uint32_t run,lumi,event; 
  evtInfo(uint32_t r, uint32_t l, uint32_t e){run=r; lumi=l; event=e;}
  evtInfo(int32_t r, int32_t l, int32_t e){run=static_cast<uint32_t>(r); lumi=static_cast<uint32_t>(l); event=static_cast<uint32_t>(e);}
  bool operator== (const evtInfo& e) const{
    return ((e.run==run) && (e.lumi==lumi || e.lumi==NOLUMI || lumi==NOLUMI) && (e.event==event));
  }

};

struct eiHash{
  size_t operator() (const evtInfo& ei) const {
    //return ei.run+size_t(ei.lumi)<<18+size_t(ei.event)<<31;
    return ei.run+size_t(ei.event)<<18;
  }
};

class PassTrigger : public Vecbos {
public:
  PassTrigger(TTree* vecbostree=0,TTree* analysisTree=0);

  void addTrigger(string t){ triggers.push_back(t); }

  void Loop(string outputFileName);


  void setRunLabel(string s) {runLabel=s;}
  void setLumiLabel(string s) {lumiLabel=s;}
  void setEvtLabel(string s) {evtLabel=s;}

  void setDoOutputStep(string s){inputText=s; outputStep=true;}

protected:
  void LoopOnVecbos(string outputFileName);
  void LoopOnOutput(string outputFileName);

  vector<string> triggers;

  bool outputStep=false;
  string inputText = "";

  string runLabel = "run";
  string lumiLabel= "lumi";
  string evtLabel="evt";

  unordered_set<evtInfo, eiHash> selectedEvents;
  unordered_set<evtInfo, eiHash> passTriggerEvents;
  TTree* analysisTree;

  //from the analysis tree
  int32_t run,evt,lumi;

};


#endif

