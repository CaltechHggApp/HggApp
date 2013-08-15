#!/usr/bin/env python

from ROOT import gROOT,TTree,TH2F,TFile
from samples import *

def dump(ifn,ofn,hn):
    fi = TFile(ifn)
    t  = fi.Get("output")
    print type(t)
    d = TH2F(hn,"",600,-3,3,100,0,0.1)
    t.Project(hn,"se:etaSC")    
    fo = TFile(ofn,"RECREATE")
    fo.cd()
    d.Write()
    fo.Close()
    fi.Close()

gROOT.Reset()

for h in samples:
    dump("output/"+h+".root",h+"_hist_dump.root",h)
