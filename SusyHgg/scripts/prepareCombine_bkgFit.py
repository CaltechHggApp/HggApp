#!/usr/bin/env python

import ROOT as rt
import sys

def prepare(inputFile, SMSFile,SMSTag, outputFile):
    inf = rt.TFile(inputFile)
    w = inf.Get("susy_hgg_workspace")
    cats = w.obj("evtcat")
    catNames = []
    i=0
    while cats.setIndex(i)==False:
        catNames.append(cats.getLabel())
        i+=1

    xvar = w.obj("signal_region")
    n = xvar.numBins()

    sideband = rt.TH1D("sideband","",n * len(catNames),-0.5,n*len(catNames)+2)
    smhiggs = rt.TH1D("smhiggs","",n * len(catNames),-0.5,n*len(catNames)+2)
    signal = rt.TH1D("data_obs","",n * len(catNames),-0.5,n*len(catNames)+2)

    smsFile = rt.TFile(SMSFile)
    smsw    = smsFile.Get("susy_hgg_workspace")
    sms_signal = rt.TH1D("signal","",n * len(catNames),-0.5,n*len(catNames)+2)

    outf = rt.TFile(outputFile+".root","RECREATE")
    N_total_Higgs=0
    N_total_bkg  =0

    for iCat in range(len(catNames)):
        cat = catNames[iCat]
        #N_total_Higgs+=inf.Get("SMTot_%s_Signal_SigRegions").Integral()
        #N_total_Higgs+=inf.Get("SMTot_%s_Signal_SigRegions").Integral()
        #N_total_bkg+=   w.obj( "data_"+cat+"_Sideband_SigRegions_norm_hist").sumEntries()
        #N_total_Higgs+= w.obj( "SMTot_"+cat+"_Signal_SigRegions_norm_hist").sumEntries()

    print N_total_bkg
    print N_total_Higgs
    
    for iCat in range(len(catNames)):
        cat = catNames[iCat]
        smh  = inf.Get("SMTot_%s_Signal_SigRegions" % cat)
        sig  = inf.Get("data_%s_Signal_sidebandSub_SigRegions"%cat)
        side = inf.Get("data_%s_Sideband_SigRegions"%cat)
        sms  = smsFile.Get(SMSTag+"_"+cat+"_Signal_SigRegions")

        xbins = smh.GetNbinsX()

        for i in range(n):
            sideband.SetBinContent( n*iCat+i+1, side.GetBinContent( (i+1)%xbins, (i+1)/xbins ) )
            smhiggs.SetBinContent( n*iCat+i+1, smh.GetBinContent( (i+1)%xbins, (i+1)/xbins ) )
            signal.SetBinContent( n*iCat+i+1, sig.GetBinContent( (i+1)%xbins, (i+1)/xbins ) )
            sms_signal.SetBinContent( n*iCat+i+1, sms.GetBinContent( (i+1)%xbins, (i+1)/xbins ) )

    
    print "Total Events: "
    print "Observed:  "+str(signal.Integral())
    print "Expected:  "+str(smhiggs.Integral())
    print "SMS:       "+str(sms_signal.Integral())

    writeScript(outputFile,signal.Integral(),["signal","smhiggs"],[sms_signal.Integral(),smhiggs.Integral()])

    outf.cd()
    sideband.Write()
    smhiggs.Write()
    signal.Write()
    sms_signal.Write()
    outf.Close()


def writeScript(outputName,obs,names,rates):
    f = open(outputName+".txt",'w')

    f.write("imax   1\n")
    f.write("jmax   *\n")
    f.write("kmax   *\n")

    f.write("-------------------------------------\n")
    f.write("bin 1\n")
    f.write("observation   "+str(obs)+"\n")
    f.write("-------------------------------------\n")
    f.write("shapes * * "+str(outputName)+".root $PROCESS $PROCESS_$SYSTEMATIC\n")
    f.write("-------------------------------------\n")

    line1 ="bin    "
    line2 ="process"
    line3 ="process"
    line4 ="rate   "
    print zip(names,rates)
    for i,(n,r) in enumerate(zip(names,rates)):
        line1 += "%40d" % 1
        line2 += "%40s" % n
        line3 += "%40s" % i
        line4 += "%40s" % r
    f.write(line1+"\n")
    f.write(line2+"\n")
    f.write(line3+"\n")
    f.write(line4+"\n")
    f.write("-------------------------------------\n")

    f.close()

def usage(name):
    print "USAGE:  "+name+" <input data file> <input SMS> <SMS tag> <output tag>"
    sys.exit(0)
 
if __name__=="__main__":
    if len(sys.argv) <5:
        usage(sys.argv[0])

    input = sys.argv[1]
    smsfile = sys.argv[2]
    smstag  = sys.argv[3]
    output = sys.argv[4]

    prepare(input,smsfile,smstag,output)
