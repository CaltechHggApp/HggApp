#!/usr/bin/env python

import ROOT as rt
import sys

#type = "sidebandSub_"
type=""

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
        print cat
        N_total_bkg+=   w.obj( "data_"+cat+"_Sideband_SigRegions_norm_hist").sumEntries()
        N_total_Higgs+= w.obj( "SMTot_"+cat+"_Signal_SigRegions_norm_hist").sumEntries()

    print N_total_bkg
    print N_total_Higgs
    
    for iCat in range(len(catNames)):
        cat = catNames[iCat]
        side = w.obj( "data_"+cat+"_Sideband_SigRegions_hist")
        smh  = w.obj( "SMTot_"+cat+"_Signal_SigRegions_hist")
        sig  = w.obj( "data_"+cat+"_Signal_%sSigRegions_hist"%type)
        sms  = smsw.obj(SMSTag+"_"+cat+"_Signal_SigRegions_hist")

    
        th_side = side.createHistogram("sideband_"+cat,xvar)
        th_smh  = smh.createHistogram("smhiggs_"+cat,xvar)
        th_sig  = sig.createHistogram("signal_"+cat,xvar)
        th_sms   = sms.createHistogram("sms_"+cat,xvar)

        sig_region_integral = w.var("signal_fit_int_"+cat)
        side_band_sum = th_side.Integral()

        for i in range(n):
            sideband.SetBinContent( n*iCat+i+1, th_side.GetBinContent( i+1 ) *sig_region_integral.getVal()/side_band_sum )
            smhiggs.SetBinContent( n*iCat+i+1, th_smh.GetBinContent( i+1 ) )
            signal.SetBinContent( n*iCat+i+1, th_sig.GetBinContent( i+1 ) )
            sms_signal.SetBinContent( n*iCat+i+1, th_sms.GetBinContent( i+1 ) )

            

    '''
    N_Higgs = w.var("N_FITTED_%sSMTot"%type)
    N_bkg   = w.var("N_FITTED_%sSideband"%type)

    print N_Higgs.getVal()
    print N_bkg.getVal()
    
    sideband.Scale( N_bkg.getVal()/N_total_bkg )
    smhiggs.Scale(  N_Higgs.getVal()/N_total_Higgs )
    '''
    background = sideband.Clone("background")
    #background.Add(smhiggs)

    #background = smhiggs.Clone("background")

    print "Total Events: "
    print "Observed:  "+str(signal.Integral())
    print "Expected:  "+str(background.Integral())
    print "SMS:       "+str(sms_signal.Integral())

    writeScript(outputFile,signal.Integral(),["signal","background","smhiggs"],[sms_signal.Integral(),background.Integral(),smhiggs.Integral()])

    outf.cd()
    sideband.Write()
    smhiggs.Write()
    signal.Write()
    background.Write()
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

    if smstag == "ALL_W_0":
        for mcharg in range(125,501,25):
            # scan mchargino with mlsp=0
            print mcharg
            prepare(input,smsfile,"sms_ChiWH_0_%d"%mcharg,output+"/sms_ChiWH_0_%d"%mcharg)
    else:
        prepare(input,smsfile,smstag,output)
