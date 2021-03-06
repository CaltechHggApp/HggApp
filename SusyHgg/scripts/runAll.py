#!/usr/bin/env python

import os
import sys

m2_pts = range(125,501,25)
m1_pts = range(0,376,25)

combine_dir = "for_combine/"

def runAll():
    if not os.path.exists(combine_dir):
        os.mkdir(combine_dir)
    
    for m1 in m1_pts:
        for m2 in m2_pts:
            if m2-m1 < 125: continue
            ## build the script for combine

            #cmd = "python scripts/prepareCombine.py fits/data_ABCD_noR9_iso_NewRegions_smhiggs.root fits/sms_19.78_noR9_iso_NewRegions_ChiWH.root sms_ChiWH_%d_%d %s/sms_ChiWH_%d_%d" % ( m1,m2,combine_dir,m1,m2) 
            cmd = "python scripts/prepareCombine.py fits/signalInj_ABCD_ChiWH_0_175_1pb__noR9_iso_NewRegions_SigEoE_PtGG110.root fits/sms_19.78_noR9_iso_NewRegions_SigEoE_PtGG110_ChiWH.root sms_ChiWH_%d_%d %s/sms_ChiWH_%d_%d" % ( m1,m2,combine_dir,m1,m2) 
            print cmd
            os.system(cmd)  

            ## run combine
            cmd = "combine -M Asymptotic %s/sms_ChiWH_%d_%d.txt -n ChiWH_%d_%d" % (combine_dir,m1,m2,m1,m2)
            print cmd
            os.system(cmd)
            ##

            #move the output into the combine dir
            cmd = "mv higgsCombineChiWH_%d_%d.Asymptotic.mH120.root %s" %(m1,m2,combine_dir)
            print cmd
            os.system(cmd)

if __name__=='__main__':
    runAll()
