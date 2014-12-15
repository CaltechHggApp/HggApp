#!/usr/bin/env python

import os
import sys
import ROOT as rt

m2_pts = range(125,501,25)
m1_pts = range(0,376,25)

combine_dir = "for_combine/"

def make1DLimit():
    error = rt.TGraphAsymmErrors(len(m2_pts))
    graph = rt.TGraph(len(m2_pts))

    for i in range(len(m2_pts)):
        m = m2_pts[i]
        limit_file = rt.TFile(combine_dir+"higgsCombineChiWH_0_%d.Asymptotic.mH120.root" % (m,))
        limit_tree = limit_file.Get("limit")
        limit_form = rt.TTreeFormula("get_limit","limit",limit_tree)

        limit_tree.GetEntry(1)
        down = limit_form.EvalInstance()
        limit_tree.GetEntry(2)
        exp = limit_form.EvalInstance()
        limit_tree.GetEntry(3)
        up = limit_form.EvalInstance()
        limit_tree.GetEntry(5)
        obs = limit_form.EvalInstance()

        graph.SetPoint(i,m,exp)
        error.SetPoint(i,m,exp)
        error.SetPointError(i,0,0,exp-down,up-exp);

    c = rt.TCanvas()
    c.SetLogy()
    error.SetFillColor(3)
    error.Draw("AP4")
    graph.Draw("PC")
    c.SaveAs(combine_dir+"expected_exclusion_1D.png")

make1DLimit()
