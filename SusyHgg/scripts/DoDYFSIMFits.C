{
  TFile *fastSIM_File = TFile::Open("/home/amott/raid4/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball__Summer12-START53_V19_FSIM_PU_S12-v1__WRONGSMEAR__SUSY.root");

  TTree *fastSIM_Tree = (TTree*)fastSIM_File->Get("SusyHggTree");

  RooRealVar mgg("mgg","",80,100);

  RooDataSet fastSIM_Data("fastSIM_Data","",fastSIM_Tree,mgg);

  RooRealVar fastSIM_m0("fastSIM_m0","",88,92);
  RooRealVar fastSIM_sig("fastSIM_sig","",2.495);
  RooRealVar fastSIM_alpha("fastSIM_alpha","",0,10);
  RooRealVar fastSIM_n("fastSIM_n","",0,10);

  RooRealVar fastSIM_smear_mean("fastSIM_smear_mean","",0);
  RooRealVar fastSIM_smear_sigma("fastSIM_smear_sigma","",0,5);

  RooCBShape fastSIM_cb("fastSIM_cb","",mgg,fastSIM_m0,fastSIM_sig,fastSIM_alpha,fastSIM_n);
  RooGaussian fastSIM_smear("fastSIM_smear","",mgg,fastSIM_smear_mean,fastSIM_smear_sigma);

  RooFFTConvPdf fastSIM_conv("fastSIM_conv","",mgg,fastSIM_cb,fastSIM_smear);

  RooRealVar fastSIM_bkg_alpha("fastSIM_bkg_alpha","",-3,0);

  RooExponential fastSIM_bkg("fastSIM_bkg","",mgg,fastSIM_bkg_alpha);

  RooRealVar fastSIM_N_sig("fastSIM_N_sig","",100,1e9);
  RooRealVar fastSIM_N_bkg("fastSIM_N_bkg","",10,500);

  RooAddPdf fastSIM_fitModel("fastSIM_fitModel","",RooArgList(fastSIM_conv,fastSIM_bkg),RooArgList(fastSIM_N_sig,fastSIM_N_bkg));

  fastSIM_fitModel.fitTo(fastSIM_Data,RooFit::Strategy(0));
  fastSIM_fitModel.fitTo(fastSIM_Data,RooFit::Strategy(2));

  RooPlot *fmgg = mgg.frame();
  fastSIM_Data.plotOn(fmgg);
  fastSIM_fitModel.plotOn(fmgg);

  TCanvas fastSIM_cv("fastSIM_cv","FastSIM");
  fmgg->Draw();


  TFile *fullSIM_File = TFile::Open("/home/amott/raid4/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball__SUSY.root");

  TTree *fullSIM_Tree = (TTree*)fullSIM_File->Get("SusyHggTree");

  RooRealVar mgg("mgg","",80,100);

  RooDataSet fullSIM_Data("fullSIM_Data","",fullSIM_Tree,mgg);

  RooRealVar fullSIM_m0("fullSIM_m0","",90,92);
  RooRealVar fullSIM_sig("fullSIM_sig","",2.495);
  RooRealVar fullSIM_alpha("fullSIM_alpha","",0,10);
  RooRealVar fullSIM_n("fullSIM_n","",0,10);

  RooRealVar fullSIM_smear_mean("fullSIM_smear_mean","",0);
  RooRealVar fullSIM_smear_sigma("fullSIM_smear_sigma","",0,5);

  RooCBShape fullSIM_cb("fullSIM_cb","",mgg,fullSIM_m0,fullSIM_sig,fullSIM_alpha,fullSIM_n);
  RooGaussian fullSIM_smear("fullSIM_smear","",mgg,fullSIM_smear_mean,fullSIM_smear_sigma);

  RooFFTConvPdf fullSIM_conv("fullSIM_conv","",mgg,fullSIM_cb,fullSIM_smear);

  RooRealVar fullSIM_bkg_alpha("fullSIM_bkg_alpha","",-3,0);

  RooExponential fullSIM_bkg("fullSIM_bkg","",mgg,fullSIM_bkg_alpha);

  RooRealVar fullSIM_N_sig("fullSIM_N_sig","",100,1e9);
  RooRealVar fullSIM_N_bkg("fullSIM_N_bkg","",10,500);

  RooAddPdf fullSIM_fitModel("fullSIM_fitModel","",RooArgList(fullSIM_conv,fullSIM_bkg),RooArgList(fullSIM_N_sig,fullSIM_N_bkg));

  fullSIM_fitModel.fitTo(fullSIM_Data,RooFit::Strategy(0));
  fullSIM_fitModel.fitTo(fullSIM_Data,RooFit::Strategy(2));

  RooPlot *fmgg = mgg.frame();
  fullSIM_Data.plotOn(fmgg);
  fullSIM_fitModel.plotOn(fmgg);
  
  TCanvas fullSIM_cv("fullSIM_cv","FullSIM");
  fmgg->Draw();


  cout << "Smearing: \nFastSIM: " << fastSIM_smear_sigma.getVal() << std::endl
       <<"FullSIM: " << fullSIM_smear_sigma.getVal() << std::endl;

}
