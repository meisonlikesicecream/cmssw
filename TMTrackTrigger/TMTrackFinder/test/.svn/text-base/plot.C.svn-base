{
// In unnamed scripts, variables not forgotten at end, so must delete them before rerunning script, so ...
  gROOT->Reset("a");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat("");
  //gStyle->SetOptStat("emr");
  //  gStyle->SetOptStat("euom");
  gStyle->SetStatFontSize(0.035);
  gStyle->SetHistFillColor(kBlue);
  gStyle->SetHistFillStyle(1001);

  gStyle->SetStatFormat("5.3f");
  gStyle->SetStatFontSize(0.04);
  gStyle->SetOptFit(0111);
  gStyle->SetStatW(0.30);
  gStyle->SetStatH(0.02);
  gStyle->SetStatX(0.5);
  gStyle->SetStatY(0.9);
  gStyle->SetTitleYOffset(1.25);
  gStyle->SetTitleSize(0.05, "XYZ");

  TCanvas d1("d1");

  TH1F* his;
  TProfile* prof;
  TEfficiency* teffi;

  TFile fileSLR("out_LR_davidenew_noDigiSLR_noDigi_20171012_153941/Hist.root");
  TFile fileKF("out_KF_ChosenR612_HTcellMod_20170920_112548/Hist.root");

  fileSLR.GetObject("TMTrackProducer/Effi_SimpleLR/AlgEffFitVsInvPt_SimpleLR", teffi);
  teffi->Draw();
  teffi->SetTitle(";1/Pt (1/GeV);Algorithmic Efficiency");
  d1.Draw(); d1.Update(); 
  d1.Print("Plots/AlgEffFitVsInvPt_SLR.pdf");
  cin.get(); 

  fileKF.GetObject("TMTrackProducer/Effi_KF4ParamsComb/AlgEffFitVsInvPt_KF4ParamsComb", teffi);
  teffi->Draw();
  teffi->SetTitle(";1/Pt (1/GeV);Algorithmic Efficiency");
  d1.Draw(); d1.Update(); 
  d1.Print("Plots/AlgEffFitVsInvPt_KF.pdf");
  cin.get(); 

  fileSLR.GetObject("TMTrackProducer/Effi_SimpleLR/AlgEffFitVsEta_SimpleLR", teffi);
  teffi->Draw();
  teffi->SetTitle(";#eta;Algorithmic Efficiency");
  d1.Draw(); d1.Update(); 
  d1.Print("Plots/AlgEffFitVsEta_SLR.pdf");
  cin.get(); 

  fileKF.GetObject("TMTrackProducer/Effi_KF4ParamsComb/AlgEffFitVsEta_KF4ParamsComb", teffi);
  teffi->Draw();
  teffi->SetTitle(";#eta;Algorithmic Efficiency");
  d1.Draw(); d1.Update(); 
  d1.Print("Plots/AlgEffFitVsEta_KF.pdf");
  cin.get(); 

  fileSLR.GetObject("TMTrackProducer/SimpleLR/Z0ResVsTrueEta_SimpleLR", prof);
  prof->Draw();
  prof->SetTitle(";#eta;z0 resolution (cm)");
  d1.Draw(); d1.Update(); 
  d1.Print("Plots/Z0resVsEta_SLR.pdf");
  cin.get(); 

  fileKF.GetObject("TMTrackProducer/KF4ParamsComb/Z0ResVsTrueEta_KF4ParamsComb", prof);
  prof->Draw();
  prof->SetTitle(";#eta;z0 resolution (cm)");
  d1.Draw(); d1.Update(); 
  d1.Print("Plots/Z0resVsEta_KF.pdf");
  cin.get(); 

  //float ymax = BendResStub.GetMaximum();
  //BendResStub.SetMaximum(1.9*ymax);

  //his.SetXTitle("fun");

  fileSLR.Close();
  fileKF.Close();
}
