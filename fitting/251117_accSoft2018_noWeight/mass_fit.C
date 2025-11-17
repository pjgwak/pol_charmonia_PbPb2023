using namespace RooFit;

void mass_fit(float ptLow = 20, float ptHigh = 50, float yLow = 0, float yHigh = 1.6, bool isSkipFit = false)
{
	cout << "=== start mass_fit() ===\n";
  TStopwatch time;
  time.Start();

  // time check helper
  auto output_current_time = []() {
    auto now = chrono::system_clock::now();
    time_t current_time = chrono::system_clock::to_time_t(now);
    tm* local_tm = localtime(&current_time);
    cout << "Current time: " << put_time(local_tm, "%c") << "\n";
  };

	// HI data taking started from 2 p.m today

  // === set kinematics ===
  // --- kinematics ---
  // float ptLow = 3, ptHigh = 50;
  // float yLow = 1.6, yHigh = 2.4;
  float massLow = 2.6, massHigh = 3.5;
  int cLow = 0, cHigh = 180;

  // --- initial observable range ---
  // double ctLow = -10, ctHigh = 10;
  // float ctErrLow = 0, ctErrHigh = 0.99;

  // --- mass bkg parameters ---
  int bkgMassOrder = 2; // must be same with the order or Chebychev (check codes defining mass bkg model)
  double sbL_lo = 2.6, sbL_hi = 2.9;
  double sbR_lo = 3.3, sbR_hi = 3.5;
  double sig_lo = 2.9, sig_hi = 3.3;


  // === silence ===
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
  // RooMsgService::instance().getStream(0).removeTopic(InputArguments);
  // RooMsgService::instance().getStream(1).removeTopic(InputArguments);
  // RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Integration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Integration);
  // RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING); // only print from WARNING to FATAL


  // === set output properies ===
  // --- make output folders ---
  TString userLabel = "";
  TString figDir = Form("figs%s/pT%.1f_%.1f_y%.1f_%.1f", userLabel.Data(), ptLow, ptHigh, yLow, yHigh);
  gSystem->mkdir(figDir, kTRUE);
  TString rootDir = Form("roots%s/pT%.1f_%.1f_y%.1f_%.1f", userLabel.Data(), ptLow, ptHigh, yLow, yHigh);
  gSystem->mkdir(rootDir, kTRUE);
  gSystem->mkdir(Form("logs%s", userLabel.Data()), kTRUE);


  // === set cosmetics - rootlogon ===
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");


  // === read inputs ===
  cout << "\n=== import inputs ===\n";

  // --- data ---
  string fileNameData = "/data/users/pjgwak/work/pol_PbPb2023/skimming/roodataset_files/RooDataSet_miniAOD_isMC0_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_PbPb23_MinBias.root";
  TFile fInData(fileNameData.c_str());
  // cout << fileNameData.c_str() << endl;
  RooDataSet *data = (RooDataSet *)fInData.Get("dataset");
  data->SetName("data");


  // === created workspace ===
  RooWorkspace *ws = new RooWorkspace("ws");


  // === set cuts 1 ===
  cout << "\n=== make cuts 1 ===\n";
  string reduceDS_woCtErr = Form("(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f)",  ptLow, ptHigh, yLow, yHigh);

	// string reduceDS_woCtErr = Form("(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f && ctau3D >= %.3f && ctau3D < %.3f)",  ptLow, ptHigh, yLow, yHigh, massLow, massHigh, ctLow, ctHigh);


  // === reduce dataset 1 ===
  cout << "\n=== reduce dataset 1 ===\n";
  RooDataSet *redData = (RooDataSet *)data->reduce(reduceDS_woCtErr.c_str()); // smae cut with Data
  redData->SetName("redData");
  ws->import(*redData);


  // === mc mass fit ===
  // --- set observables ---
  ws->var("mass")->setRange(2.6, 3.5);
  ws->var("mass")->setRange("massRange", 2.6, 3.5);

  // --- build mass signal model ---
  RooRealVar massMean("massMean", "", 3.096, 3.07, 3.102);
  RooRealVar massSigma1("massSigma1", "CB sigma", 0.030, 0.001, 0.200);
  RooRealVar massSigma12("massSigma12", "ratio of sigma 1 vs 2", 1.5, 1, 10);
  RooFormulaVar massSigma2("massSigma2", "@0*@1", {massSigma1, massSigma12});

  RooRealVar alphaL("alphaL", "", 1.5, 0.001, 10.0);
  RooRealVar nL("nL", "nL", 3.0, 1, 20.0);

	// RooRealVar alphaR("alphaR", "", 1.0, 0, 10.0);
  // RooRealVar nR("nR", "", 2.0, 1, 10.0);
  
  RooRealVar alphaLR("alphaLR", "", 1.0, 0, 20.0);
  RooRealVar nLR("nLR", "", 2.0, 1.0, 20.0);
  RooFormulaVar alphaR("alphaR", "@0*@1", {alphaL, alphaLR});
  RooFormulaVar nR("nR", "@0*@1", {nL, nLR});

  // RooFormulaVar alphaR("alphaR", "@0", {alphaL});
  // RooFormulaVar nR("nR", "@0", {nL});

  // RooRealVar alphaR("alphaR", "alphaR", -1.5, -5.0, -0.2)
  // RooRealVar nR("nR", "nR", 3.0, 1.0, 20.0);

  RooRealVar fCB1("fCB1", "frac(Gauss in total)", 0.50, 0.00, 1.00);
	
	RooCBShape CB1("CB1", "left-tail CB", *ws->var("mass"), massMean, massSigma1, alphaL, nL);

  // Gauss
  RooRealVar massSigma1G("massSigma1G", "ratio of sigma1 vs sigmaG", 1.1, 1, 10);
  RooFormulaVar massSigmaG("massSigmaG", "@0*@1", {massSigma1, massSigma1G});
  RooGaussian massG("massG", "core Gauss", *ws->var("mass"), massMean, massSigmaG);
  
  RooCrystalBall DCB("DCB", "", *ws->var("mass"), massMean, massSigma1, alphaL, nL, alphaR, nR);

  RooAddPdf massSigPdf("massSigPdf", "", RooArgList(DCB, massG), RooArgList(fCB1)); // Recursive fracton


  // --- build mass bkg model ---
  RooRealVar massSl1("massSl1", "massSl1", 0.01, -1.0, 1.0);
  RooRealVar massSl2("massSl2", "massSl2", 0.01, -1.0, 1.0);
  RooRealVar massSl3("massSl3", "massSl3", 0.01, -1.0, 1.0);
  // RooChebychev massBkgPdf("massBkgPdf", "", *ws->var("mass"), RooArgList(massSl1, massSl2));
  // RooChebychev massBkgPdf("massBkgPdf", "", *ws->var("mass"), RooArgList(massSl1, massSl2, massSl3)); //7.5-8.5
  RooChebychev massBkgPdf("massBkgPdf", "", *ws->var("mass"), RooArgList(massSl1)); // low pt: 9.5-11, 11-13

  // --- fit ---
  RooRealVar nSigMass("nSigMass", "mass signal yield", 8e4, 1, 1e8);
  RooRealVar nBkgMass("nBkgMass", "mass bkg yield", 1e5, 1, 1e8);
  RooAddPdf massFitModel("massFitModel", "",
                  RooArgList(massSigPdf, massBkgPdf),
                  RooArgList(nSigMass, nBkgMass));
  
  alphaL.setConstant();
  nL.setConstant();
  nLR.setConstant();
  alphaLR.setConstant();
  massSigma1G.setConstant();

  // --- helper function ---
    auto syncVar = [](RooRealVar& v, const RooFitResult& fr) {
      auto findIn = [&](const RooArgList& lst) -> const RooRealVar* {
        if (auto* a = lst.find(v.GetName())) return dynamic_cast<const RooRealVar*>(a);
        return nullptr;
      };
      const RooRealVar* src = nullptr;
      if (!(src = findIn(fr.floatParsFinal())))
          src = findIn(fr.constPars());
      if (!src) return false;
      v.setVal(src->getVal());
      v.setError(src->getError());
      return true;
    };

  // --- use mc result ---
  TFile finMC(Form("%s/mc_mass_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "READ");
  RooFitResult* fitMcMass = nullptr; finMC.GetObject("fitMcMass", fitMcMass);
  if (fitMcMass) fitMcMass = (RooFitResult*) fitMcMass->Clone("fitMcMass");
  
  syncVar(alphaL, *fitMcMass); 
  syncVar(alphaLR, *fitMcMass);
  syncVar(nL, *fitMcMass);
  syncVar(nLR, *fitMcMass);
  syncVar(massSigma1G, *fitMcMass);

  // --- fit ---
  RooFitResult *fitMass;
  if (isSkipFit && !gSystem->AccessPathName(Form("%s/mass_fit_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), kReadPermission)) {
    TFile fin(Form("%s/mass_fit_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "READ");
    RooFitResult* tmp = nullptr; fin.GetObject("fitMass", tmp);
    if (tmp) fitMass = (RooFitResult*) tmp->Clone("fitMass");

    syncVar(alphaL, *fitMass); 
    syncVar(alphaLR, *fitMass);
    syncVar(fCB1, *fitMass);
    syncVar(massMean, *fitMass);
    syncVar(massSigma1, *fitMass);
    syncVar(massSigma1G, *fitMass);
    syncVar(nL, *fitMass);
    syncVar(nLR, *fitMass);
    syncVar(nSigMass, *fitMass);
    syncVar(nBkgMass, *fitMass);
    syncVar(massSl1, *fitMass);
    syncVar(massSl2, *fitMass);
    syncVar(massSl3, *fitMass);
  } 
  else
    fitMass = massFitModel.fitTo(*redData, Offset(true), Save(), Extended(), PrintLevel(-1), PrintEvalErrors(-1), RecoverFromUndefinedRegions(2), NumCPU(32), EvalBackend("legacy")); //, EvalBackend("legacy")

	// alphaL.setConstant();

	// fitMass = massFitModel.fitTo(*redData, Offset(true), Save(), Extended(), PrintLevel(-1), PrintEvalErrors(-1), RecoverFromUndefinedRegions(2), NumCPU(4));
	
  // --- draw plots ---
	auto findObj = [&](RooPlot *fr, const char *n) -> TObject *
  { return fr ? fr->findObject(n) : nullptr; };
  {
    TCanvas c("c_mass", "c_mass", 800, 800);
    TPad pad1("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
    pad1.SetBottomMargin(0.00001);
    pad1.SetLogy();
    pad1.Draw();
    pad1.cd();

    // --- frame & plot ---
    RooPlot *fr = ws->var("mass")->frame(Range(massLow, massHigh), Title("")); // , Bins(nBins)
    redData->plotOn(fr, DataError(RooAbsData::SumW2), Name("data"));
    massFitModel.plotOn(fr, Range("massRange"), Name("model"));
    massFitModel.plotOn(fr, Range("massRange"), Components("DCB"), Name("DCB"), LineStyle(kDotted), LineColor(kRed));
    massFitModel.plotOn(fr, Range("massRange"), Components("massG"), Name("G"), LineStyle(kDotted), LineColor(kGreen));
    massFitModel.plotOn(fr, Range("massRange"), Components("massBkgPdf"), Name("bkg"), LineStyle(kDotted), LineColor(kAzure));

    // --- dynamic y-range for log scale ---
    double ymin = 1e300, ymax = -1e300;
    if (auto *hdata = dynamic_cast<RooHist *>(fr->getHist("data")))
    {
      for (int i = 0; i < hdata->GetN(); ++i)
      {
        double x, y;
        hdata->GetPoint(i, x, y);
        if (y > 0 && y < ymin)
          ymin = y;
        if (y > ymax)
          ymax = y;
      }
    }
    if (ymin <= 0 || ymin == 1e300)
      ymin = 1e-3;
    fr->SetMinimum(ymin * 0.5);
    fr->SetMaximum(std::max(ymax, ymin) * 1e3);

    fr->GetYaxis()->SetTitle("Events");
    fr->GetXaxis()->SetTitle("");
    fr->Draw("e");

    // --- legend ---
    TLegend leg(0.49, 0.66, 0.70, 0.94);
    {
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextSize(0.03);
      if (auto *o = findObj(fr, "data"))
        leg.AddEntry(o, "Data", "lep");
      if (auto *o = findObj(fr, "model"))
        leg.AddEntry(o, "Model", "pe");
      if (auto *o = findObj(fr, "DCB"))
        leg.AddEntry(o, "DCB", "pe");
      if (auto *o = findObj(fr, "G"))
        leg.AddEntry(o, "Gauss", "pe");
      if (auto *o = findObj(fr, "bkg"))
        leg.AddEntry(o, "Bkg", "pe");

      leg.Draw("same");
    }

    // --- CMS/info latex ---
    {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.03);
      tx.SetTextFont(42);
      double x = 0.19, y0 = 0.90, dy = -0.06;
      int k = 0;
      tx.DrawLatex(x, y0 + dy * k++, "CMS PbPb #sqrt{s_{NN}} = 5.32 TeV");
      tx.DrawLatex(x, y0 + dy * k++, "2023, J/#psi #rightarrow #mu^{+}#mu^{-}");
      if (yLow == 0)
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
      else
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));
      
      // fit status
      int st = fitMass->status(); // 0 = success
      if (st != 0)
      {
        tx.DrawLatex(x, y0 + dy * k++, Form("Status=%d", st));
        std::ofstream flog(Form("logs%s/mass_status_pT%.1f_%.1f.txt", userLabel.Data(), ptLow, ptHigh), std::ios::app);
        flog.close();
      }
    }

    // --- parameter latex ---
    {
      TLatex tp;
      tp.SetNDC();
      tp.SetTextSize(0.024);
      tp.SetTextFont(42);
      double x = 0.71, y0 = 0.91, dy = -0.04;
      int k = 0;

      // lambda function for printing
      auto print = [&](const char *title, const char *vname, RooAddPdf &model)
      {
        auto *v = dynamic_cast<RooRealVar *>(model.getVariables()->find(vname));
        if (!v)
          return;

        const double val = v->getVal(), err = v->getError();
        const double eps = 1e-9 * (1.0 + std::fabs(val));
        const bool fixed = v->isConstant();
        const bool atMin = v->hasMin() && std::fabs(val - v->getMin()) <= eps;
        const bool atMax = v->hasMax() && std::fabs(val - v->getMax()) <= eps;
        const bool atBound = atMin || atMax;

        TString note;
        if (fixed)
          note += "(fixed)";
        if (atBound)
          note += note.IsNull() ? (atMin ? "(at min)" : "(at max)")
                                : (atMin ? ", (at min)" : ", (at max)");

        const Int_t oldColor = tp.GetTextColor();
        if (atBound)
          tp.SetTextColor(kRed + 1);

        if (fixed)
          tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g %s", title, val, note.Data()));
        else
          tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g #pm %.3g %s", title, val, err, note.Data()));

        tp.SetTextColor(oldColor);
      };

      // print
      print("N_{Sig}", "nSigMass", massFitModel);
      print("N_{Bkg}", "nBkgMass", massFitModel);

      print("mean", "massMean", massFitModel);
      print("#alpha_{L}", "alphaL", massFitModel);
      print("n_{L}", "nL", massFitModel);
      print("#alpha_{L/R}", "alphaLR", massFitModel);
      print("n_{L/R}", "nLR", massFitModel);

      print("#sigma_{1}", "massSigma1", massFitModel);
      // print("#sigma_{2/1}", "massSigma12", massFitModel);
      print("#sigma_{G/1}", "massSigma1G", massFitModel);
      
      print("f_{CB1}", "fCB1", massFitModel);
      // print("f_{CB2}", "fCB2", massFitModel);
      
      print("sl1", "massSl1", massFitModel);
      print("sl2", "massSl2", massFitModel);
      print("sl3", "massSl3", massFitModel);
    }

    // --- pull pad ---
    c.cd();
    TPad pad2("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
    pad2.SetTopMargin(0.00001);
    pad2.SetBottomMargin(0.4);
    pad2.Draw();
    pad2.cd();

    RooHist *hpull = fr->pullHist("data", "model");
    RooPlot *fpull = ws->var("mass")->frame(Range(massLow, massHigh), Title(""));
    fpull->addPlotable(hpull, "P");
    fpull->GetYaxis()->SetTitle("Pull");
    fpull->GetXaxis()->SetTitle("mass^{inv}_{#mu#mu} [GeV/c^{2}]");
    fpull->GetXaxis()->CenterTitle();
    fpull->SetMinimum(-8);
    fpull->SetMaximum(8);
    fpull->GetYaxis()->SetNdivisions(505);
    fpull->GetYaxis()->SetTitleSize(0.12);
    fpull->GetYaxis()->SetLabelSize(0.10);
    fpull->GetXaxis()->SetTitleSize(0.15);
    fpull->GetXaxis()->SetLabelSize(0.10);
    fpull->Draw();

    TLine line(massLow, 0.0, massHigh, 0.0);
    line.SetLineStyle(2);
    line.Draw("same");

    // --- chi2/ndf ---
    if (fitMass)
    {
      int npar = fitMass->floatParsFinal().getSize();
      double chi2ndf = fr->chiSquare("model", "data", npar);
      TLatex tc;
      tc.SetNDC();
      tc.SetTextSize(0.10);
      tc.DrawLatex(0.82, 0.86, Form("#chi^{2}/ndf = %.2f", chi2ndf));
      cout << "\nchi2/ndf = " << chi2ndf << "\n\n";
    }
    c.SaveAs(Form("%s/mass_fit_pT%.1f_%.1f.png", figDir.Data(), ptLow, ptHigh)); 
  }

  // --- save results ---
  ws->import(massFitModel);

  fitMass->Print("V");
  TFile outMass(Form("%s/mass_fit_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "RECREATE");
  fitMass->Write("fitMass");
  ws->Write("wsMass");
	// massFitModel.Write();
  outMass.Close();

	cout << "\n=== finish mass_fit() ===\n";
  time.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", time.RealTime(), time.CpuTime());
}