

void lambda_fit(float ptLow = 20, float ptHigh = 50, float yLow = 0, float yHigh = 1.6, bool isSkipFit = true)
{
  using namespace RooFit;

  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  TString userLabel = "";
  TString figDir = Form("figs%s/pT%.1f_%.1f_y%.1f_%.1f", userLabel.Data(), ptLow, ptHigh, yLow, yHigh);
  TString rootDir = Form("roots%s/pT%.1f_%.1f_y%.1f_%.1f", userLabel.Data(), ptLow, ptHigh, yLow, yHigh);

  TString outFile  = "nSigMass_vs_cosHX.root";


  TH1D* h_nSig = new TH1D("h_nSig", "nSigMass vs cosHX;cosHX;N_{sig}", 10, 0.0, 1.0);

  for (int i = 0; i <= 8; ++i) {
    double cosLow  = 0.1 * i;
    double cosHigh = cosLow + 0.1;

    TString fileName = Form("%s/mass_fit_pT%.1f_%.1f_cosHX%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh, cosLow, cosHigh);

    std::cout << "[INFO] opening: " << fileName << std::endl;

    TFile *f = TFile::Open(fileName, "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "[WARN] failed to open " << fileName << std::endl;
      if (f) f->Close();
      continue;
    }

    RooFitResult *fitMass = nullptr;
    f->GetObject("fitMass", fitMass);

    if (!fitMass) {
      std::cerr << "[WARN] fitMass not found in " << fileName << std::endl;
      f->Close();
      delete f;
      continue;
    }

    // fitMass 안에서 floating parameters 중 nSigMass 찾기
    RooRealVar *nSigMass =
      dynamic_cast<RooRealVar*>( fitMass->floatParsFinal().find("nSigMass") );

    if (!nSigMass) {
      std::cerr << "[WARN] nSigMass not found in fitMass (file: "
                << fileName << ")" << std::endl;
      f->Close();
      delete f;
      continue;
    }

    double val  = nSigMass->getVal();
    double err  = nSigMass->getError();
    // double errHi = nSigMass->getAsymErrorHi();
    // double errLo = nSigMass->getAsymErrorLo();

    // std::cout << Form("cosHX %.1f–%.1f : nSigMass = %.3f  +/- %.3f"
    //                   "  (asym: -%.3f / +%.3f)",
    //                   cosLow, cosHigh, val, err)
    //           << std::endl;

    // find bin posiion
    double cosMid = 0.5 * (cosLow + cosHigh);
    int bin = h_nSig->FindBin(cosMid);

    h_nSig->SetBinContent(bin, val);
    h_nSig->SetBinError(bin,   err);

    f->Close();
    delete f;
  }

  RooRealVar  costh("costh", "cos#theta", -1.0, 1.0);
  RooRealVar  lamth("lamth", "#lambda_{#theta}", 0.0, -1.0, 1.0);
  RooGenericPdf thetaPdf(
    "thetaPdf",
    "1.0 + lamth*costh*costh",
    RooArgList(costh, lamth)
  );

  RooDataHist dh("dh", "dh", RooArgList(costh), Import(*h_nSig));

  costh.setRange("fitRange", 0, 0.7);


  // === fit ===
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
  
  RooFitResult *fitResult;
  if (isSkipFit && !gSystem->AccessPathName(Form("%s/lambda_fit_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), kReadPermission)) {
    TFile fin(Form("%s/lambda_fit_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "READ");
    RooFitResult* tmp = nullptr; fin.GetObject("fitResult", tmp);
    if (tmp) fitResult = (RooFitResult*) tmp->Clone("fitResult");

    syncVar(lamth, *fitResult); 
  } 
  else fitResult = thetaPdf.fitTo(dh, Range("fitRange"), Save());
  fitResult->Print("V");


  auto findObj = [&](RooPlot *fr, const char *n) -> TObject *
    { return fr ? fr->findObject(n) : nullptr; };
    {
      TCanvas c("c_mass", "c_mass", 800, 800);
      TPad pad1("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
      pad1.SetBottomMargin(0.00001);
      pad1.Draw();
      pad1.cd();

      // --- frame & plot ---
      RooPlot *fr = costh.frame(Title("")); // , Bins(nBins)
      dh.plotOn(fr, DataError(RooAbsData::SumW2), Name("data"));
      thetaPdf.plotOn(fr, Name("model"));

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
      // fr->SetMinimum(0.0000001); // small number for osmetics
      fr->SetMaximum(std::max(ymax, ymin) * 2.5);

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
        // if (auto *o = findObj(fr, "DCB"))
        //   leg.AddEntry(o, "DCB", "pe");

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
        int st = fitResult->status(); // 0 = success
        if (st != 0)
        {
          tx.DrawLatex(x, y0 + dy * k++, Form("Status=%d", st));
          std::ofstream flog(Form("logs%s/lambda_status_pT%.1f_%.1f.txt", userLabel.Data(), ptLow, ptHigh), std::ios::app);
          flog.close();
        }
      }

      // --- parameter latex ---
      {
        TLatex tp;
        tp.SetNDC();
        tp.SetTextSize(0.03);
        tp.SetTextFont(42);
        double x = 0.71, y0 = 0.91, dy = -0.04;
        int k = 0;

        // lambda function for printing
        auto print = [&](const char *title, const char *vname, RooAbsPdf &model)
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
        print("#lambda_{#theta, HX}", "lamth", thetaPdf);
      }

      // --- pull pad ---
      c.cd();
      TPad pad2("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
      pad2.SetTopMargin(0.00001);
      pad2.SetBottomMargin(0.4);
      pad2.Draw();
      pad2.cd();

      RooHist *hpull = fr->pullHist("data", "model");
      RooPlot *fpull = costh.frame(Title(""));//Range("fitRange"), 
      fpull->addPlotable(hpull, "P");
      fpull->GetYaxis()->SetTitle("Pull");
      fpull->GetXaxis()->SetTitle("cos#theta_{HX}");
      fpull->GetXaxis()->CenterTitle();
      fpull->SetMinimum(-8);
      fpull->SetMaximum(8);
      fpull->GetYaxis()->SetNdivisions(505);
      fpull->GetYaxis()->SetTitleSize(0.12);
      fpull->GetYaxis()->SetLabelSize(0.10);
      fpull->GetXaxis()->SetTitleSize(0.15);
      fpull->GetXaxis()->SetLabelSize(0.10);
      fpull->Draw();

      TLine line(0, 0.0, 1, 0.0);
      line.SetLineStyle(2);
      line.Draw("same");

      // --- chi2/ndf ---
      if (fitResult)
      {
        int npar = fitResult->floatParsFinal().getSize();
        double chi2ndf = fr->chiSquare("model", "data", npar);
        TLatex tc;
        tc.SetNDC();
        tc.SetTextSize(0.10);
        tc.DrawLatex(0.82, 0.86, Form("#chi^{2}/ndf = %.2f", chi2ndf));
        cout << "\nchi2/ndf = " << chi2ndf << "\n\n";
      }
      c.SaveAs(Form("%s/lambda_fit_pT%.1f_%.1f_cosHX.png", figDir.Data(), ptLow, ptHigh));
    }

  TFile* fout = TFile::Open(Form("%s/lambda_fit_pT%.1f_%.1f_cosHX.root", rootDir.Data(), ptLow, ptHigh), "RECREATE");
  h_nSig->Write();
  fitResult->Write("fitResult");
  fout->Close();
}

// lambda fit
  // print