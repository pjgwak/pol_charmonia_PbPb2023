// For PbPb2023 data

#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <string>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "cutsAndBins.h"
using std::cout; using std::string;

static const long MAXTREESIZE = 1000000000000;

TVector3 MuPlusVector_Helicity(const TLorentzVector &QQLV_Lab, const TLorentzVector &MuPlusLV_Lab);
TVector3 MuPlusVector_CollinsSoper(const TLorentzVector &QQLV_Lab, const TLorentzVector &MuPlusLV_Lab);

// 
bool muonAccSoft2018(double pt, double eta)
{
  return ((fabs(eta) < 0.3 && pt > 3.4) ||
          (0.3 < fabs(eta) && fabs(eta) < 1.1 && pt > 3.3) ||
          (1.1 < fabs(eta) && fabs(eta) < 1.5 && pt > 9.08 - 5.25*fabs(eta)) ||
          (1.5 < fabs(eta) && fabs(eta) < 2.4 && pt > 0.8 && pt>2.4-0.8*fabs(eta)));
}

void onia_to_skim_data(int myTrig=24, bool isMC = false, long nEvt = -1)
{
  cout << "Start onia_to_skim_data()\n";

  // ===== read input =====
  TFile *fInput = TFile::Open("/data/Oniatree/polarization/oniatree_5p36/PbPb2023_Data_noTrack/Oniatree_2023PbPbPromptRecoData_132X_miniAOD_noTrack.root", "read");
  TTree *oniaTree = (TTree *)fInput->Get("hionia/myTree");

  // user defined kinematic cuts
  int cLow = 0, cHigh = 180;
  double massLow = 2.6, massHigh = 3.5;

  // ===== labeling =====
  // bool isPr=true
  // string mcLabel = "";
  // if (isMC==true, isPr==true) mcLabel = "PR";
  // else if (isMC==true, isPr==false) mcLabel = "NP";
  // skip

  // ===== set Oniatree branch address (input) =====
  const long int maxBranchSize = 1000;

  // event-level scalars
  UInt_t runNb;
  UInt_t eventNb;
  UInt_t LS;
  float zVtx;
  Int_t Centrality;
  ULong64_t HLTriggers;
  Float_t SumET_HF;
  Float_t Gen_weight; // MC

  // collection sizes
  Short_t Reco_QQ_size;
  Short_t Reco_mu_size;

  // TClonesArray pointer
  TClonesArray *Reco_QQ_4mom = nullptr;
  TClonesArray *Reco_mu_4mom = nullptr;

  // per-dimuon (size = Reco_QQ_size)
  ULong64_t Reco_QQ_trig[maxBranchSize];
  Float_t Reco_QQ_VtxProb[maxBranchSize];
  Short_t Reco_QQ_mupl_idx[maxBranchSize];
  Short_t Reco_QQ_mumi_idx[maxBranchSize];
  Short_t Reco_QQ_sign[maxBranchSize];
  Float_t Reco_QQ_ctau3D[maxBranchSize];
  Float_t Reco_QQ_ctauErr3D[maxBranchSize];

  // per-muon (size = Reco_mu_size)
  ULong64_t Reco_mu_trig[maxBranchSize];
  Bool_t Reco_mu_highPurity[maxBranchSize];
  Int_t Reco_mu_nTrkHits[maxBranchSize];
  Float_t Reco_mu_normChi2_global[maxBranchSize];
  Int_t Reco_mu_nMuValHits[maxBranchSize];
  Int_t Reco_mu_StationsMatched[maxBranchSize];
  Float_t Reco_mu_dxy[maxBranchSize];
  Float_t Reco_mu_dxyErr[maxBranchSize];
  Float_t Reco_mu_dz[maxBranchSize];
  Float_t Reco_mu_dzErr[maxBranchSize];
  Int_t Reco_mu_nTrkWMea[maxBranchSize];
  Int_t Reco_mu_nPixWMea[maxBranchSize];
  Int_t Reco_mu_nPixValHits[maxBranchSize];
  Int_t Reco_mu_SelectionType[maxBranchSize];
  Int_t Reco_mu_whichGen[maxBranchSize]; // MC

  Bool_t Reco_mu_isSoftCutBased[maxBranchSize];
  Bool_t Reco_mu_isGlobal[maxBranchSize];
  Bool_t Reco_mu_isTracker[maxBranchSize];

  // ----- SetBranchAddress -----
  oniaTree->SetBranchAddress("runNb", &runNb);
  oniaTree->SetBranchAddress("eventNb", &eventNb);
  oniaTree->SetBranchAddress("LS", &LS);
  oniaTree->SetBranchAddress("zVtx", &zVtx);
  oniaTree->SetBranchAddress("Centrality", &Centrality);
  oniaTree->SetBranchAddress("HLTriggers", &HLTriggers);
  oniaTree->SetBranchAddress("SumET_HF", &SumET_HF);

  // sizes
  oniaTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
  oniaTree->SetBranchAddress("Reco_mu_size", &Reco_mu_size);

  // TLorentzVecotr
  oniaTree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom);
  oniaTree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom);

  // dimuon
  oniaTree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig);
  oniaTree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb);
  oniaTree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
  oniaTree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);
  oniaTree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign);
  oniaTree->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D);
  oniaTree->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D);

  // single-muon
  oniaTree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig);
  oniaTree->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity);
  oniaTree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits);
  oniaTree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global);
  oniaTree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits);
  oniaTree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched);
  oniaTree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy);
  oniaTree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr);
  oniaTree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz);
  oniaTree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr);
  oniaTree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea);
  oniaTree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea);
  oniaTree->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits);
  oniaTree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType);

  oniaTree->SetBranchAddress("Reco_mu_isSoftCutBased", Reco_mu_isSoftCutBased);
  oniaTree->SetBranchAddress("Reco_mu_isGlobal", Reco_mu_isGlobal);
  oniaTree->SetBranchAddress("Reco_mu_isTracker", Reco_mu_isTracker);

  // MC only
  if (isMC)
  {
    oniaTree->SetBranchAddress("Gen_weight", &Gen_weight);
    oniaTree->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen);
  }

  // ===== set output file and  FlowSkim branch address (output) =====
    // ----- bit mask - change here -----
  // int trigIndx = 0;
  // int kTrigSel = myTrig; // 24:HLT_HIMinimumBiasHF1ANDZDC1nOR_v
  // if (kTrigSel == kTrigJpsipp) trigIndx = 0; // PbPb : 0
  // else if (kTrigSel == kTrigUps) trigIndx = 1;
  // else if (kTrigSel == kTrigL1DBOS40100) trigIndx = 2;
  // else if (kTrigSel == kTrigL1DB50100) trigIndx = 3;

  // int kL2filter = 16;
  // int kL3filter = 17;

  // const ULong64_t HLT_MASK = (1ULL << kTrigSel);
  // const ULong64_t L2_MASK = (1ULL << kL2filter);
  // const ULong64_t L3_MASK = (1ULL << kL3filter);

  // ----- declare local variables ----
  const static long long int nMaxDimu = 1000;
  
  // event-level scalar
  Int_t event = 0;
  Int_t runN = 0;
  Int_t lumi = 0;
  Float_t vz = 0;

  // dimuon size
  Int_t nDimu = 0;

  // dimuon candidate-level variables (size = nDimu)
  Float_t mass[nMaxDimu];
  Float_t y[nMaxDimu];
  Float_t pt[nMaxDimu], pt1[nMaxDimu], pt2[nMaxDimu];
  Float_t eta[nMaxDimu], eta1[nMaxDimu], eta2[nMaxDimu];
  Float_t phi[nMaxDimu], phi1[nMaxDimu], phi2[nMaxDimu];
  Int_t recoQQsign[nMaxDimu], cBin[nMaxDimu];
  Float_t ctau3D[nMaxDimu], ctau3DErr[nMaxDimu], ctau3DRes[nMaxDimu];
  Float_t phiHX[nMaxDimu], cosHX[nMaxDimu], phiCS[nMaxDimu], cosCS[nMaxDimu];

  // weights
  Double_t weight = 1.0; // event-level weight
  Double_t TnPweight[nMaxDimu]; // per-candidate level weight

  // TLorentzVector
  TLorentzVector *JP_Reco = new TLorentzVector;
  TLorentzVector *mupl_Reco = new TLorentzVector;
  TLorentzVector *mumi_Reco = new TLorentzVector;

  TFile *fFlowSkim = new TFile(Form("skim_files/flowSkim_PbPb2023_isMC%d_MinBias.root", isMC), "recreate");

  // fFlowSkim->SetCompressionSettings(207); // LZMA: 207, ZLIB: 1xx
  TTree *flowTree = new TTree("myTree", "");
  flowTree->SetMaxTreeSize(MAXTREESIZE);
  // flowTree->SetAutoSave(50 * 1024 * 1024); // save every 50 MB

  // ----- SetBranchAddress -----
  // event values
  flowTree->Branch("eventNb", &event, "eventNb/I");
  flowTree->Branch("runNb", &runN, "runNb/I");
  flowTree->Branch("LS", &lumi, "LS/I");
  flowTree->Branch("zVtx", &vz, "zVtx/F");

  // dimuon size
  flowTree->Branch("nDimu", &nDimu, "nDimu/I");

  // per-candidate (size = [nDimu])
  flowTree->Branch("mass", mass, "mass[nDimu]/F");
  flowTree->Branch("y", y, "y[nDimu]/F");
  flowTree->Branch("cBin", cBin, "cBin[nDimu]/I");
  flowTree->Branch("pt", pt, "pt[nDimu]/F");
  flowTree->Branch("pt1", pt1, "pt1[nDimu]/F");
  flowTree->Branch("pt2", pt2, "pt2[nDimu]/F");
  flowTree->Branch("eta", eta, "eta[nDimu]/F");
  flowTree->Branch("eta1", eta1, "eta1[nDimu]/F");
  flowTree->Branch("eta2", eta2, "eta2[nDimu]/F");
  flowTree->Branch("phi", phi, "phi[nDimu]/F");
  flowTree->Branch("phi1", phi1, "phi1[nDimu]/F");
  flowTree->Branch("phi2", phi2, "phi2[nDimu]/F");
  flowTree->Branch("recoQQsign", recoQQsign, "recoQQsign[nDimu]/I");
  flowTree->Branch("ctau3D", ctau3D, "ctau3D[nDimu]/F");
  flowTree->Branch("ctau3DErr", ctau3DErr, "ctau3DErr[nDimu]/F");
  flowTree->Branch("ctau3DRes", ctau3DRes, "ctau3DRes[nDimu]/F");

  // polarization angles
  flowTree->Branch("cosHX", cosHX, "cosHX[nDimu]/F");
  flowTree->Branch("phiHX", phiHX, "phiHX[nDimu]/F");
  flowTree->Branch("cosCS", cosCS, "cosCS[nDimu]/F");
  flowTree->Branch("phiCS", phiCS, "phiCS[nDimu]/F");

  // weights
  flowTree->Branch("weight", &weight, "weight/D");
  flowTree->Branch("TnPweight", TnPweight, "TnPweight[nDimu]/D");

  // ===== event loop =====
  // counters
  long count_dimuon = 0;

  // ---- preaparation before loop -----
  const Long64_t nTot = oniaTree->GetEntries();
  if (nEvt == -1)
    nEvt = nTot;
  nEvt = std::min<Long64_t>(nEvt, nTot);

  // ----- verbose option -----
  const bool VERBOSE_EVENT = false;  // event-level
  const bool VERBOSE_CANDID = false; // dimuon-level

  // ----- start event loop -----
  for (long iev = 0; iev < nEvt; ++iev)
  {
    // print progress
    if (iev % 1000000 == 0)
    {
      cout << ">>>>> EVENT " << iev << " / " << nTot
                << " (" << int(100. * iev / max<Long64_t>(1, nTot)) << "%)\n";
      cout << "Saved dimuon: " << count_dimuon << "\n";
    }
    if (VERBOSE_EVENT) {
      cout << "Total Event: " << iev << "\n";
      cout << "Total saved dimuon: " << count_dimuon << "\n";
    }

    // read event in Oniatree
    oniaTree->GetEntry(iev);

    // event-level values
    runN = runNb;
    event = eventNb;
    lumi = LS;
    vz = zVtx;

    // centrality
    int cBin_ = -999;
    if(isMC) cBin_ = Centrality;
    else if(!isMC) {
      cBin_ = getHiBinFromhiHF(SumET_HF);
      // if(hiHFBinEdge ==0) cBin_ = getHiBinFromhiHF(SumET_HF);
      // else if(hiHFBinEdge == 1) cBin_ = getHiBinFromhiHF_Up(SumET_HF);
      // else if(hiHFBinEdge == -1) cBin_ = getHiBinFromhiHF_Down(SumET_HF);
    } 
    if(cBin_==-999){ cout << "ERROR!!! No HF Centrality Matching!!" << endl; return;}
    if(cBin_ < cLow || cBin_ > cHigh) continue;
    
    // NColl - MC only, Galuber model
    double Ncoll_weight = 0, Gen_weight_ = 0;
    if(isMC){
      weight = findNcoll(Centrality) * Gen_weight;
      Ncoll_weight = findNcoll(Centrality);
      Gen_weight_ = Gen_weight;
    }

    if (TMath::Abs(vz) > 15) continue;

    // // apply HLT
    // bool HLTPass=false;
    // // if((HLTriggers&((ULong64_t)pow(2, 24))) == ((ULong64_t)pow(2, 24))
    // //   || (HLTriggers&((ULong64_t)pow(2, 25))) == ((ULong64_t)pow(2, 25))
    // //   || (HLTriggers&((ULong64_t)pow(2, 26))) == ((ULong64_t)pow(2, 26)) ) HLTPass=true;

    // if((HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel))) HLTPass=true;
    // if(HLTPass==false) continue;
    
    // check dimoun number
    if (Reco_QQ_size<0) continue;

    // check TClonesArray
    if (!Reco_QQ_4mom || !Reco_mu_4mom) continue;
    
    // ----- dimuon loop -----
    // reset number of dimuon
    nDimu = 0;

    for (Int_t irqq = 0; irqq < Reco_QQ_size; ++irqq)
    {
      if (VERBOSE_CANDID) cout << "  irqq: " << irqq << "\n";

      // check maxNDimuon
      if (nDimu >= nMaxDimu)
      {
        std::cerr << "[ERROR] nDimu reached nMaxDimu at event " << iev << "\n";
        break;
      }

      // apply HLT to dimuon candidate - not used??? - check later
      // if ((Reco_QQ_trig[irqq] & HLT_MASK) == 0ULL) continue;

      // check TLorentzVector pointers
      const int iMuPl = Reco_QQ_mupl_idx[irqq], iMuMi = Reco_QQ_mumi_idx[irqq]; // index guard
      if (iMuPl < 0 || iMuMi < 0 || iMuPl >= Reco_mu_size || iMuMi >= Reco_mu_size) continue;
      if (Reco_QQ_4mom->GetEntriesFast() <= irqq) continue;
      if (Reco_mu_4mom->GetEntriesFast() <= iMuPl || Reco_mu_4mom->GetEntriesFast() <= iMuMi) continue;

      // bring TLorentzVector
      auto JP = *static_cast<TLorentzVector *>(Reco_QQ_4mom->At(irqq));
      auto mupl = *static_cast<TLorentzVector *>(Reco_mu_4mom->At(iMuPl));
      auto mumi = *static_cast<TLorentzVector *>(Reco_mu_4mom->At(iMuMi));

      // user-defined kinematic cuts
      double mass_ = JP.M();
      if(mass_ < massLow || mass_ > massHigh) continue;

      // check MC Reco-Gen match 
      if (isMC)
      {
        if (Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1) continue;
        if (Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1) continue;
      }

      // isSoftCutBased
      if (!Reco_mu_isSoftCutBased[iMuPl] || !Reco_mu_isSoftCutBased[iMuMi]) continue;

      // dimuon y < 2.4 and single muon acceptance cut - 2018 soft
      if ( !(TMath::Abs(JP.Rapidity()) < 2.4) ||
           !muonAccSoft2018(mupl.Pt(), mupl.Eta()) ||
           !muonAccSoft2018(mumi.Pt(), mumi.Eta())) continue;

      // isGlobal - all muons in the test file pass isTracker
      if (!Reco_mu_isGlobal[iMuPl] || !Reco_mu_isGlobal[iMuMi]) continue;
      if (!Reco_mu_isTracker[iMuPl] || !Reco_mu_isTracker[iMuMi]) continue;

      // vertex probability cut
      if (Reco_QQ_VtxProb[irqq] < 0.01f) continue;

      // // selection bits
      // const bool passMuonTypePl =
      //     ((Reco_mu_SelectionType[iMuPl] & (1 << 1)) != 0) &&
      //     ((Reco_mu_SelectionType[iMuPl] & (1 << 3)) != 0);
      // const bool passMuonTypeMi =
      //     ((Reco_mu_SelectionType[iMuMi] & (1 << 1)) != 0) &&
      //     ((Reco_mu_SelectionType[iMuMi] & (1 << 3)) != 0);
      
      // // soft id cut
      // const bool muplSoft =
      //     (Reco_mu_nTrkWMea[iMuPl] > 5) &&
      //     (Reco_mu_nPixWMea[iMuPl] > 0) &&
      //     (std::fabs(Reco_mu_dxy[iMuPl]) < 0.3f) &&
      //     (std::fabs(Reco_mu_dz[iMuPl]) < 20.f) &&
      //     passMuonTypePl && Reco_mu_highPurity[Reco_QQ_mupl_idx[iMuPl]]==true;
      // const bool mumiSoft =
      //     (Reco_mu_nTrkWMea[iMuMi] > 5) &&
      //     (Reco_mu_nPixWMea[iMuMi] > 0) &&
      //     (std::fabs(Reco_mu_dxy[iMuMi]) < 0.3f) &&
      //     (std::fabs(Reco_mu_dz[iMuMi]) < 20.f) &&
      //     passMuonTypeMi && Reco_mu_highPurity[Reco_QQ_mumi_idx[iMuMi]]==true;
      // if (!(muplSoft && mumiSoft)) continue;

      // vertex probability cut
      if (Reco_QQ_VtxProb[irqq] < 0.01f) continue;

      // init weights
      weight = 1.;
      // if(isMC) weight = findNcoll(Centrality) * Gen_weight; // need it for PbPb23
      Double_t tnp_weight = 1.0;
      Double_t tnp_trig_w_pl = -1.0;
      Double_t tnp_trig_w_mi = -1.0;

      recoQQsign[irqq] = Reco_QQ_sign[irqq];

      if (Reco_QQ_sign[irqq] != 0) continue;

      // ----- HX, CS transformation -----
      // HX
      TVector3 muPlus_HX = MuPlusVector_Helicity(JP, mupl);
			float cosHXVar = muPlus_HX.CosTheta();
			float phiHXVar = muPlus_HX.Phi() * 180 / TMath::Pi();

			if (cosHXVar < 0) {
				// if phi value is smaller than -pi, add 2pi
				if ((phiHXVar - 135) < -180)
					phiHXVar += 225;
				else
					phiHXVar -= 135;
			}

			else if (cosHXVar > 0) {
				// if phi value is smaller than -pi, add 2pi
				if ((phiHXVar - 45) < -180)
					phiHXVar += 315;
				else
					phiHXVar -= 45;
			}

      // CS
      TVector3 muPlus_CS = MuPlusVector_CollinsSoper(JP, mupl);
      float cosCSVar = muPlus_CS.CosTheta();
			float phiCSVar = muPlus_CS.Phi() * 180 / TMath::Pi();

			if (cosCSVar < 0) {
				// if phi value is smaller than -pi, add 2pi
				if ((phiCSVar - 135) < -180)
					cosCSVar += 225;
				else
					cosCSVar -= 135;
			}

			else if (cosCSVar > 0) {
				// if phi value is smaller than -pi, add 2pi
				if ((phiCSVar - 45) < -180)
					cosCSVar += 315;
				else
					cosCSVar -= 45;
			}

      // ----- fill per-candidate outputs -----
      phi[nDimu] = JP.Phi();
      phi1[nDimu] = mupl.Phi();
      phi2[nDimu] = mumi.Phi();

      if (isMC) TnPweight[nDimu] = tnp_weight;
      else TnPweight[nDimu] = 1.0;

      mass[nDimu] = JP.M();
      
      y[nDimu] = JP.Rapidity();

      cBin[nDimu] = cBin_;
      
      pt[nDimu] = JP.Pt();
      pt1[nDimu] = mupl.Pt();
      pt2[nDimu] = mumi.Pt();

      eta[nDimu] = JP.Eta();
      eta1[nDimu] = mupl.Eta();
      eta2[nDimu] = mumi.Eta();

      phi[nDimu] = JP.Phi();
      phi1[nDimu] = mupl.Phi();
      phi2[nDimu] = mumi.Phi();

      recoQQsign[nDimu] = Reco_QQ_sign[irqq];

      ctau3D[nDimu] = Reco_QQ_ctau3D[irqq];
      ctau3DErr[nDimu] = Reco_QQ_ctauErr3D[irqq];
      ctau3DRes[nDimu] = (std::fabs(ctau3DErr[nDimu]) > 0.f)
                             ? (ctau3D[nDimu] / ctau3DErr[nDimu])
                             : 0.f;
      cosHX[nDimu] = cosHXVar;
      phiHX[nDimu] = phiHXVar;
      cosCS[nDimu] = cosCSVar;
      phiCS[nDimu] = phiCSVar;

      ++nDimu;
    } // end of dimuon loop

    if (nDimu > 0)
    {
      ++count_dimuon;
      flowTree->Fill();
      if (VERBOSE_EVENT) cout << "  -> Fill (" << nDimu << " dimuons)\n";
    }
  }

  // ===== save results =====
  fFlowSkim->cd();
  flowTree->Write("myTree");
  fFlowSkim->Close();

  // ===== release memory =====
  delete JP_Reco;
  delete mupl_Reco;
  delete mumi_Reco;

  // ----- print -----
  cout << "Total saved dimuon: " << count_dimuon << "\n";

  cout << "Finish onia_to_skim_data()\n";
}

TVector3 MuPlusVector_Helicity(const TLorentzVector &QQLV_Lab, const TLorentzVector &MuPlusLV_Lab) {
	// ******** Transform variables of muons from the lab frame to the upsilon's rest frame ******** //
	TVector3 QQVector_Lab = QQLV_Lab.Vect();
	TLorentzVector MuPlusLV_QQRestFrame(MuPlusLV_Lab);

	//(Note. TLorentzVector.BoostVector() gives beta(=px/E,py/E,pz/E) of the parents)
	//(TLorentzVector.Boost() boosts from the rod frame to the lab frame, so plug in -beta to get lab to rod)
	MuPlusLV_QQRestFrame.Boost(-QQLV_Lab.BoostVector());

	// ******** Rotate the coordinates ******** //
	TVector3 MuPlusVec_Boosted = MuPlusLV_QQRestFrame.Vect();

	//Note: TVector3.Rotate() rotates the vectors, not the coordinates, so should rotate -phi and -theta

	MuPlusVec_Boosted.RotateZ(-QQVector_Lab.Phi());

	MuPlusVec_Boosted.RotateY(-QQVector_Lab.Theta());

	return MuPlusVec_Boosted;
}

// Lab to Collins-Soper
// requires the beam parameters
TVector3 MuPlusVector_CollinsSoper(const TLorentzVector &QQLV_Lab, const TLorentzVector &MuPlusLV_Lab) {
	// ******** Set beam energy for the Collins-Soper reference frame ******** //
	double sqrt_S_NN = 5.32;                    //(Center of mass Energy per nucleon pair in TeV)
	double beamEnergy = sqrt_S_NN * 1000. / 2.; //(in GeV) (Note. sqrt_S_NN = sqrt(2*E1*E2+2*p1*p2) = 2E1 when two beams have the same E)

	// ******** HX to CS (rotation from HX frame to CS frame) ******** //
	// (1. Boost two beams to upsilon's rest frame)
	// (2. Rotate the coordinates)
	// (3. Get angle between two beams(b1 and -b2), and between b1 and ZHX in the upsilon's rest frame)
	// (4. Calculate delta (angle btw ZHX and ZCS))

	// ******** Transform variables of beams from the lab frame to the upsilon's rest frame ******** //
	TLorentzVector Beam1LV_Boosted(0., 0., beamEnergy, beamEnergy);
	TLorentzVector Beam2LV_Boosted(0., 0., -beamEnergy, beamEnergy); // mind the sign!!

	Beam1LV_Boosted.Boost(-QQLV_Lab.BoostVector());
	Beam2LV_Boosted.Boost(-QQLV_Lab.BoostVector());

	// ******** Rotate the coordinates ******** //
	TVector3 Beam1Vector_QQRestFrame(Beam1LV_Boosted.Vect());
	TVector3 Beam2Vector_QQRestFrame(Beam2LV_Boosted.Vect());

	TVector3 QQVector_Lab = QQLV_Lab.Vect();
	Beam1Vector_QQRestFrame.RotateZ(-QQVector_Lab.Phi());
	Beam1Vector_QQRestFrame.RotateY(-QQVector_Lab.Theta());

	Beam2Vector_QQRestFrame.RotateZ(-QQVector_Lab.Phi());
	Beam2Vector_QQRestFrame.RotateY(-QQVector_Lab.Theta());

	// ******** Calculate the angle between z_HX and z_CS ******** //
	TVector3 ZHXunitVec(0, 0, 1.);                                                  //(define z_HX unit vector)
	double Angle_B1ZHX = Beam1Vector_QQRestFrame.Angle(ZHXunitVec);                //(angle between beam1 and z_HX)
	double Angle_B2ZHX = Beam2Vector_QQRestFrame.Angle(-ZHXunitVec);               //(angle between beam2 and -z_HX =(-beam2 and z_HX) )
	double Angle_B1miB2 = Beam1Vector_QQRestFrame.Angle(-Beam2Vector_QQRestFrame); //(angle between beam1 and -beam2)

	double delta = 0; //(define and initialize the angle between z_HX and z_CS)

	// Maths for calculating the angle between z_HX and z_CS is different depending on the sign of the beam1's z-coordinate)
	if (Angle_B1ZHX > Angle_B2ZHX)
		delta = Angle_B2ZHX + Angle_B1miB2 / 2.;
	else if (Angle_B1ZHX < Angle_B2ZHX)
		delta = Angle_B1ZHX + Angle_B1miB2 / 2.;
	else
		std::cout << "beam1PvecBoosted.Pz() = 0?" << std::endl;

	// ******** Rotate the coordinates along the y-axis by the angle between z_HX and z_CS ******** //
	TVector3 MuPlusVec_CS(MuPlusVector_Helicity(QQLV_Lab, MuPlusLV_Lab));

	MuPlusVec_CS.RotateY(delta);

	return MuPlusVec_CS;
}