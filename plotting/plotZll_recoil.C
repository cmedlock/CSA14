#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TVector2.h>
#include <TF1.h>
#include <TMath.h>
#include <TDirectory.h>
#include "Math/GenVector/LorentzVector.h"

// All of this code is a reimplementation of Kevin Sung's code from the 8 TeV analysis.

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

// Convert int to binary value
bool Convert(unsigned int val,bool print=kFALSE)
{
   unsigned int mask = 1 << (sizeof(int) * 8 - 1);
   bool lastDigit;
   for(unsigned int i = 0; i < sizeof(int) * 8; i++)
   {
      if( (val & mask) == 0 ) {
        if(print) cout << "0";
        lastDigit=0;
      } else {
        if(print) cout << "1";
        lastDigit=1;
      }
      mask  >>= 1;
   }
   if(print) cout << endl;
   return lastDigit;
}

void plotZll_recoil(const TString inputFileName = "/scratch/cmedlock/PHYS14/DYJetsToLL_M-50_13TeV-madgraph-pythia8_PU20bx25_PHYS14_25_V1-v1_00000.root") {

  //
  // Setup input ntuple
  //
  TFile* inputFile = new TFile(inputFileName);
  TTree* inputTree = (TTree*)inputFile->Get("Events");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t MASS_LOW  = 40;
  const Double_t MASS_HIGH = 200;
  
  Int_t NVTXBINS = 70; // 70 for Wenu_p, 35 for Wmunu_p

  //
  // Declare variables to read in ntuple
  //
  Int_t nVtx, nEvents;
  TVector2 *vtype1pfMET=0, *vrawpfMET=0, *vgenMET=0;
  LorentzVector *lep1=0, *lep2=0, *dilep=0;

  inputTree->SetBranchAddress("nVtx",         &nVtx);        // number of vertices
  inputTree->SetBranchAddress("vtype1pfMET",  &vtype1pfMET); // type-1 corrected pf MET
  inputTree->SetBranchAddress("vrawpfMET",    &vrawpfMET);   // raw pf MET
  inputTree->SetBranchAddress("vgenMET",      &vgenMET);     // generated MET
  inputTree->SetBranchAddress("nEvents",      &nEvents);
  inputTree->SetBranchAddress("lep1",         &lep1);
  inputTree->SetBranchAddress("lep2",         &lep2);
  inputTree->SetBranchAddress("dilep",        &dilep);
           
  //
  // Declare histograms
  //
  //
  // MET, phi(MET)
  TH1D *htype1corr  = new TH1D("htype1corr","",100,0,150);
        htype1corr->SetStats(0); htype1corr->SetLineColor(4);
  TH1D *htype1phicorr = new TH1D("htype1phicorr","",100,-3.5,3.5);
        htype1phicorr->SetStats(0); htype1phicorr->SetLineColor(4);
  TH1D *hgen = new TH1D("hgen","MET",100,0,150);
        hgen->GetYaxis()->SetTitle("Events / 1.5 GeV"); hgen->GetYaxis()->SetTitleOffset(1.5);
        hgen->GetXaxis()->SetTitle("MET [GeV]"); hgen->GetXaxis()->SetTitleOffset(1.2);
        hgen->SetStats(0); hgen->SetLineColor(2);
  TH1D *hgenphi = new TH1D("hgenphi","",100,-3.5,3.5);
        hgenphi->SetStats(0); hgenphi->SetLineColor(2);
  TH2D *hmetx = new TH2D("hmetx","MET_{x} v. Number of Vertices",NVTXBINS,0,NVTXBINS,100,-150,150);
        hmetx->GetYaxis()->SetTitle("MET_{x} [GeV]");
        hmetx->GetXaxis()->SetTitle("Number of vertices");
  TH2D *hmety = new TH2D("hmety","MET_{y} v. Number of Vertices",NVTXBINS,0,NVTXBINS,100,-150,150);
        hmety->GetYaxis()->SetTitle("MET_{y} [GeV]");
        hmety->GetXaxis()->SetTitle("Number of vertices");
  // Dilepton mass and pT
  TH1D *hzmass = new TH1D("hzmass","Dilepton Mass",60,60,120);
        hzmass->GetYaxis()->SetTitle("Events / 1 GeV"); hzmass->GetYaxis()->SetTitleOffset(1.5);
        hzmass->GetXaxis()->SetTitle("M_{\el\el} [GeV]");
  TH1D *hzpt = new TH1D("hzpt","Dilepton pT",100,0,200);
        hzpt->GetYaxis()->SetTitle("Events / 4 GeV"); hzpt->GetYaxis()->SetTitleOffset(1.5);
        hzpt->GetXaxis()->SetTitle("pT_{e} [GeV]"); hzpt->GetXaxis()->SetTitleOffset(1.2);
  // Hadronic recoil
  TH2D *hescale = new TH2D("hescale","MET Energy Scale",50,40,200,100,0,1.4);
        hescale->GetYaxis()->SetTitle("-<u_{||}> / pT_{\el\el}"); hescale->GetYaxis()->SetTitleOffset(1.5);
        hescale->GetXaxis()->SetTitle("pT_{\el\el}"); hescale->GetXaxis()->SetTitleOffset(1.2);
  TH2D *hres_upar = new TH2D("hres_upar","Resolution of u_{||}",50,40,200,50,0,100);
        hres_upar->GetYaxis()->SetTitle("#sigma(u_{||}) [GeV]"); hres_upar->GetYaxis()->SetTitleOffset(1.5);
        hres_upar->GetXaxis()->SetTitle("pT_{\el\el} [GeV]"); hres_upar->GetXaxis()->SetTitleOffset(1.2);
  TH2D *hres_uperp = new TH2D("hres_uperp","Resolution of u_{#perp}",50,40,200,50,0,100);
        hres_uperp->GetYaxis()->SetTitle("#sigma(u_{#perp}) [GeV]"); hres_uperp->GetYaxis()->SetTitleOffset(1.5);
        hres_uperp->GetXaxis()->SetTitle("pT_{\el\el} [GeV]"); hres_uperp->GetXaxis()->SetTitleOffset(1.2);

  //
  // Make legend
  //
  TLegend *leg = new TLegend(0.4181034,0.6758475,0.6954023,0.8135593,NULL,"brNDC");
  leg->SetTextFont(62);
  leg->SetTextSize(0.03330866);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetBorderSize(0);

  leg->AddEntry(htype1corr,"Type 1 + xy Shift Corrected PF","l");
  leg->AddEntry(hgen,"Generated","l");

  TLegend *legphi = new TLegend(0.420977,0.3114407,0.6982759,0.4491525,NULL,"brNDC");
  legphi->SetTextFont(62);
  legphi->SetTextSize(0.03330866);
  legphi->SetLineColor(1);
  legphi->SetLineStyle(1);
  legphi->SetLineWidth(1);
  legphi->SetFillColor(0);
  legphi->SetFillStyle(1001);
  legphi->SetBorderSize(0);

  legphi->AddEntry(htype1corr,"Type 1 + xy Shift Corrected PF","l");
  legphi->AddEntry(hgen,"Generated","l");

  Int_t totalEvents=0, nsel=0;

  for(int kentry=0;kentry<inputTree->GetEntries();kentry++) {
    inputTree->GetEntry(kentry);
    totalEvents += nEvents;

    //
    // Fill histograms
    //
    hmetx->Fill(nVtx,vtype1pfMET->Px());
    hmety->Fill(nVtx,vtype1pfMET->Py());
  }

  // Loop through nVtx bins and find the mean value of metx and mety in each bin
  // Plot nVtx v. mean values of metx/mety in a separate histogram (1 for metx and 1 for mety)
  TH1D* hmetx_proj = new TH1D();
  TH1D* hmety_proj = new TH1D();
  Double_t meanmetx, meanmety;
  TH2D* hmetxfit = new TH2D("hmetxfit","MET_{x} v. Number of vertices",NVTXBINS,0,NVTXBINS,100,-25,5);
        hmetxfit->GetXaxis()->SetTitle("Number of vertices");
        hmetxfit->GetYaxis()->SetTitle("<MET_{x}> [GeV]");
        hmetxfit->SetMarkerSize(21);
  TH2D* hmetyfit = new TH2D("hmetyfit","MET_{y} v. Number of vertices",NVTXBINS,0,NVTXBINS,100,-25,5);
        hmetyfit->GetXaxis()->SetTitle("Number of vertices");
        hmetyfit->GetYaxis()->SetTitle("<MET_{y}> [GeV]");
        hmetyfit->SetMarkerSize(21);
  for(int jbin=1;jbin<hmetx->GetNbinsX()+1;jbin++) {
    hmetx_proj = hmetx->ProjectionY("metx_proj",jbin,jbin+1,"");
    hmety_proj = hmety->ProjectionY("mety_proj",jbin,jbin+1,"");
    meanmetx = hmetx_proj->GetMean();
    meanmety = hmety_proj->GetMean();
    hmetxfit->Fill(jbin,meanmetx);
    hmetyfit->Fill(jbin,meanmety);
  }

  // Fit each of these with a line to calculate the appropriate correction factor for an event with the given number of vertices
  // Estimation of initial values for parameters can be streamlined (don't hardcode in the 49 and 50)
  TF1 *flinex = new TF1("flinex","[0]+[1]*x",0,NVTXBINS);
  Double_t xpar0=0, xpar1=0;
  xpar0 = hmetxfit->ProjectionY("",2,3,"")->GetMean(); // mean of metx for events with 1 vertex (minimum)
  xpar1 = (hmetxfit->ProjectionY("",49,50,"")->GetMean()-xpar0)/hmetxfit->GetNbinsX(); // slope estimate from 2 outermost points
  flinex->SetParameter(0,xpar0);
  flinex->SetParameter(1,xpar1);
  hmetxfit->Fit(flinex);
  TF1 *fliney = new TF1("fliney","[0]+[1]*x",0,NVTXBINS);
  Double_t ypar0=0, ypar1=0;
  ypar0 = hmetyfit->ProjectionY("",2,3,"")->GetMean(); // mean of mety for events with 1 vertex (minimum)
  ypar1 = (hmetyfit->ProjectionY("",NVTXBINS-1,NVTXBINS,"")->GetMean()-ypar0)/NVTXBINS; // slope estimate from 2 outermost points
  fliney->SetParameter(0,ypar0);
  fliney->SetParameter(1,ypar1);
  hmetyfit->Fit(fliney);

  // Correct the MET
  Double_t flinex0=flinex->GetParameter(0), flinex1=flinex->GetParameter(1);
  Double_t fliney0=fliney->GetParameter(0), fliney1=fliney->GetParameter(1);
  Double_t type1x=0, type1y=0;
  Double_t type1phi=0, rawphi=0, genphi=0;
  Double_t corrMETx=0, corrMETy=0, corrMET=0, corrMETphi=0;

  for(int jentry=0;jentry<inputTree->GetEntries();jentry++) {
    inputTree->GetEntry(jentry);

    // Mass window
    if((dilep->M()<MASS_LOW) || (dilep->M()>MASS_HIGH)) continue;
    // For hadronic recoil study, only look at high pT Z's
    if(dilep->Pt() < 50) continue;
 
    TVector2 vDilepPt(dilep->Px(),dilep->Py());

    nsel++;

    type1x = vtype1pfMET->Px();
    type1y = vtype1pfMET->Py();

    type1phi = vtype1pfMET->Phi();
    if(type1phi>TMath::Pi()) type1phi -= 2*TMath::Pi();
    rawphi = vrawpfMET->Phi();
    if(rawphi>TMath::Pi())   rawphi   -= 2*TMath::Pi();
    genphi = vgenMET->Phi();
    if(genphi>TMath::Pi())   genphi   -= 2*TMath::Pi();

    corrMETx = type1x-flinex0-flinex1*nVtx;
    corrMETy = type1y-fliney0-fliney1*nVtx;
    corrMET  = TMath::Sqrt(corrMETx*corrMETx + corrMETy*corrMETy);
    corrMETphi = TMath::ATan2(corrMETy,corrMETx);
    TVector2 vcorrMET(corrMETx,corrMETy);

    // Calculate hadronic recoil
    TVector2 uT    = -1*(vDilepPt+vcorrMET);
    TVector2 uPar  = uT.Proj(vDilepPt);
    TVector2 uPerp = uT-uPar;

    //
    // Fill histograms
    //
    hzmass->Fill(dilep->M());
    hzpt->Fill(dilep->Pt());
    hescale->Fill(vDilepPt.Mod(),uPar.Mod()/vDilepPt.Mod());
    hres_upar->Fill(vDilepPt.Mod(),uPar.Mod());
    hres_uperp->Fill(vDilepPt.Mod(),uPerp.Mod());

    htype1corr->Fill(corrMET);
    htype1phicorr->Fill(corrMETphi);
    hgen->Fill(vgenMET->Mod());
    hgenphi->Fill(genphi);

  }

  hescale->FitSlicesY();
  TH1D *hescale_means = (TH1D*)gDirectory->Get("hescale_1");
        hescale_means->SetTitle("MET Energy Scale");
        hescale_means->GetYaxis()->SetTitle("-<u_{||}> / pT_{\el\el}"); hescale_means->GetYaxis()->SetTitleOffset(1.5);
        hescale_means->GetXaxis()->SetTitle("pT_{\el\el} [GeV]"); hescale_means->GetXaxis()->SetTitleOffset(1.2);
  hres_upar->FitSlicesY();
  TH1D *hres_upar_rms = (TH1D*)gDirectory->Get("hres_upar_2");
        hres_upar_rms->SetTitle("Resolution of u_{||}");
        hres_upar_rms->GetYaxis()->SetTitle("#sigma(u_{||}) [GeV]"); hres_upar_rms->GetYaxis()->SetTitleOffset(1.5);
        hres_upar_rms->GetXaxis()->SetTitle("pT_{\el\el} [GeV]"); hres_upar_rms->GetXaxis()->SetTitleOffset(1.2);
  hres_uperp->FitSlicesY();
  TH1D *hres_uperp_rms = (TH1D*)gDirectory->Get("hres_uperp_2");
        hres_uperp_rms->SetTitle("Resolution of u_{#perp}");
        hres_uperp_rms->GetYaxis()->SetTitle("#sigma(u_{#perp}) [GeV]"); hres_uperp_rms->GetYaxis()->SetTitleOffset(1.5);
        hres_uperp_rms->GetXaxis()->SetTitle("pT_{\el\el} [GeV]"); hres_uperp_rms->GetXaxis()->SetTitleOffset(1.2);

  //
  // Save plots
  //

  TCanvas* czmass = new TCanvas();
  czmass->cd(); hzmass->Draw(); czmass->Print("czmass.png"); czmass->Close();
  TCanvas* czpt = new TCanvas();
  czpt->cd(); hzpt->Draw(); czpt->Print("czpt.png"); czpt->Close();
  TCanvas* cescale = new TCanvas();
  cescale->cd(); hescale_means->Draw(); cescale->Print("cescale.png"); // cescale->Close();
  TCanvas* cres_upar = new TCanvas();
  cres_upar->cd(); hres_upar_rms->Draw(); cres_upar->Print("cres_upar.png"); // cres_upar->Close();
  TCanvas* cres_uperp = new TCanvas();
  cres_uperp->cd(); hres_uperp_rms->Draw(); cres_uperp->Print("cres_uperp.png"); // cres_uperp->Close();

}
