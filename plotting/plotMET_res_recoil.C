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

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

void plotMET_res_recoil(const TString inputFileName, const TString outputFileName) {

  //
  // Setup input ntuple
  //
  TFile* inputFile = new TFile(inputFileName);
  TTree* inputTree = (TTree*)inputFile->Get("Events");
  
  //
  // Declare variables to read in ntuple
  //
  Float_t type1pfmet, type1pfmetPhi, genmet, genmetPhi;
  Float_t u1, u2;
  LorentzVector *dilep=0;

  inputTree->SetBranchAddress("type1pfmet",    &type1pfmet);    // Type-1 corrected PF MET
  inputTree->SetBranchAddress("type1pfmetPhi", &type1pfmetPhi); // Type-1 corrected PF MET phi
  inputTree->SetBranchAddress("genmet",        &genmet);        // Generator level MET
  inputTree->SetBranchAddress("genmetPhi",     &genmetPhi);     // Generator level MET phi
  inputTree->SetBranchAddress("u1",            &u1);            // parallel component of recoil
  inputTree->SetBranchAddress("u2",            &u2);            // perpendicular component of recoil
  inputTree->SetBranchAddress("dilep",         &dilep);         // Dilepton 4-vector
           
  //
  // Declare histograms
  //
  //
  // MET
  TH1D *htype1res = new TH1D("htype1res","",50,-1,10);
        htype1res->SetStats(0);
        htype1res->SetLineColor(1);
        htype1res->SetTitle("Resolution of Type-1 Corrected PF MET");
        htype1res->GetXaxis()->SetTitle("( MET - Gen. MET ) / Gen. MET");
        htype1res->GetYaxis()->SetTitle("Events / 0.22 GeV");
  // Hadronic recoil
  TH2D *hres_upar = new TH2D("hres_upar","Resolution of u_{||}",50,50,100,100,0,100);
        hres_upar->GetYaxis()->SetTitle("#sigma(u_{||}) [GeV]");
        hres_upar->GetXaxis()->SetTitle("pT_{#mu#mu} [GeV]");
  TH2D *hres_uperp = new TH2D("hres_uperp","Resolution of u_{#perp}",50,50,100,50,0,100);
        hres_uperp->GetYaxis()->SetTitle("#sigma(u_{#perp}) [GeV]");
        hres_uperp->GetXaxis()->SetTitle("pT_{#mu#mu} [GeV]");


  for(int jentry=0;jentry<inputTree->GetEntries();jentry++) {
    inputTree->GetEntry(jentry);

    // Mass window
    if((dilep->M()<40) || (dilep->M()>200)) continue;
    // For hadronic recoil study, only look at high pT Z's
    if(dilep->Pt() < 50) continue;

    if(type1pfmetPhi>TMath::Pi()) type1pfmetPhi -= 2*TMath::Pi();
    if(genmetPhi>TMath::Pi())   genmetPhi -= 2*TMath::Pi();

    //
    // Fill histograms
    //
    hres_upar->Fill(dilep->Pt(),fabs(u1));
    hres_uperp->Fill(dilep->Pt(),fabs(u2));

    htype1res->Fill((type1pfmet-genmet)/genmet);
  }

  hres_upar->FitSlicesY();
  TH1D *hres_upar_rms = (TH1D*)gDirectory->Get("hres_upar_2");
        hres_upar_rms->SetStats(0);
        hres_upar_rms->SetTitle("Resolution of u_{||}");
        hres_upar_rms->GetYaxis()->SetTitle("#sigma(u_{||}) [GeV]"); hres_upar_rms->GetYaxis()->SetTitleOffset(1.5);
        hres_upar_rms->GetXaxis()->SetRangeUser(0,100);
        hres_upar_rms->GetXaxis()->SetTitle("pT_{#mu#mu} [GeV]"); hres_upar_rms->GetXaxis()->SetTitleOffset(1.2);
  hres_uperp->FitSlicesY();
  TH1D *hres_uperp_rms = (TH1D*)gDirectory->Get("hres_uperp_2");
        hres_uperp_rms->SetStats(0);
        hres_uperp_rms->SetTitle("Resolution of u_{#perp}");
        hres_uperp_rms->GetYaxis()->SetTitle("#sigma(u_{#perp}) [GeV]"); hres_uperp_rms->GetYaxis()->SetTitleOffset(1.5);
        hres_uperp_rms->GetXaxis()->SetRangeUser(0,100);
        hres_uperp_rms->GetXaxis()->SetTitle("pT_{#mu#mu} [GeV]"); hres_uperp_rms->GetXaxis()->SetTitleOffset(1.2);

  //
  // Save plots
  //
  TFile* outputFile = new TFile(outputFileName,"RECREATE");
  hres_upar_rms->Write();
  hres_uperp_rms->Write();
  htype1res->Write();
  outputFile->Close();
}
