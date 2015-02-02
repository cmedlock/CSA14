#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TMath.h>

void plotMET(const TString inputFileName, const TString outputFileName) {

  //
  // Setup input ntuple
  //
  TFile* inputFile = new TFile(inputFileName);
  TTree* inputTree = (TTree*)inputFile->Get("Events");

  //
  // Declare variables to read in ntuple
  //
  Float_t type1pfmet, type1pfmetPhi;

  inputTree->SetBranchAddress("type1pfmet",    &type1pfmet);   // Type-1 corrected PF MET
  inputTree->SetBranchAddress("type1pfmetPhi", &type1pfmetPhi);// Type-1 corrected PF MET phi

  //
  // Declare histograms
  //
  //
  TH1D *htype1  = new TH1D("htype1","Type-1 Corrected PF MET",100,0,150);
        htype1->SetStats(0);
        htype1->SetLineColor(1);
        htype1->GetXaxis()->SetTitle("MET  [ GeV ]");
        htype1->GetXaxis()->SetTitleSize(0.05);
        htype1->GetXaxis()->SetTitleOffset(0.75);
        htype1->GetYaxis()->SetTitle("Events / 1.5 GeV");
        htype1->GetYaxis()->SetTitleSize(0.055);
        htype1->GetYaxis()->SetTitleOffset(0.95);
  TH1D *htype1phi = new TH1D("htype1phi","Type-1 Corrected PF MET Phi",100,-3.5,3.5);
        htype1phi->SetStats(0);
        htype1phi->SetLineColor(1);
        htype1phi->GetXaxis()->SetTitle("#phi(MET)  [ GeV ]");
        htype1phi->GetXaxis()->SetTitleSize(0.05);
        htype1phi->GetXaxis()->SetTitleOffset(0.75);
        htype1phi->GetYaxis()->SetTitle("Events / 0.07");
        htype1phi->GetYaxis()->SetTitleSize(0.055);
        htype1phi->GetYaxis()->SetTitleOffset(0.8);

  for(int jentry=0;jentry<inputTree->GetEntries();jentry++) {
    inputTree->GetEntry(jentry);

    if(type1pfmetPhi>TMath::Pi()) type1pfmetPhi -= 2*TMath::Pi();

    //
    // Fill histograms
    //
    htype1->Fill(type1pfmet);
    htype1phi->Fill(type1pfmetPhi);
  }

  //
  // Normalize
  //
  htype1->Scale(1/htype1->Integral());
  htype1->Scale(19551.0*1.0/9.0); // Cross section [pb] * branching ratio

  htype1phi->Scale(1/htype1phi->Integral());
  htype1phi->Scale(19551.0*1.0/9.0);  // Cross section [pb] * branching ratio

  //
  // Save plots
  //
  TFile* outputFile = new TFile(outputFileName,"RECREATE");
  htype1->Write();
  htype1phi->Write();
  outputFile->Close();

}
