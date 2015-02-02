#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLatex.h>
#include "Math/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

void plotDileptonMass(const TString infile1, const TString infile2) {

  TFile* inputFile = new TFile();
  TTree* inputTree = new TTree();

  //
  // Declare variables to read in ntuple
  //
  LorentzVector *dilep=0;

  //
  // Declare histograms
  //
  //
  TH1D *hmass_pu4bx50  = new TH1D("hmass_pu4bx50","",60,60,120);
        hmass_pu4bx50->SetStats(0);
        hmass_pu4bx50->SetLineColor(1);
        hmass_pu4bx50->GetXaxis()->SetTitle("Dilepton Mass  [ GeV ]");
        hmass_pu4bx50->GetXaxis()->SetTitleSize(0.05);
        hmass_pu4bx50->GetXaxis()->SetTitleOffset(0.75);
        hmass_pu4bx50->GetYaxis()->SetTitle("Events / 1.0 GeV");
        hmass_pu4bx50->GetYaxis()->SetTitleSize(0.055);
        hmass_pu4bx50->GetYaxis()->SetTitleOffset(0.9);
  TH1D *hmass_pu20bx25  = new TH1D("hmass_pu20bx25","",60,60,120);
        hmass_pu20bx25->SetStats(0);
        hmass_pu20bx25->SetLineColor(2);
        hmass_pu20bx25->GetXaxis()->SetTitle("Dilepton Mass  [ GeV ]");
        hmass_pu20bx25->GetXaxis()->SetTitleSize(0.05);
        hmass_pu20bx25->GetXaxis()->SetTitleOffset(0.75);
        hmass_pu20bx25->GetYaxis()->SetTitle("Events / 1.0 GeV");
        hmass_pu20bx25->GetYaxis()->SetTitleSize(0.055);
        hmass_pu20bx25->GetYaxis()->SetTitleOffset(0.9);

  for(int jfile=0;jfile<2;jfile++) {

    //
    // Setup input ntuple
    //
    if(jfile==0) inputFile = new TFile(infile1);
    else if(jfile==1) inputFile = new TFile(infile2);
    inputTree = (TTree*)inputFile->Get("Events");

    inputTree->SetBranchAddress("dilep", &dilep);   // Dilepton 4-vector
           
    for(int jentry=0;jentry<inputTree->GetEntries();jentry++) {
    //for(int jentry=0;jentry<10000;jentry++) {
      inputTree->GetEntry(jentry);

      //
      // Fill histograms
      //
      if(jfile==0) hmass_pu4bx50->Fill(dilep->M());
      else if (jfile==1) hmass_pu20bx25->Fill(dilep->M());
    }

  }

  TLegend *leg = new TLegend(0.5373563,0.7097458,0.8577586,0.8474576,NULL,"brNDC");
  leg->SetBorderSize(1);
  leg->SetTextFont(62);
  leg->SetTextSize(0.05084746);
  leg->SetLineColor(0);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->AddEntry(hmass_pu4bx50,"PU 4 BX 50","lp");
  leg->AddEntry(hmass_pu20bx25,"PU 20 BX 25","lp");

  hmass_pu4bx50->Scale(1/hmass_pu4bx50->Integral());
  hmass_pu20bx25->Scale(1/hmass_pu20bx25->Integral());

  hmass_pu4bx50->Scale(1953.0*0.1);
  hmass_pu20bx25->Scale(1953.0*0.1);

  //
  // Save plots
  //
  TCanvas *cmass = new TCanvas("cmass", "cmass",18,40,700,500);
  cmass->Range(18.78366,0.5855855,56.2702,0.9308669);
  cmass->SetFillColor(0);
  cmass->SetBorderMode(0);
  cmass->SetBorderSize(2);
  cmass->SetFrameBorderMode(0);
  cmass->cd();
  hmass_pu4bx50->Draw();
  hmass_pu20bx25->Draw("same");
  leg->Draw("same");
  cmass->Print("dileptonmass.png");
}
