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

void mergeMetPlots(const TString inputFile_PU4bx50, const TString inputFile_PU20bx25) {

  //
  // Setup input ntuple
  //
  TFile* file_pu4bx50  = new TFile(inputFile_PU4bx50);
  TFile* file_pu20bx25 = new TFile(inputFile_PU20bx25);

  TH1D* htype1_pu4bx50 = (TH1D*)file_pu4bx50->Get("htype1");
        htype1_pu4bx50->SetStats(0);
        htype1_pu4bx50->SetLineColor(1);
        htype1_pu4bx50->GetXaxis()->SetTitleSize(0.05);
        htype1_pu4bx50->GetXaxis()->SetTitleOffset(0.75);
        htype1_pu4bx50->GetYaxis()->SetTitleSize(0.055);
        htype1_pu4bx50->GetYaxis()->SetTitleOffset(0.95);
  TH1D* htype1_pu20bx25 = (TH1D*)file_pu20bx25->Get("htype1");
        htype1_pu20bx25->SetStats(0);
        htype1_pu20bx25->SetLineColor(2);
        htype1_pu20bx25->GetXaxis()->SetTitleSize(0.05);
        htype1_pu20bx25->GetXaxis()->SetTitleOffset(0.75);
        htype1_pu20bx25->GetYaxis()->SetTitleSize(0.055);
        htype1_pu20bx25->GetYaxis()->SetTitleOffset(0.95);
  TH1D* htype1phi_pu4bx50 = (TH1D*)file_pu4bx50->Get("htype1phi");
        htype1phi_pu4bx50->SetStats(0);
        htype1phi_pu4bx50->SetLineColor(1);
        htype1phi_pu4bx50->GetXaxis()->SetTitleSize(0.05);
        htype1phi_pu4bx50->GetXaxis()->SetTitleOffset(0.75);
        htype1phi_pu4bx50->GetYaxis()->SetTitleSize(0.055);
        htype1phi_pu4bx50->GetYaxis()->SetTitleOffset(0.95);
  TH1D* htype1phi_pu20bx25 = (TH1D*)file_pu20bx25->Get("htype1phi");
        htype1phi_pu20bx25->SetStats(0);
        htype1phi_pu20bx25->SetLineColor(2);
        htype1phi_pu20bx25->GetXaxis()->SetTitleSize(0.05);
        htype1phi_pu20bx25->GetXaxis()->SetTitleOffset(0.75);
        htype1phi_pu20bx25->GetYaxis()->SetTitleSize(0.055);
        htype1phi_pu20bx25->GetYaxis()->SetTitleOffset(0.95);

  TLegend *leg = new TLegend(0.5373563,0.7097458,0.8577586,0.8474576,NULL,"brNDC");
  leg->SetBorderSize(1);
  leg->SetTextFont(62);
  leg->SetTextSize(0.05084746);
  leg->SetLineColor(0);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->AddEntry(htype1_pu4bx50,"PU 4 BX 50","lp");
  leg->AddEntry(htype1_pu20bx25,"PU 20 BX 25","lp");

  TLegend *legphi = new TLegend(0.5258621,0.2944915,0.8146552,0.4322034,NULL,"brNDC");
  legphi->SetBorderSize(1);
  legphi->SetTextFont(62);
  legphi->SetTextSize(0.05084746);
  legphi->SetLineColor(0);
  legphi->SetLineStyle(1);
  legphi->SetLineWidth(1);
  legphi->SetFillColor(0);
  legphi->SetFillStyle(1001);
  legphi->AddEntry(htype1_pu4bx50,"PU 4 BX 50","lp");
  legphi->AddEntry(htype1_pu20bx25,"PU 20 BX 25","lp");

  //
  // Normalize
  //
  htype1_pu4bx50->Scale(1/htype1_pu4bx50->Integral());
  htype1_pu20bx25->Scale(1/htype1_pu20bx25->Integral());
  htype1phi_pu4bx50->Scale(1/htype1phi_pu4bx50->Integral());
  htype1phi_pu20bx25->Scale(1/htype1phi_pu20bx25->Integral());

  //
  // Save plots
  //
  TCanvas *ctype1 = new TCanvas("ctype1", "ctype1",18,40,700,500);
  ctype1->Range(18.78366,0.5855855,56.2702,0.9308669);
  ctype1->SetFillColor(0);
  ctype1->SetBorderMode(0);
  ctype1->SetBorderSize(2);
  ctype1->SetFrameBorderMode(0);
  ctype1->SetFrameBorderMode(0);
  ctype1->cd();
  htype1_pu4bx50->Draw("");
  htype1_pu20bx25->Draw("same");
  leg->Draw("same");
  ctype1->Print("ctype1pfmet.png");
  TCanvas *ctype1phi = new TCanvas("ctype1phi", "ctype1phi",18,40,700,500);
  ctype1phi->Range(18.78366,0.5855855,56.2702,0.9308669);
  ctype1phi->SetFillColor(0);
  ctype1phi->SetBorderMode(0);
  ctype1phi->SetBorderSize(2);
  ctype1phi->SetFrameBorderMode(0);
  ctype1phi->SetFrameBorderMode(0);
  ctype1phi->cd();
  htype1phi_pu20bx25->Draw("");
  htype1phi_pu4bx50->Draw("same");
  legphi->Draw("same");
  ctype1phi->Print("ctype1phi.png");

}
