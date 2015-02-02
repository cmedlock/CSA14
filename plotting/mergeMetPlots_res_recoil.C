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

void mergeMetPlots_res_recoil(const TString inputFile_PU4bx50, const TString inputFile_PU20bx25) {

  //
  // Setup input ntuple
  //
  TFile* file_pu4bx50  = new TFile(inputFile_PU4bx50);
  TFile* file_pu20bx25 = new TFile(inputFile_PU20bx25);

  TH1D* hres_upar_rms_pu4bx50 = (TH1D*)file_pu4bx50->Get("hres_upar_2");
        hres_upar_rms_pu4bx50->SetStats(0);
        hres_upar_rms_pu4bx50->SetLineColor(1);
        hres_upar_rms_pu4bx50->GetXaxis()->SetTitleSize(0.05);
        hres_upar_rms_pu4bx50->GetXaxis()->SetTitleOffset(0.9);
        hres_upar_rms_pu4bx50->GetYaxis()->SetTitleSize(0.065);
        hres_upar_rms_pu4bx50->GetYaxis()->SetTitleOffset(0.6);
  TH1D* hres_uperp_rms_pu4bx50 = (TH1D*)file_pu4bx50->Get("hres_uperp_2");
        hres_uperp_rms_pu4bx50->SetStats(0);
        hres_uperp_rms_pu4bx50->SetLineColor(1);
        hres_uperp_rms_pu4bx50->GetXaxis()->SetTitleSize(0.05);
        hres_uperp_rms_pu4bx50->GetXaxis()->SetTitleOffset(0.9);
        hres_uperp_rms_pu4bx50->GetYaxis()->SetTitleSize(0.065);
        hres_uperp_rms_pu4bx50->GetYaxis()->SetTitleOffset(0.6);
  TH1D* htype1res_pu4bx50 = (TH1D*)file_pu4bx50->Get("htype1res");
        htype1res_pu4bx50->SetStats(0);
        htype1res_pu4bx50->SetLineColor(1);
        htype1res_pu4bx50->GetXaxis()->SetTitleSize(0.05);
        htype1res_pu4bx50->GetXaxis()->SetTitleOffset(0.75);
        htype1res_pu4bx50->GetYaxis()->SetTitleSize(0.065);
        htype1res_pu4bx50->GetYaxis()->SetTitleOffset(0.8);

  TH1D* hres_upar_rms_pu20bx25 = (TH1D*)file_pu20bx25->Get("hres_upar_2");
        hres_upar_rms_pu20bx25->SetStats(0);
        hres_upar_rms_pu20bx25->SetLineColor(2);
        hres_upar_rms_pu20bx25->GetXaxis()->SetTitleSize(0.05);
        hres_upar_rms_pu20bx25->GetXaxis()->SetTitleOffset(0.9);
        hres_upar_rms_pu20bx25->GetYaxis()->SetTitleSize(0.065);
        hres_upar_rms_pu20bx25->GetYaxis()->SetTitleOffset(0.6);
  TH1D* hres_uperp_rms_pu20bx25 = (TH1D*)file_pu20bx25->Get("hres_uperp_2");
        hres_uperp_rms_pu20bx25->SetStats(0);
        hres_uperp_rms_pu20bx25->SetLineColor(2);
        hres_uperp_rms_pu20bx25->GetXaxis()->SetTitleSize(0.05);
        hres_uperp_rms_pu20bx25->GetXaxis()->SetTitleOffset(0.9);
        hres_uperp_rms_pu20bx25->GetYaxis()->SetTitleSize(0.065);
        hres_uperp_rms_pu20bx25->GetYaxis()->SetTitleOffset(0.6);
  TH1D* htype1res_pu20bx25 = (TH1D*)file_pu20bx25->Get("htype1res");
        htype1res_pu20bx25->SetStats(0);
        htype1res_pu20bx25->SetLineColor(2);
        htype1res_pu20bx25->GetXaxis()->SetTitleSize(0.05);
        htype1res_pu20bx25->GetXaxis()->SetTitleOffset(0.75);
        htype1res_pu20bx25->GetYaxis()->SetTitleSize(0.065);
        htype1res_pu20bx25->GetYaxis()->SetTitleOffset(0.8);

  TLegend *leg = new TLegend(0.5431034,0.3347458,0.8635057,0.4724576,NULL,"brNDC");
  leg->SetBorderSize(1);
  leg->SetTextFont(62);
  leg->SetTextSize(0.05084746);
  leg->SetLineColor(0);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->AddEntry(htype1res_pu4bx50,"PU 4 BX 50","lp");
  leg->AddEntry(htype1res_pu20bx25,"PU 20 BX 25","lp");

  TLegend *legu = new TLegend(0.1408046,0.7055085,0.4612069,0.8432203,NULL,"brNDC");
  legu->SetBorderSize(1);
  legu->SetTextFont(62);
  legu->SetTextSize(0.05084746);
  legu->SetLineColor(0);
  legu->SetLineStyle(1);
  legu->SetLineWidth(1);
  legu->SetFillColor(0);
  legu->SetFillStyle(1001);
  legu->AddEntry(htype1res_pu4bx50,"PU 4 BX 50","lp");
  legu->AddEntry(htype1res_pu20bx25,"PU 20 BX 25","lp");

  //
  // Normalize
  //
  htype1res_pu4bx50->Scale(1/htype1res_pu4bx50->Integral());
  htype1res_pu20bx25->Scale(1/htype1res_pu20bx25->Integral());

  //
  // Save plots
  //
  TCanvas *cres_upar = new TCanvas("cres_upar", "cres_upar",18,40,700,500);
  cres_upar->Range(18.78366,0.5855855,56.2702,0.9308669);
  cres_upar->SetFillColor(0);
  cres_upar->SetBorderMode(0);
  cres_upar->SetBorderSize(2);
  cres_upar->SetFrameBorderMode(0);
  cres_upar->cd();
  hres_upar_rms_pu20bx25->GetYaxis()->SetRangeUser(0,45);
  hres_upar_rms_pu20bx25->Draw("");
  hres_upar_rms_pu4bx50->Draw("same");
  legu->Draw("same");
  cres_upar->Print("cres_upar.png");

  TCanvas *cres_uperp = new TCanvas("cres_uperp", "cres_uperp",18,40,700,500);
  cres_uperp->Range(18.78366,0.5855855,56.2702,0.9308669);
  cres_uperp->SetFillColor(0);
  cres_uperp->SetBorderMode(0);
  cres_uperp->SetBorderSize(2);
  cres_uperp->SetFrameBorderMode(0);
  cres_uperp->cd();
  hres_uperp_rms_pu20bx25->GetYaxis()->SetRangeUser(0,45);
  hres_uperp_rms_pu20bx25->Draw("");
  hres_uperp_rms_pu4bx50->Draw("same");
  legu->Draw("same");
  cres_uperp->Print("cres_uperp.png");

  TCanvas *cmetres = new TCanvas("cmetres", "cmetres",18,40,700,500);
  cmetres->Range(18.78366,0.5855855,56.2702,0.9308669);
  cmetres->SetFillColor(0);
  cmetres->SetBorderMode(0);
  cmetres->SetBorderSize(2);
  cmetres->SetFrameBorderMode(0);
  cmetres->cd();
  htype1res_pu4bx50->Scale(1/htype1res_pu4bx50->Integral());
  htype1res_pu20bx25->Scale(1/htype1res_pu20bx25->Integral());
  htype1res_pu4bx50->Draw("");
  htype1res_pu20bx25->Draw("same");
  leg->Draw("same");
  cmetres->Print("metres.png");

}
