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
#include "Math/GenVector/LorentzVector.h"

// All of this code is a reimplementation of Kevin Sung's code from the 8 TeV analysis.

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

// Convert int to binary value
bool Convert(unsigned int val,bool print=kFALSE)
{
   unsigned int mask = 1 << (sizeof(int) * 8 - 1);
   bool lastDigit;
   for(int i = 0; i < sizeof(int) * 8; i++)
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

void plotWe(const TString inputFileName = "Wenu_p_select.root") {

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
  const Double_t PT_CUT    = 20;
  const Double_t ETA_CUT   = 2.5;
  const Double_t ELE_MASS = 0.000511;
  
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;

  Int_t NVTXBINS = 70; // 70 for Wenu_p, 35 for Wmunu_p

  //
  // Declare variables to read in ntuple
  //
  Int_t nVtx, nEvents;
  TVector2 *vtype1pfMET=0, *vrawpfMET=0, *vgenMET=0;
  LorentzVector *lep0=0, *lep1=0;
  LorentzVector *sc0=0, *sc1=0;
  Float_t isLooseEle0=0, isTightEle0=0;
  Float_t isLooseEle1=0, isTightEle1=0;

  inputTree->SetBranchAddress("nVtx",         &nVtx);        // number of vertices
  inputTree->SetBranchAddress("vtype1pfmet",  &vtype1pfMET); // type-1 corrected pf MET
  inputTree->SetBranchAddress("vrawpfmet",    &vrawpfMET);   // raw pf MET
  inputTree->SetBranchAddress("vgenmet",      &vgenMET);     // generated MET
  inputTree->SetBranchAddress("nEvents",      &nEvents);
  inputTree->SetBranchAddress("lep0",         &lep0);
  inputTree->SetBranchAddress("lep1",         &lep1);
  inputTree->SetBranchAddress("sc0",          &sc0);
  inputTree->SetBranchAddress("sc1",          &sc1);
  inputTree->SetBranchAddress("isLooseEle0",  &isLooseEle0);
  inputTree->SetBranchAddress("isTightEle0",  &isTightEle0);
  inputTree->SetBranchAddress("isLooseEle1",  &isLooseEle1);
  inputTree->SetBranchAddress("isTightEle1",  &isTightEle1);
           
  //
  // Declare histograms
  //
  //
  TH1D *hzmass = new TH1D("hzmass","Dilepton Mass",60,60,120);
        hzmass->GetYaxis()->SetTitle("Events / 1 GeV");
        hzmass->GetYaxis()->SetTitleOffset(1.5);
        hzmass->GetXaxis()->SetTitle("M_{\el\el} [GeV]");
  TH1D *hele0pt = new TH1D("hele0pt","pT of Electron Supercluster Passing Tight Selection",100,0,400);
        hele0pt->GetYaxis()->SetTitle("Events / 4 GeV");
        hele0pt->GetYaxis()->SetTitleOffset(1.5);
        hele0pt->GetXaxis()->SetTitle("pT_{e} [GeV]");
        hele0pt->GetXaxis()->SetTitleOffset(1.2);
  TH1D *hele0eta = new TH1D("hele0eta","Eta of Electron Supercluster Passing Tight Selection",100,-3.0,3.0);
        hele0eta->GetYaxis()->SetTitle("Events / 0.6");
        hele0eta->GetYaxis()->SetTitleOffset(1.5);
        hele0eta->GetXaxis()->SetTitle("#eta_{e}");
        hele0eta->GetXaxis()->SetTitleOffset(1.2);
  TH1D *hele0phi = new TH1D("hele0phi","Phi of Electron Supercluster Passing Tight Selection",100,-3.5,3.5);
        hele0phi->GetYaxis()->SetTitle("Events / 0.7");
        hele0phi->GetYaxis()->SetTitleOffset(1.5);
        hele0phi->GetXaxis()->SetTitle("#phi_{e}");
        hele0phi->GetXaxis()->SetTitleOffset(1.2);
  TH1D *htype1  = new TH1D("htype1","",100,0,150);
        htype1->SetStats(0);
        htype1->SetLineColor(1);
  TH1D *htype1phi = new TH1D("htype1phi","",100,-3.5,3.5);
        htype1phi->SetStats(0);
        htype1phi->SetLineColor(1);
  TH1D *htype1corr  = new TH1D("htype1corr","",100,0,150);
        htype1corr->SetStats(0);
        htype1corr->SetLineColor(4);
  TH1D *htype1phicorr = new TH1D("htype1phicorr","",100,-3.5,3.5);
        htype1phicorr->SetStats(0);
        htype1phicorr->SetLineColor(4);
  TH1D *hraw  = new TH1D("hraw","",100,0,150);
        hraw->GetYaxis()->SetTitle("Events / 1.5 GeV");
        hraw->SetStats(0);
        hraw->SetLineColor(3);
  TH1D *hrawphi = new TH1D("hrawphi","MET(Phi)",100,-3.5,3.5);
        hrawphi->GetYaxis()->SetTitle("Events / 0.7");
        hrawphi->GetYaxis()->SetTitleOffset(1.5);
        hrawphi->GetXaxis()->SetTitle("#phi_{MET}");
        hrawphi->GetXaxis()->SetTitleOffset(1.2);
        hrawphi->SetStats(0);
        hrawphi->SetLineColor(3);
  TH1D *hgen = new TH1D("hgen","MET",100,0,150);
        hgen->GetYaxis()->SetTitle("Events / 1.5 GeV");
        hgen->GetYaxis()->SetTitleOffset(1.5);
        hgen->GetXaxis()->SetTitle("MET [GeV]");
        hgen->GetXaxis()->SetTitleOffset(1.2);
        hgen->SetStats(0);
        hgen->SetLineColor(2);
  TH1D *hgenphi = new TH1D("hgenphi","",100,-3.5,3.5);
        hgenphi->SetStats(0);
        hgenphi->SetLineColor(2);
  TH2D *hmetx = new TH2D("hmetx","MET_{x} v. Number of Vertices",NVTXBINS,0,NVTXBINS,100,-150,150);
        hmetx->GetXaxis()->SetTitle("Number of vertices");
        hmetx->GetYaxis()->SetTitle("MET_{x} [GeV]");
  TH2D *hmety = new TH2D("hmety","MET_{y} v. Number of Vertices",NVTXBINS,0,NVTXBINS,100,-150,150);
        hmety->GetXaxis()->SetTitle("Number of vertices");
        hmety->GetYaxis()->SetTitle("MET_{y} [GeV]");

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

  leg->AddEntry(htype1,"Type 1 Corrected PF","l");
  leg->AddEntry(htype1corr,"Type 1 + xy Shift Corrected PF","l");
  leg->AddEntry(hraw,"Raw PF","l");
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

  legphi->AddEntry(htype1,"Type 1 Corrected PF","l");
  legphi->AddEntry(htype1corr,"Type 1 + xy Shift Corrected PF","l");
  legphi->AddEntry(hraw,"Raw PF","l");
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

  cout << "totalEvents is " << totalEvents << endl;

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

    //
    // SELECTION PROCEDURE:
    //  (1) Find a "tag" electron that passes the tight selection
    //  (2) Find a Supercluster "probe" electron which gives a dilepton mass along with the tag inside the Z-mass window
    //
    Bool_t foundTag=kFALSE, foundProbe=kFALSE;
    Bool_t inECALgap0=kFALSE, inECALgap1=kFALSE;
    // Find tag electron
    if(  fabs(sc0->Eta())>=ECAL_GAP_LOW && fabs(sc0->Eta())<=ECAL_GAP_HIGH  ) inECALgap0=kTRUE; // check ECAL gap
    if(  !inECALgap0 && fabs(sc0->Eta())<=ETA_CUT && sc0->Pt()>=PT_CUT && Convert(isLooseEle0)  ) foundTag=kTRUE;
    if(  foundTag==kFALSE  ) continue; // event not interesting if there is no tag electron
    LorentzVector vTag(lep0->pt(),lep0->eta(),lep0->phi(),ELE_MASS);
    LorentzVector vTagSC(sc0->pt(),sc0->eta(),sc0->phi(),ELE_MASS);
    // Find probe electron
    if(  fabs(sc1->Eta())>=ECAL_GAP_LOW && fabs(sc1->Eta())<=ECAL_GAP_HIGH  ) inECALgap1=kTRUE;
    if(  !inECALgap1 && fabs(sc1->Eta())<=ETA_CUT && sc1->Pt()>=PT_CUT  ) foundProbe=kTRUE;
    if(  foundProbe==kFALSE  ) continue; // event not interesting if there is no probe electron
    LorentzVector vProbe(lep1->pt(),lep1->eta(),lep1->phi(),ELE_MASS);
    LorentzVector vProbeSC(sc1->pt(),sc1->eta(),sc1->phi(),ELE_MASS);
    // Mass window
    LorentzVector vDilep = vTag + vProbe;
    if((vDilep.M()<MASS_LOW) || (vDilep.M()>MASS_TIGHT)) continue;
 
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

    //
    // Fill histograms
    //
    hzmass->Fill(vDilep.M());
/*
    if(passSel0 && !passSel1) {
      hele0pt->Fill(sc0->Pt());
      hele0eta->Fill(sc0->Eta());
      hele0phi->Fill(sc0->Phi());
    } else if (!passSel0 && passSel1) {
      hele0pt->Fill(sc1->Pt());
      hele0eta->Fill(sc1->Eta());
      hele0phi->Fill(sc1->Phi());
    } */
    htype1->Fill(vtype1pfMET->Mod());
    htype1phi->Fill(type1phi);
    htype1corr->Fill(corrMET);
    htype1phicorr->Fill(corrMETphi);
    hraw->Fill(vrawpfMET->Mod());
    hrawphi->Fill(rawphi);
    hgen->Fill(vgenMET->Mod());
    hgenphi->Fill(genphi);

  }

  cout << "nsel is " << nsel << endl;

  //
  // Save plots
  //
  TCanvas* czmass = new TCanvas();
  czmass->cd();
  hzmass->Draw();
  czmass->Print("czmass.png");
/*
  TCanvas* cele0pt = new TCanvas();
  cele0pt->cd();
  cele0pt->SetLogy();
  hele0pt->Draw();
  cele0pt->Print("ele0pt.png");
  cele0pt->Close();
  TCanvas* cele0eta = new TCanvas();
  cele0eta->cd();
  hele0eta->Draw();
  cele0eta->Print("ele0eta.png");
  cele0eta->Close();
  TCanvas* cele0phi = new TCanvas();
  cele0phi->cd();
  hele0phi->Draw();
  cele0phi->Print("ele0phi.png");
  cele0phi->Close();
  TCanvas* cmet = new TCanvas();
  cmet->cd();
  hgen->Draw();
  htype1->Draw("same");
  htype1corr->Draw("same");
  hraw->Draw("same");
  leg->Draw("same");
  cmet->Print("met.png");
  cmet->Close();
  TCanvas* cmetphi = new TCanvas();
  cmetphi->cd();
  hrawphi->Draw();
  htype1phi->Draw("same");
  htype1phicorr->Draw("same");
  hgenphi->Draw("same");
  legphi->Draw("same");
  cmetphi->Print("metphi.png");
//  cmetphi->Close();
  TCanvas* cmetx = new TCanvas();
  cmetx->cd();
  hmetx->Draw();
  cmetx->Print("metx.png");
  cmetx->Close();
  TCanvas* cmety = new TCanvas();
  cmety->cd();
  hmety->Draw();
  cmety->Print("mety.png");
  cmety->Close();
  TCanvas* cmetxfit = new TCanvas();
  cmetxfit->cd();
  hmetxfit->Draw();
  cmetxfit->Print("metxfit.png");
  cmetxfit->Close();
  TCanvas* cmetyfit = new TCanvas();
  cmetyfit->cd();
  hmetyfit->Draw();
  cmetyfit->Print("metyfit.png");
  cmetyfit->Close();
*/
}
