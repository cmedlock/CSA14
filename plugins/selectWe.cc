// -*- C++ -*-
//
// Package:    Test/MiniAnalyzer
// Class:      MiniAnalyzer
// 
/**\class MiniAnalyzer MiniAnalyzer.cc Test/MiniAnalyzer/plugins/MiniAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Catherine Aiko Medlock
//         Created:  Tue, 22 Jul 2014 11:26:17 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TVector2.h>               // 2D vector class

//
// class declaration
//

class selectWe : public edm::EDAnalyzer {
   public:
      explicit selectWe(const edm::ParameterSet&);
      ~selectWe();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      edm::EDGetTokenT<pat::METCollection> metToken_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
};

//
// constants, enums and typedefs
//

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

//
// static data member definitions
//

//--------------------------------------------------------------------------------------------------------------
// Settings 
//============================================================================================================== 

const Double_t PT_CUT    = 25;
const Double_t ETA_CUT   = 2.5;
const Double_t ELE_MASS  = 0.000511;

Double_t nsel=0, nselvar=0;

TString outfilename = TString("WplusToENu_CT10_13TeV-powheg-pythia8_SELECT.root");
TFile *outFile = new TFile(outfilename,"RECREATE");
TTree *outTree = new TTree("Events","Events");

//
// Declare output ntuple variables
//
// Initialze everything to -10 so that it is obvious if there is only 1 electron
//
Int_t   nVtx;
Float_t genVPt0=-10, genVPhi0=-10;
Float_t genVPt1=-10, genVPhi1=-10;
Float_t scale1fb;
Float_t rawpfMETpx, rawpfMETpy; // Will be used in vrawpfMET
TVector2 vtype1pfMET, vrawpfMET, vgenMET;
Int_t   q0=-10, q1=-10;
Int_t   dummynEvents=0, nEvents=0;
Float_t Ele0pt=-10, Ele0eta=-10, Ele0phi=-10; // Will be used in lep0
Float_t Ele1pt=-10, Ele1eta=-10, Ele1phi=-10; // Will be used in lep1
Float_t Elesc0pt=-10, Elesc0eta=-10, Elesc0phi=-10; // Will be used in sc0
Float_t Elesc1pt=-10, Elesc1eta=-10, Elesc1phi=-10; // Will be used in sc1
LorentzVector *lep0=0, *lep1=0;
LorentzVector *sc0=0, *sc1=0;
///// electron specific /////
Float_t pfChIso0=-10, pfGamIso0=-10, pfNeuIso0=-10;
Float_t pfChIso1=-10, pfGamIso1=-10, pfNeuIso1=-10;
Float_t isLooseEle0=-10, isTightEle0=-10;
Float_t isLooseEle1=-10, isTightEle1=-10;

// Compute MC event weight_sel per 1/fb
Double_t weight = 1;
//const Double_t xsec = 1;
//if(xsec>0) weight = 1000.*xsec/(Double_t)eventTree->GetEntries();

//
// constructors and destructor
//
selectWe::selectWe(const edm::ParameterSet& iConfig):
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
   electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
   metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
   pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands")))
{
   //now do what ever initialization is needed

}


selectWe::~selectWe()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
selectWe::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   dummynEvents++;

   isLooseEle0=-10; isTightEle0=-10;
   isLooseEle1=-10; isTightEle1=-10;

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   // good vertex requirement
   if (vertices->empty()) return; // skip the event if no PV found
//   const reco::Vertex &PV = vertices->front();
   nVtx = vertices->size();

   Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);
   if(electrons->size()==0) return;   // skip the event if no electrons found

   unsigned int idxEle0=-1, idxEle1=-1;
   float maxElept=-1, secondmaxElept=-1;
   // find the electron with the highest pt
   for (unsigned int jElectron=0;jElectron < electrons->size();jElectron++) {
     const pat::Electron &ele = (*electrons)[jElectron];
     if     (maxElept==-1)                      { maxElept=ele.pt(); idxEle0=jElectron; }
     else if(maxElept!=-1 && ele.pt()>maxElept) { maxElept=ele.pt(); idxEle0=jElectron; }
   }
   const pat::Electron &Ele0 = (*electrons)[idxEle0];
   q0           = Ele0.charge();
   Ele0pt       = Ele0.pt();
   Ele0eta      = Ele0.eta();
   Ele0phi      = Ele0.phi();
   Elesc0pt     = Ele0.superCluster()->energy()*(Ele0.pt()/Ele0.p());
   Elesc0eta    = Ele0.superCluster()->eta();
   Elesc0phi    = Ele0.superCluster()->phi();
   pfChIso0     = Ele0.chargedHadronIso();
   pfGamIso0    = Ele0.photonIso();
   pfNeuIso0    = Ele0.neutralHadronIso();
   isLooseEle0  = Ele0.electronID("eidLoose");
   isTightEle0  = Ele0.electronID("eidTight");
   // get mother particle information
   const reco::GenParticle* gen = Ele0.genLepton();
   Bool_t foundGenV0 = kFALSE; // some events do not have a generated Z or W, this protects against a seg fault
   if(gen!=NULL) {
     const reco::Candidate* genCand = gen;
     while(genCand!=NULL && genCand->numberOfMothers()==1) {
       genCand = genCand->mother(0);
//       if      (fabs(genCand->pdgId())==23) std::cout << "mother particle of first lepton is a Z with pt " << genCand->pt() << std::endl;
//       else if (fabs(genCand->pdgId())==24) std::cout << "mother particle of first lepton is a W with pt " << genCand->pt() << std::endl;
//       else                                 std::cout << "pdgID of mother particle is " << fabs(genCand->pdgId()) << std::endl;
       if( (  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && !foundGenV0 ) foundGenV0=kTRUE;
       genVPt0  = genCand->pt();
       genVPhi0 = genCand->phi();
     }
     // mother particle should always be either a Z or a W
     // if it isn't, trace back to the first daughter particle that is either a Z or a W
     if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && foundGenV0) {
//       std::cout << "***ERROR - pdgID of mother particle is " << fabs(genCand->pdgId()) << std::endl;
//       std::cout << "********tracing back to daughter particles..." << std::endl;
       while( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) ) {
         genCand = genCand->daughter(0);
//         if      (fabs(genCand->pdgId())==23) std::cout << "daughter particle is a Z with pt " << genCand->pt() << std::endl;
//         else if (fabs(genCand->pdgId())==24) std::cout << "daughter particle is a W with pt " << genCand->pt() << std::endl;
//         else                                 std::cout << "pdgID of mother particle is " << fabs(genCand->pdgId()) << std::endl;
         genVPt0  = genCand->pt();
         genVPhi0 = genCand->phi();
       } // end of while
     } // end of if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) )
     else if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && !foundGenV0) {
//       std::cout << "***ERROR - pdgID of mother particle is " << fabs(genCand->pdgId()) << " and no Z or W found" << std::endl;
//       std::cout << "********resetting genVPt0 and genVPhi0..." << std::endl;
       genVPt0  = -10;
       genVPhi0 = -10;
     } // end of else if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && foundGenV==kFALSE)
   } // end of gen!=NULL
//   std::cout << "Ele0pt is " << Ele0pt<< " and genVPt0 is " << genVPt0 << std::endl;
   // if there are > 1 muons, find the 2 with the highest pt's
   if(electrons->size() > 1) {
     for (unsigned int kElectron=0; kElectron < electrons->size();kElectron++) {
       const pat::Electron &ele = (*electrons)[kElectron];
       if(kElectron==idxEle0)                                                       continue;
       else if(kElectron!=idxEle0 && secondmaxElept==-1)                            { secondmaxElept=ele.pt(); idxEle1=kElectron; }
       else if(kElectron!=idxEle0 && secondmaxElept!=-1 && ele.pt()>secondmaxElept) { secondmaxElept=ele.pt(); idxEle1=kElectron; }
     }
     const pat::Electron &Ele1 = (*electrons)[idxEle1];
     q1           = Ele1.charge();
     Ele1pt       = Ele1.pt();
     Ele1eta      = Ele1.eta();
     Ele1phi      = Ele1.phi();
     Elesc1pt     = Ele1.superCluster()->energy()*(Ele1.pt()/Ele1.p());
     Elesc1eta    = Ele1.superCluster()->eta();
     Elesc1phi    = Ele1.superCluster()->phi();
     pfChIso1     = Ele1.chargedHadronIso();
     pfGamIso1    = Ele1.photonIso();
     pfNeuIso1    = Ele1.neutralHadronIso();
     isLooseEle1  = Ele1.electronID("eidLoose");
     isTightEle1  = Ele1.electronID("eidTight");
     // get mother particle information
     const reco::GenParticle* gen = Ele1.genLepton();
     Bool_t foundGenV1 = kFALSE; // some events do not have a generated Z or W, this protects against a seg fault
     if(gen!=NULL) {
       const reco::Candidate* genCand = gen;
       while(genCand!=NULL && genCand->numberOfMothers()==1) {
         genCand = genCand->mother(0);
//         if      (fabs(genCand->pdgId())==23) std::cout << "mother particle of second lepton is a Z with pt " << genCand->pt() << std::endl;
//         else if (fabs(genCand->pdgId())==24) std::cout << "mother particle of second lepton is a W with pt " << genCand->pt() << std::endl;
//         else                                 std::cout << "pdgID of mother particle is " << fabs(genCand->pdgId()) << std::endl;
         if( (  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && !foundGenV1 ) foundGenV1=kTRUE;
         genVPt1  = genCand->pt();
         genVPhi1 = genCand->phi();
       }
       // mother particle should always be either a Z or a W
       // if it isn't, trace back to the first daughter particle that is either a Z or a W
       if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && foundGenV1) {
//         std::cout << "***ERROR - pdgID of mother particle is " << fabs(genCand->pdgId()) << std::endl;
//         std::cout << "********tracing back to daughter particles..." << std::endl;
         while( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) ) {
           genCand = genCand->daughter(0);
//           if      (fabs(genCand->pdgId())==23) std::cout << "daughter particle is a Z with pt " << genCand->pt() << std::endl;
//           else if (fabs(genCand->pdgId())==24) std::cout << "daughter particle is a W with pt " << genCand->pt() << std::endl;
//           else                                 std::cout << "pdgID of mother particle is " << fabs(genCand->pdgId()) << std::endl;
           genVPt1  = genCand->pt();
           genVPhi1 = genCand->phi();
         } // end of while
       } // end of if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) )
       else if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && !foundGenV1) {
//         std::cout << "***ERROR - pdgID of mother particle is " << fabs(genCand->pdgId()) << " and no Z or W found" << std::endl;
//         std::cout << "********resetting genVPt1 and genVPhi1..." << std::endl;
         genVPt1  = -10;
         genVPhi1 = -10;
       } // end of else if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && foundGenV==kFALSE)
     } // end of gen!=NULL
   }

   nsel    += weight;
   nselvar += weight*weight;
   scale1fb = weight;

   LorentzVector vLep0(Ele0pt, Ele0eta, Ele0phi, ELE_MASS); lep0 = &vLep0;
   LorentzVector vLep1(Ele1pt, Ele1eta, Ele1phi, ELE_MASS); lep1 = &vLep1;
   LorentzVector vSC0(Elesc0pt, Elesc0eta, Elesc0phi, ELE_MASS); sc0 = &vSC0;
   LorentzVector vSC1(Elesc1pt, Elesc1eta, Elesc1phi, ELE_MASS); sc1 = &vSC1;

   edm::Handle<pat::METCollection> mets;
   iEvent.getByToken(metToken_, mets);
   const pat::MET &met = mets->front();

   edm::Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);
   rawpfMETpx=0;
   rawpfMETpy=0;
   // loop on pf candidates to calculate sum of Et's
   for(unsigned int jcand=0; jcand < pfs->size(); jcand++) {
     const pat::PackedCandidate &pf = (*pfs)[jcand];
     rawpfMETpx -= pf.px();
     rawpfMETpy -= pf.py();
   }

//   pat::MET::UncorrectionType ix;
//   ix=pat::MET::uncorrALL;
//   pfMETPhi   = met.uncorrectedPhi(ix);

   vtype1pfMET.Set(met.px(),met.py());
   vrawpfMET.Set(rawpfMETpx,rawpfMETpy);
   vgenMET.Set(met.genMET()->px(),met.genMET()->py());

   //
   // Fill tree
   //
   outTree->Fill();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
selectWe::beginJob()
{
   // Create output directory
//   gSystem->mkdir(outputDir,kTRUE); // AIKO - THIS NEEDS TO BE FIXED (NOT URGENT)
//   gSystem->mkdir(ntupDir,kTRUE);

  //
  // Set up output ntuple
  //

  outTree->Branch("nVtx",          &nVtx,          "nVtx/I");         // number of vertices
  outTree->Branch("genVPt0",       &genVPt0,       "genVPt0/F");      // GEN boson pT  (signal MC)  (mother of first lepton )
  outTree->Branch("genVPhi0",      &genVPhi0,      "genVPhi0/F");     // GEN boson phi (signal MC)  (mother of first lepton )
  outTree->Branch("genVPt1",       &genVPt1,       "genVPt1/F");      // GEN boson pT  (signal MC)  (mother of second lepton)
  outTree->Branch("genVPhi1",      &genVPhi1,      "genVPhi1/F");     // GEN boson phi (signal MC)  (mother of second lepton)
  outTree->Branch("scale1fb",      &scale1fb,      "scale1fb/F");     // event weight per 1/fb (MC)
  outTree->Branch("vtype1pfmet",   "TVector2",     &vtype1pfMET);     // type-1 corrected pf MET
  outTree->Branch("vrawpfmet",     "TVector2",     &vrawpfMET);       // raw pf MET
  outTree->Branch("vgenmet",       "TVector2",     &vgenMET);         // generated MET
  outTree->Branch("q0",            &q0,            "q0/I");           // lepton charge (first lepton )
  outTree->Branch("q1",            &q1,            "q1/I");           // lepton charge (second lepton)
  outTree->Branch("nEvents",       &nEvents,       "nEvents/I");      // events in MC file
  outTree->Branch("lep0", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep0);   // lepton 4-vector (first lepton )
  outTree->Branch("lep1", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep1);   // lepton 4-vector (second lepton)
  outTree->Branch("sc0",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sc0);    // supercluster 4-vector (first lepton )
  outTree->Branch("sc1",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sc1);    // supercluster 4-vector (second lepton)
  ///// electron specific /////
  outTree->Branch("pfChIso0",     &pfChIso0,    "pfChIso0/F");      // PF charged hadron isolation of lepton (first lepton )
  outTree->Branch("pfGamIso0",    &pfGamIso0,   "pfGamIso0/F");     // PF photon isolation of lepton         (first lepton )
  outTree->Branch("pfNeuIso0",    &pfNeuIso0,   "pfNeuIso0/F");     // PF neutral hadron isolation of lepton (first lepton )
  outTree->Branch("pfChIso1",     &pfChIso1,    "pfChIso1/F");      // PF charged hadron isolation of lepton (second lepton)
  outTree->Branch("pfGamIso1",    &pfGamIso1,   "pfGamIso1/F");     // PF photon isolation of lepton         (second lepton)
  outTree->Branch("pfNeuIso1",    &pfNeuIso1,   "pfNeuIso1/F");     // PF neutral hadron isolation of lepton (second lepton)
  outTree->Branch("isLooseEle0",  &isLooseEle0, "isLooseEle0/F");   // loose electron ID (first lepton )
  outTree->Branch("isTightEle0",  &isTightEle0, "isTightEle0/F");   // tight electron ID (first lepton )
  outTree->Branch("isLooseEle1",  &isLooseEle1, "isLooseEle1/F");   // loose electron ID (second lepton)
  outTree->Branch("isTightEle1",  &isTightEle1, "isTightEle1/F");   // tight electron ID (second lepton)
}

// ------------ method called once each job just after ending the event loop  ------------
void 
selectWe::endJob() 
{

   nVtx          = -1;
   genVPt0       = -1; genVPhi0 = -1;
   genVPt1       = -1; genVPhi1 = -1;
   scale1fb      = -1;
   vtype1pfMET.Set(-1.0,-1.0);
   vrawpfMET.Set(-1.0,-1.0);
   vgenMET.Set(-1.0,-1.0);
   q0            = -1; q1 = -1;
   nEvents   = dummynEvents;
   LorentzVector dummyLep(-1, -1, -1, -1);
   lep0      = &dummyLep;
   lep1      = &dummyLep;
   sc0       = &dummyLep;
   sc1       = &dummyLep;
   pfChIso0  = -1; pfGamIso0 = -1; pfNeuIso0 = -1;
   pfChIso1  = -1; pfGamIso1 = -1; pfNeuIso1 = -1;
   isLooseEle0 = kFALSE; isTightEle0 = kFALSE;
   isLooseEle1 = kFALSE; isTightEle1 = kFALSE;
   
   outTree->Fill();
   std::cout << nsel << " +/- " << sqrt(nselvar) << " per 1/fb" << std::endl;
   std::cout << "endJob: nEvents is " << nEvents << std::endl;
   outFile->Write();
   outFile->Close();

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================

  std::cout << "*" << std::endl;
  std::cout << "* SUMMARY" << std::endl;
  std::cout << "*--------------------------------------------------" << std::endl;
  std::cout << "W -> e nu" << std::endl;
  std::cout << " pT > " << PT_CUT << std::endl;
  std::cout << " |eta| < " << ETA_CUT << std::endl;
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "  <> Output saved in " << outfilename << "/" << std::endl;
  std::cout << std::endl;
}

// ------------ method called when starting to processes a run  ------------
/*
void 
selectWe::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
selectWe::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
selectWe::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
selectWe::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
selectWe::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(selectWe);
