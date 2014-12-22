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
#include "DataFormats/PatCandidates/interface/Muon.h"
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

class selectWm : public edm::EDAnalyzer {
   public:
      explicit selectWm(const edm::ParameterSet&);
      ~selectWm();

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
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
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

const Double_t PT_CUT    = 20;
const Double_t ETA_CUT   = 2.4;
const Double_t MUON_MASS = 0.105658369;

Double_t nselWm=0, nselWmvar=0;

TString outFilename_Wm = TString("Wmunu_select.root");
TFile *outFile_Wm = new TFile();
TTree *outTree_Wm = new TTree();

//
// Declare output ntuple variables
//
// Initialize everything to -10 so that it is obvious if an object (such
// as a Z or W generator level mother particle) doesn't exist
//
Int_t   nVtxWm;
Float_t genVPt0Wm=-10, genVPhi0Wm=-10;
Float_t genVPt1Wm=-10, genVPhi1Wm=-10;
Float_t scale1fbWm;
Float_t rawpfMETpxWm, rawpfMETpyWm; // Will be used in vrawpfMETWm
TVector2 vtype1pfMETWm, vrawpfMETWm, vgenMETWm;
Int_t   q0Wm=-10, q1Wm=-10;
Int_t   dummynEventsWmWm=0, nEventsWm=0;
Float_t Mu0pt=-10, Mu0eta=-10, Mu0phi=-10; // Will be used in lep0Wm
Float_t Mu1pt=-10, Mu1eta=-10, Mu1phi=-10; // Will be used in lep1Wm
LorentzVector *lep0Wm=0, *lep1Wm=0;
///// muon specific /////
Float_t pfChIso0Wm=-10, pfGamIso0Wm=-10, pfNeuIso0Wm=-10;
Float_t pfChIso1Wm=-10, pfGamIso1Wm=-10, pfNeuIso1Wm=-10;
Bool_t  isLooseMuon0=kFALSE, isSoftMuon0=kFALSE, isTightMuon0=kFALSE;
Bool_t  isLooseMuon1=kFALSE, isSoftMuon1=kFALSE, isTightMuon1=kFALSE;

// Compute MC event weightWm_sel per 1/fb
Double_t weightWm = 1;
//const Double_t xsec = 1;
//if(xsec>0) weightWm = 1000.*xsec/(Double_t)eventTree->GetEntries();

//
// constructors and destructor
//
selectWm::selectWm(const edm::ParameterSet& iConfig):
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
   muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
   metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
   pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands")))
{
   //now do what ever initialization is needed

}


selectWm::~selectWm()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
selectWm::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   dummynEventsWmWm++;

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   // good vertex requirement
   if (vertices->empty()) return; // skip the event if no PV found
   const reco::Vertex &PV = vertices->front();
   nVtxWm = vertices->size();

   Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);
   if(muons->size()==0) return;   // skip the event if no muons found
   unsigned int idxMu0=-1, idxMu1=-1;
   float maxMupt=-1 , secondmaxMupt=-1;
   // find the muon with the highest pT
   for (unsigned int jMuon=0;jMuon < muons->size();jMuon++) {
     const pat::Muon &mu = (*muons)[jMuon];
     if     (maxMupt==-1)                    { maxMupt=mu.pt(); idxMu0=jMuon; }
     else if(maxMupt!=-1 && mu.pt()>maxMupt) { maxMupt=mu.pt(); idxMu0=jMuon; }
   }
   const pat::Muon &Mu0 = (*muons)[idxMu0];
   q0Wm           = Mu0.charge();
   Mu0pt        = Mu0.pt();
   Mu0eta       = Mu0.eta();
   Mu0phi       = Mu0.phi();
   pfChIso0Wm     = Mu0.chargedHadronIso();
   pfGamIso0Wm    = Mu0.photonIso();
   pfNeuIso0Wm    = Mu0.neutralHadronIso();
   isLooseMuon0 = Mu0.isLooseMuon();
   isSoftMuon0  = Mu0.isSoftMuon(PV);
   isTightMuon0 = Mu0.isTightMuon(PV);
   // get mother particle information
   const reco::GenParticle* gen = Mu0.genLepton();
   Bool_t foundGenV0 = kFALSE; // some events do not have a generated Z or W, this protects against a seg fault
   if(gen!=NULL) {
     const reco::Candidate* genCand = gen;
     while(genCand!=NULL && genCand->numberOfMothers()==1) {
       genCand = genCand->mother(0);
       if( (  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && !foundGenV0 ) foundGenV0=kTRUE;
       genVPt0Wm  = genCand->pt();
       genVPhi0Wm = genCand->phi();
     }
     // mother particle should always be either a Z or a W
     // if it isn't, trace back to the first daughter particle that is either a Z or a W
     if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && foundGenV0) {
       while( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) ) {
         genCand = genCand->daughter(0);
         genVPt0Wm  = genCand->pt();
         genVPhi0Wm = genCand->phi();
       }
     }
     // if there is no mother particle that is a Z or a W, mark this by saving the mother particle pT and phi both as -10
     else if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && !foundGenV0) {
       genVPt0Wm  = -10;
       genVPhi0Wm = -10;
     }
   }
   // if there is > 1 muon, find the one with the second highest pT
   if(muons->size() > 1) {
     for (unsigned int kMuon=0; kMuon < muons->size();kMuon++) {
       const pat::Muon &mu = (*muons)[kMuon];
       if(kMuon==idxMu0)                                                      continue;
       else if(kMuon!=idxMu0 && secondmaxMupt==-1)                          { secondmaxMupt=mu.pt(); idxMu1=kMuon; }
       else if(kMuon!=idxMu0 && secondmaxMupt!=-1 && mu.pt()>secondmaxMupt) { secondmaxMupt=mu.pt(); idxMu1=kMuon; }
     }
     const pat::Muon &Mu1 = (*muons)[idxMu1];
     q1Wm           = Mu1.charge();
     Mu1pt        = Mu1.pt();
     Mu1eta       = Mu1.eta();
     Mu1phi       = Mu1.phi();
     pfChIso1Wm     = Mu1.chargedHadronIso();
     pfGamIso1Wm    = Mu1.photonIso();
     pfNeuIso1Wm    = Mu1.neutralHadronIso();
     isLooseMuon1 = Mu1.isLooseMuon();
     isSoftMuon1  = Mu1.isSoftMuon(PV);
     isTightMuon1 = Mu1.isTightMuon(PV);
     // get mother particle information
     const reco::GenParticle* gen = Mu1.genLepton();
     Bool_t foundGenV1 = kFALSE; // some events do not have a generated Z or W, this protects against a seg fault
     if(gen!=NULL) {
       const reco::Candidate* genCand = gen;
       while(genCand!=NULL && genCand->numberOfMothers()==1) {
         genCand = genCand->mother(0);
         if( (  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && !foundGenV1 ) foundGenV1=kTRUE;
         genVPt1Wm  = genCand->pt();
         genVPhi1Wm = genCand->phi();
       }
       // mother particle should always be either a Z or a W
       // if it isn't, trace back to the first daughter particle that is either a Z or a W
       if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && foundGenV1) {
         while( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) ) {
           genCand = genCand->daughter(0);
           genVPt1Wm  = genCand->pt();
           genVPhi1Wm = genCand->phi();
         }
       }
       // if there is no mother particle that is a Z or a W, mark this by saving the mother particle pT and phi both as -10
       else if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && !foundGenV1) {
         genVPt1Wm  = -10;
         genVPhi1Wm = -10;
       }
     }
   }

   nselWm    += weightWm;
   nselWmvar += weightWm*weightWm;
   scale1fbWm = weightWm;

   LorentzVector vlep0Wm(Mu0pt, Mu0eta, Mu0phi, MUON_MASS); lep0Wm = &vlep0Wm;
   LorentzVector vlep1Wm(Mu1pt, Mu1eta, Mu1phi, MUON_MASS); lep1Wm = &vlep1Wm;

   edm::Handle<pat::METCollection> mets;
   iEvent.getByToken(metToken_, mets);
   const pat::MET &met = mets->front();

   edm::Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);
   rawpfMETpxWm=0;
   rawpfMETpyWm=0;
   // loop on pf candidates to calculate sum of Et's
   for(unsigned int jcand=0; jcand < pfs->size(); jcand++) {
     const pat::PackedCandidate &pf = (*pfs)[jcand];
     rawpfMETpxWm -= pf.px();
     rawpfMETpyWm -= pf.py();
   }

   vtype1pfMETWm.Set(met.px(),met.py());
   vrawpfMETWm.Set(rawpfMETpxWm,rawpfMETpyWm);
   vgenMETWm.Set(met.genMET()->px(),met.genMET()->py());

   //
   // Fill tree
   //
   outTree_Wm->Fill();

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
selectWm::beginJob()
{
   // Create output directory
//   gSystem->mkdir(outputDir,kTRUE);
//   gSystem->mkdir(ntupDir,kTRUE);

  //
  // Set up output ntuple
  //

  outFile_Wm = new TFile(outFilename_Wm,"RECREATE");
  outTree_Wm = new TTree("Events","Events");

  outTree_Wm->Branch("nVtx",        &nVtxWm,         "nVtx/I");     // number of vertices
  outTree_Wm->Branch("genVPt0",     &genVPt0Wm,      "genVPt0/F");  // GEN boson pT  (signal MC)  (mother of first lepton )
  outTree_Wm->Branch("genVPhi0",    &genVPhi0Wm,     "genVPhi0/F"); // GEN boson phi (signal MC)  (mother of first lepton )
  outTree_Wm->Branch("genVPt1",     &genVPt1Wm,      "genVPt1/F");  // GEN boson pT  (signal MC)  (mother of second lepton)
  outTree_Wm->Branch("genVPhi1",    &genVPhi1Wm,     "genVPhi1/F"); // GEN boson phi (signal MC)  (mother of second lepton)
  outTree_Wm->Branch("scale1fb",    &scale1fbWm,     "scale1fb/F"); // event weightWm per 1/fb (MC)
  outTree_Wm->Branch("vtype1pfMET", "TVector2",    &vtype1pfMETWm); // type-1 corrected pf MET
  outTree_Wm->Branch("vrawpfMET",   "TVector2",    &vrawpfMETWm);   // raw pf MET
  outTree_Wm->Branch("vgenMET",     "TVector2",    &vgenMETWm);     // generated MET
  outTree_Wm->Branch("q0",          &q0Wm,           "q0/I");       // lepton charge (first lepton )
  outTree_Wm->Branch("q1",          &q1Wm,           "q1/I");       // lepton charge (second lepton)
  outTree_Wm->Branch("nEvents",     &nEventsWm,      "nEvents/I");  // events in MC file
  outTree_Wm->Branch("lep0", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep0Wm);   // lepton 4-vector (first lepton )
  outTree_Wm->Branch("lep1", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep1Wm);   // lepton 4-vector (second lepton)
  ///// muon specific /////
  outTree_Wm->Branch("pfChIso0",     &pfChIso0Wm,    "pfChIso0/F");      // PF charged hadron isolation of lepton (first lepton )
  outTree_Wm->Branch("pfGamIso0",    &pfGamIso0Wm,   "pfGamIso0/F");     // PF photon isolation of lepton         (first lepton )
  outTree_Wm->Branch("pfNeuIso0",    &pfNeuIso0Wm,   "pfNeuIso0/F");     // PF neutral hadron isolation of lepton (first lepton )
  outTree_Wm->Branch("pfChIso1",     &pfChIso1Wm,    "pfChIso1/F");      // PF charged hadron isolation of lepton (second lepton)
  outTree_Wm->Branch("pfGamIso1",    &pfGamIso1Wm,   "pfGamIso1/F");     // PF photon isolation of lepton         (second lepton)
  outTree_Wm->Branch("pfNeuIso1",    &pfNeuIso1Wm,   "pfNeuIso1/F");     // PF neutral hadron isolation of lepton (second lepton)
  outTree_Wm->Branch("isLooseMuon0", &isLooseMuon0, "isLooseMuon0/O"); // loose muon ID (first lepton )
  outTree_Wm->Branch("isSoftMuon0",  &isSoftMuon0,  "isSoftMuon0/O");  // soft  muon ID (first lepton )
  outTree_Wm->Branch("isTightMuon0", &isTightMuon0, "isTightMuon0/O"); // tight muon ID (first lepton )
  outTree_Wm->Branch("isLooseMuon1", &isLooseMuon1, "isLooseMuon1/O"); // loose muon ID (second lepton)
  outTree_Wm->Branch("isSoftMuon1",  &isSoftMuon1,  "isSoftMuon1/O");  // soft  muon ID (second lepton)
  outTree_Wm->Branch("isTightMuon1", &isTightMuon1, "isTightMuon1/O"); // tight muon ID (second lepton)

}

// ------------ method called once each job just after ending the event loop  ------------
void 
selectWm::endJob() 
{

   // The only information in the last entry of the tree is the total number of events
   // that were processed (not the same as the total number of events selected).
   std::cout << nselWm << " +/- " << sqrt(nselWmvar) << " per 1/fb" << std::endl;
   std::cout << "endJob: nEventsWm is " << nEventsWm << std::endl;

   nVtxWm      = -1;
   genVPt0Wm   = -1; genVPhi0Wm = -1;
   genVPt1Wm   = -1; genVPhi1Wm = -1;
   scale1fbWm  = -1;
   vtype1pfMETWm.Set(-1.0,-1.0);
   vrawpfMETWm.Set(-1.0,-1.0);
   vgenMETWm.Set(-1.0,-1.0);
   q0Wm        = -1; q1Wm = -1;
   nEventsWm   = dummynEventsWmWm;
   LorentzVector dummyLep(-1, -1, -1, -1);
   lep0Wm      = &dummyLep;
   lep1Wm      = &dummyLep;
   pfChIso0Wm  = -1; pfGamIso0Wm = -1; pfNeuIso0Wm = -1;
   pfChIso1Wm  = -1; pfGamIso1Wm = -1; pfNeuIso1Wm = -1;
   isLooseMuon0 = kFALSE; isSoftMuon0 = kFALSE; isTightMuon0 = kFALSE;
   isLooseMuon1 = kFALSE; isSoftMuon0 = kFALSE; isTightMuon0 = kFALSE;
   
   outTree_Wm->Fill();

   std::cout << nselWm << " +/- " << sqrt(nselWmvar) << " per 1/fb" << std::endl;
   std::cout << "endJob: nEventsWm is " << nEventsWm << std::endl;
   outFile_Wm->Write();
   outFile_Wm->Close();

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================

  std::cout << "*" << std::endl;
  std::cout << "* SUMMARY" << std::endl;
  std::cout << "*--------------------------------------------------" << std::endl;
  std::cout << "W -> mu nu" << std::endl;
  std::cout << " pT > " << PT_CUT << std::endl;
  std::cout << " |eta| < " << ETA_CUT << std::endl;
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "  <> Output saved in " << outFilename_Wm << "/" << std::endl;
  std::cout << std::endl;
}

// ------------ method called when starting to processes a run  ------------
/*
void 
selectWm::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
selectWm::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
selectWm::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
selectWm::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
selectWm::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(selectWm);
