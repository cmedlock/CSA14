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

class selectZee : public edm::EDAnalyzer {
   public:
      explicit selectZee(const edm::ParameterSet&);
      ~selectZee();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
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

// Convert int to binary value
bool Convert_Zee(unsigned int val,bool print=kFALSE)
{
   unsigned int mask = 1 << (sizeof(int) * 8 - 1);
   bool lastDigit;
   for(unsigned int i = 0; i < sizeof(int) * 8; i++)
   {
      if( (val & mask) == 0 ) {
        lastDigit=0;
      } else {
        lastDigit=1;
      }
      mask >>= 1;
   }
   return lastDigit;
}

//
// static data member definitions
//

//--------------------------------------------------------------------------------------------------------------
// Settings 
//============================================================================================================== 

const Double_t MASS_LOW  = 40;
const Double_t MASS_HIGH = 200;
const Double_t PT_CUT    = 20;
const Double_t ETA_CUT   = 2.5;
const Double_t ELE_MASS  = 0.000511;

const Double_t ECAL_GAP_LOW  = 1.4442;
const Double_t ECAL_GAP_HIGH = 1.566;

Double_t nsel_Zee=0, nsel_Zeevar_Zee=0;

TString outFilename_Zee = TString("Zee_select.root");
TFile *outFile_Zee = new TFile();
TTree *outTree_Zee = new TTree();

//
// Declare output ntuple variables
//
// Initialize everything to -10 so that it is obvious if an object (such
// as a Z or W generator level mother particle) doesn't exist
//
Int_t   nVtx_Zee;
Float_t genVPt1_Zee=-10, genVPhi1_Zee=-10, genVy1_Zee=-10, genVMass1_Zee=-10;
Float_t genVPt2_Zee=-10, genVPhi2_Zee=-10, genVy2_Zee=-10, genVMass2_Zee=-10;
Float_t scale1fb_Zee;
Float_t rawpfMETpx_Zee, rawpfMETpy_Zee; // Will be used in vrawpfMET_Zee
TVector2 vtype1pfMET_Zee, vrawpfMET_Zee, vgenMET_Zee;
Float_t u1_Zee=-10, u2_Zee=-10;
Int_t   q1_Zee=-10, q2_Zee=-10;
Int_t   dummynEvents_Zee=0, nEvents_Zee=0;
Float_t Ele1pt_Zee=-10, Ele1eta_Zee=-10, Ele1phi_Zee=-10; // Will be used in lep1_Zee
Float_t Ele2pt_Zee=-10, Ele2eta_Zee=-10, Ele2phi_Zee=-10; // Will be used in lep2_Zee
Float_t Elesc1pt_Zee=-10, Elesc1eta_Zee=-10, Elesc1phi_Zee=-10; // Will be used in sc1_Zee
Float_t Elesc2pt_Zee=-10, Elesc2eta_Zee=-10, Elesc2phi_Zee=-10; // Will be used in sc2_Zee
LorentzVector *dilep_Zee=0;
LorentzVector *lep1_Zee=0, *lep2_Zee=0;
LorentzVector *sc1_Zee=0, *sc2_Zee=0;
///// electron specific /////
Float_t pfChIso1_Zee=-10, pfGamIso1_Zee=-10, pfNeuIso1_Zee=-10;
Float_t pfChIso2_Zee=-10, pfGamIso2_Zee=-10, pfNeuIso2_Zee=-10;
Float_t isLooseEle1_Zee=-10, isTightEle1_Zee=-10;
Float_t isLooseEle2_Zee=-10, isTightEle2_Zee=-10;

// Compute MC event weight_Zee_sel per 1/fb
Double_t weight_Zee = 1;
//const Double_t xsec = 1;
//if(xsec>0) weight_Zee = 1000.*xsec/(Double_t)eventTree->GetEntries();

//
// constructors and destructor
//
selectZee::selectZee(const edm::ParameterSet& iConfig):
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
   electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
   metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
   pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands")))
{
   //now do what ever initialization is needed
}


selectZee::~selectZee()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
selectZee::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   isLooseEle1_Zee=-10; isTightEle1_Zee=-10;
   isLooseEle2_Zee=-10; isTightEle2_Zee=-10;

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   // good vertex requirement
   if (vertices->empty()) return; // skip the event if no PV found
//   const reco::Vertex &PV = vertices->front();
   nVtx_Zee = vertices->size();

   Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);
   if(electrons->size()==0) return;   // skip the event if no electrons found

   dummynEvents_Zee++;

   unsigned int idxEle1=-1, idxEle2=-1;
   float maxElept=-1, secondmaxElept=-1;
   // find the electron with the highest pT
   for (unsigned int jElectron=0;jElectron < electrons->size();jElectron++) {
     const pat::Electron &ele = (*electrons)[jElectron];
     if     (maxElept==-1)                      { maxElept=ele.pt(); idxEle1=jElectron; }
     else if(maxElept!=-1 && ele.pt()>maxElept) { maxElept=ele.pt(); idxEle1=jElectron; }
   }
   const pat::Electron &Ele1 = (*electrons)[idxEle1];
   q1_Zee           = Ele1.charge();
   Ele1pt_Zee       = Ele1.pt();
   Ele1eta_Zee      = Ele1.eta();
   Ele1phi_Zee      = Ele1.phi();
   Elesc1pt_Zee     = Ele1.superCluster()->energy()*(Ele1.pt()/Ele1.p());
   Elesc1eta_Zee    = Ele1.superCluster()->eta();
   Elesc1phi_Zee    = Ele1.superCluster()->phi();
   pfChIso1_Zee     = Ele1.chargedHadronIso();
   pfGamIso1_Zee    = Ele1.photonIso();
   pfNeuIso1_Zee    = Ele1.neutralHadronIso();
   isLooseEle1_Zee  = Ele1.electronID("eidLoose");
   isTightEle1_Zee  = Ele1.electronID("eidTight");
   // get mother particle information
   const reco::GenParticle* gen = Ele1.genLepton();
   Bool_t foundGenV0 = kFALSE; // some events do not have a generated Z or W, this protects against a seg fault
   if(gen!=NULL) {
     const reco::Candidate* genCand = gen;
     while(genCand!=NULL && genCand->numberOfMothers()==1) {
       genCand = genCand->mother(0);
       if( (  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && !foundGenV0 ) foundGenV0=kTRUE;
       genVPt1_Zee   = genCand->pt();
       genVPhi1_Zee  = genCand->phi();
       genVy1_Zee    = genCand->y();
       genVMass1_Zee = genCand->mass();
     }
     // mother particle should always be either a Z or a W
     // if it isn't, trace back to the first daughter particle that is either a Z or a W
     if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && foundGenV0) {
       while( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) ) {
         genCand = genCand->daughter(0);
         genVPt1_Zee   = genCand->pt();
         genVPhi1_Zee  = genCand->phi();
         genVy1_Zee    = genCand->y();
         genVMass1_Zee = genCand->mass();
       }
     }
     // if there is no mother particle that is a Z or a W, mark this by saving the mother particle pT and phi both as -10
     else if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && !foundGenV0) {
       genVPt1_Zee   = -10;
       genVPhi1_Zee  = -10;
       genVy1_Zee    = -10;
       genVMass1_Zee = -10;
     }
   }
   // if there is > 1 electron, find the one with the second highest pT
   if(electrons->size() > 1) {
     for (unsigned int kElectron=0; kElectron < electrons->size();kElectron++) {
       const pat::Electron &ele = (*electrons)[kElectron];
       if(kElectron==idxEle1)                                                       continue;
       else if(kElectron!=idxEle1 && secondmaxElept==-1)                            { secondmaxElept=ele.pt(); idxEle2=kElectron; }
       else if(kElectron!=idxEle1 && secondmaxElept!=-1 && ele.pt()>secondmaxElept) { secondmaxElept=ele.pt(); idxEle2=kElectron; }
     }
     const pat::Electron &Ele2 = (*electrons)[idxEle2];
     q2_Zee           = Ele2.charge();
     Ele2pt_Zee       = Ele2.pt();
     Ele2eta_Zee      = Ele2.eta();
     Ele2phi_Zee      = Ele2.phi();
     Elesc2pt_Zee     = Ele2.superCluster()->energy()*(Ele2.pt()/Ele2.p());
     Elesc2eta_Zee    = Ele2.superCluster()->eta();
     Elesc2phi_Zee    = Ele2.superCluster()->phi();
     pfChIso2_Zee     = Ele2.chargedHadronIso();
     pfGamIso2_Zee    = Ele2.photonIso();
     pfNeuIso2_Zee    = Ele2.neutralHadronIso();
     isLooseEle2_Zee  = Ele2.electronID("eidLoose");
     isTightEle2_Zee  = Ele2.electronID("eidTight");
     // get mother particle information
     const reco::GenParticle* gen = Ele2.genLepton();
     Bool_t foundGenV1 = kFALSE; // some events do not have a generated Z or W, this protects against a seg fault
     if(gen!=NULL) {
       const reco::Candidate* genCand = gen;
       while(genCand!=NULL && genCand->numberOfMothers()==1) {
         genCand = genCand->mother(0);
         if( (  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && !foundGenV1 ) foundGenV1=kTRUE;
         genVPt2_Zee   = genCand->pt();
         genVPhi2_Zee  = genCand->phi();
         genVy2_Zee    = genCand->y();
         genVMass2_Zee = genCand->mass();
       }
       // mother particle should always be either a Z or a W
       // if it isn't, trace back to the first daughter particle that is either a Z or a W
       if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && foundGenV1) {
         while( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) ) {
           genCand = genCand->daughter(0);
           genVPt2_Zee   = genCand->pt();
           genVPhi2_Zee  = genCand->phi();
           genVy2_Zee    = genCand->y();
           genVMass2_Zee = genCand->mass();
         }
       }
       // if there is no mother particle that is a Z or a W, mark this by saving the mother particle pT and phi both as -10
       else if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && !foundGenV1) {
         genVPt2_Zee   = -10;
         genVPhi2_Zee  = -10;
         genVy2_Zee    = -10;
         genVMass2_Zee = -10;
       }
     }
   }

   LorentzVector vlep1_Zee(Ele1pt_Zee, Ele1eta_Zee, Ele1phi_Zee, ELE_MASS);
   LorentzVector vlep2_Zee(Ele2pt_Zee, Ele2eta_Zee, Ele2phi_Zee, ELE_MASS);
   LorentzVector vsc1_Zee(Elesc1pt_Zee, Elesc1eta_Zee, Elesc1phi_Zee, ELE_MASS);
   LorentzVector vsc2_Zee(Elesc2pt_Zee, Elesc2eta_Zee, Elesc2phi_Zee, ELE_MASS);

   // Save MET information

   edm::Handle<pat::METCollection> mets;
   iEvent.getByToken(metToken_, mets);
   const pat::MET &met = mets->front();

   edm::Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);
   rawpfMETpx_Zee=0;
   rawpfMETpy_Zee=0;
   // loop on PF candidates to calculate sum of Et's
   for(unsigned int jcand=0; jcand < pfs->size(); jcand++) {
     const pat::PackedCandidate &pf = (*pfs)[jcand];
     rawpfMETpx_Zee -= pf.px();
     rawpfMETpy_Zee -= pf.py();
   }

   vtype1pfMET_Zee.Set(met.px(),met.py());
   vrawpfMET_Zee.Set(rawpfMETpx_Zee,rawpfMETpy_Zee);
   vgenMET_Zee.Set(met.genMET()->px(),met.genMET()->py());

   //
   // SELECTION PROCEDURE:
   //  (1) Find a "tag" electron that passes the tight selection
   //  (2) Find a Supercluster "probe" electron which gives a dilepton mass along with the tag inside the Z-mass window
   //
   Bool_t foundTag=kFALSE, foundProbe=kFALSE;
   Bool_t inECALgap0=kFALSE, inECALgap1=kFALSE;
   // Find tag electron
   if(  fabs(vsc1_Zee.Eta())>=ECAL_GAP_LOW && fabs(vsc1_Zee.Eta())<=ECAL_GAP_HIGH  ) inECALgap0=kTRUE; // check ECAL gap
   if(  !inECALgap0 && fabs(vsc1_Zee.Eta())<=ETA_CUT && vsc1_Zee.Pt()>=PT_CUT && Convert_Zee(isLooseEle1_Zee)  ) foundTag=kTRUE;
   if(  foundTag==kFALSE  ) return; // event not interesting if there is no tag electron
   LorentzVector vTag_Zee(vlep1_Zee.Pt(),vlep1_Zee.Eta(),vlep1_Zee.Phi(),ELE_MASS); lep1_Zee = &vTag_Zee;
   LorentzVector vTagSC_Zee(vsc1_Zee.Pt(),vsc1_Zee.Eta(),vsc1_Zee.Phi(),ELE_MASS);  sc1_Zee  = &vTagSC_Zee;
   // Find probe electron
   if(  fabs(vsc2_Zee.Eta())>=ECAL_GAP_LOW && fabs(vsc2_Zee.Eta())<=ECAL_GAP_HIGH  ) inECALgap1=kTRUE;
   if(  !inECALgap1 && fabs(vsc2_Zee.Eta())<=ETA_CUT && vsc2_Zee.Pt()>=PT_CUT  ) foundProbe=kTRUE;
   if(  foundProbe==kFALSE  ) return; // event not interesting if there is no probe electron
   LorentzVector vProbe_Zee(vlep2_Zee.Pt(),vlep2_Zee.Eta(),vlep2_Zee.Phi(),ELE_MASS); lep2_Zee = &vProbe_Zee;
   LorentzVector vProbeSC_Zee(vsc2_Zee.Pt(),vsc2_Zee.Eta(),vsc2_Zee.Phi(),ELE_MASS);  sc2_Zee  = &vProbeSC_Zee;
   // Mass window
   LorentzVector vdilep_Zee = vTag_Zee + vProbe_Zee; dilep_Zee = &vdilep_Zee;
   TVector2 vdilepPt_Zee(vTag_Zee.px()+vProbe_Zee.px(),vTag_Zee.py()+vProbe_Zee.py());
   if((vdilep_Zee.M()<MASS_LOW) || (vdilep_Zee.M()>MASS_HIGH)) return;

   // Calculate hadronic recoil
   TVector2 uT  = -1*(vdilepPt_Zee+vtype1pfMET_Zee);
   TVector2 vu1 = uT.Proj(vdilepPt_Zee); u1_Zee = vu1.Mod();
   TVector2 vu2 = uT-u1_Zee; u2_Zee = vu2.Mod();

   nsel_Zee    += weight_Zee;
   nsel_Zeevar_Zee += weight_Zee*weight_Zee;
   scale1fb_Zee = weight_Zee;

   //
   // Fill tree
   //
   outTree_Zee->Fill();

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
selectZee::beginJob()
{
   // Create output directory
//   gSystem->mkdir(outputDir,kTRUE);
//   gSystem->mkdir(ntupDir,kTRUE);

  //
  // Set up output ntuple
  //

  outTree_Zee->Branch("nVtx",          &nVtx_Zee,          "nVtx/I");         // number of vertices
  outTree_Zee->Branch("genVPt1",       &genVPt1_Zee,       "genVPt1/F");      // GEN boson pT  (signal MC)  (mother of first lepton )
  outTree_Zee->Branch("genVPhi1",      &genVPhi1_Zee,      "genVPhi1/F");     // GEN boson phi (signal MC)  (mother of first lepton )
  outTree_Zee->Branch("genVy1",        &genVy1_Zee,        "genVy1/F");       // GEN boson rapidity (signal MC)  (mother of first lepton )
  outTree_Zee->Branch("genVMass1",     &genVMass1_Zee,     "genVMass1/F");    // GEN boson mass (signal MC) (mother of first lepton )
  outTree_Zee->Branch("genVPt2",       &genVPt2_Zee,       "genVPt2/F");      // GEN boson pT  (signal MC)  (mother of second lepton)
  outTree_Zee->Branch("genVPhi2",      &genVPhi2_Zee,      "genVPhi2/F");     // GEN boson phi (signal MC)  (mother of second lepton)
  outTree_Zee->Branch("genVy2",        &genVy2_Zee,        "genVy2/F");       // GEN boson rapidity (signal MC)  (mother of second lepton)
  outTree_Zee->Branch("genVMass2",     &genVMass2_Zee,     "genVMass2/F");    // GEN boson mass (signa MC)  (mother of second lepton)
  outTree_Zee->Branch("scale1fb",      &scale1fb_Zee,      "scale1fb/F");     // event weight_Zee per 1/fb (MC)
  outTree_Zee->Branch("vtype1pfMET",   "TVector2",         &vtype1pfMET_Zee);     // type-1 corrected PF MET
  outTree_Zee->Branch("vrawpfMET",     "TVector2",         &vrawpfMET_Zee);       // raw PF MET
  outTree_Zee->Branch("vgenMET",       "TVector2",         &vgenMET_Zee);         // generated MET
  outTree_Zee->Branch("u1",            &u1_Zee,            "u1/F");           // parallel component of recoil
  outTree_Zee->Branch("u2",            &u2_Zee,            "u2/F");           // perpendicular component of recoil
  outTree_Zee->Branch("q1",            &q1_Zee,            "q1/I");           // lepton charge (first lepton )
  outTree_Zee->Branch("q2",            &q2_Zee,            "q2/I");           // lepton charge (second lepton)
  outTree_Zee->Branch("nEvents",       &nEvents_Zee,       "nEvents/I");      // events in MC file
  outTree_Zee->Branch("dilep","ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &dilep_Zee);  // dilepton 4-vector
  outTree_Zee->Branch("lep1", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep1_Zee);   // lepton 4-vector (first lepton )
  outTree_Zee->Branch("lep2", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep2_Zee);   // lepton 4-vector (second lepton)
  outTree_Zee->Branch("sc1",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sc1_Zee);    // supercluster 4-vector (first lepton )
  outTree_Zee->Branch("sc2",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sc2_Zee);    // supercluster 4-vector (second lepton)
  ///// electron specific /////
  outTree_Zee->Branch("pfChIso1",     &pfChIso1_Zee,    "pfChIso1/F");      // PF charged hadron isolation of lepton (first lepton )
  outTree_Zee->Branch("pfGamIso1",    &pfGamIso1_Zee,   "pfGamIso1/F");     // PF photon isolation of lepton         (first lepton )
  outTree_Zee->Branch("pfNeuIso1",    &pfNeuIso1_Zee,   "pfNeuIso1/F");     // PF neutral hadron isolation of lepton (first lepton )
  outTree_Zee->Branch("pfChIso2",     &pfChIso2_Zee,    "pfChIso2/F");      // PF charged hadron isolation of lepton (second lepton)
  outTree_Zee->Branch("pfGamIso2",    &pfGamIso2_Zee,   "pfGamIso2/F");     // PF photon isolation of lepton         (second lepton)
  outTree_Zee->Branch("pfNeuIso2",    &pfNeuIso2_Zee,   "pfNeuIso2/F");     // PF neutral hadron isolation of lepton (second lepton)
  outTree_Zee->Branch("isLooseEle1",  &isLooseEle1_Zee, "isLooseEle1/F");   // loose electron ID (first lepton )
  outTree_Zee->Branch("isTightEle1",  &isTightEle1_Zee, "isTightEle1/F");   // tight electron ID (first lepton )
  outTree_Zee->Branch("isLooseEle2",  &isLooseEle2_Zee, "isLooseEle2/F");   // loose electron ID (second lepton)
  outTree_Zee->Branch("isTightEle2",  &isTightEle2_Zee, "isTightEle2/F");   // tight electron ID (second lepton)
}

// ------------ method called once each job just after ending the event loop  ------------
void 
selectZee::endJob() 
{
   // The only information in the last entry of the tree is the total number of events
   // that were processed (not the same as the total number of events selected).
   nVtx_Zee = -10;
   genVPt1_Zee = -10; genVPhi1_Zee = -10;
   genVPt2_Zee = -10; genVPhi2_Zee = -10;
   scale1fb_Zee = -10;
   vtype1pfMET_Zee.Set(-10.0,-10.0);
   vrawpfMET_Zee.Set(-10.0,-10.0);
   vgenMET_Zee.Set(-10.0,-10.0);
   q1_Zee = -10; q2_Zee = -10;
   nEvents_Zee = dummynEvents_Zee;
   LorentzVector dummyLep(-10, -10, -10, -10);
   lep1_Zee = &dummyLep;
   lep2_Zee = &dummyLep;
   sc1_Zee = &dummyLep;
   sc2_Zee = &dummyLep;
   pfChIso1_Zee = -10; pfGamIso1_Zee = -10; pfNeuIso1_Zee = -10;
   pfChIso2_Zee = -10; pfGamIso2_Zee = -10; pfNeuIso2_Zee = -10;
   isLooseEle1_Zee = kFALSE; isTightEle1_Zee = kFALSE;
   isLooseEle2_Zee = kFALSE; isTightEle2_Zee = kFALSE;
   outTree_Zee->Fill();

   std::cout << nsel_Zee << " +/- " << sqrt(nsel_Zeevar_Zee) << " per 1/fb" << std::endl;
   std::cout << "endJob: nEvents_Zee is " << nEvents_Zee << std::endl;
//   outFile_Zee->Write();
//   outFile_Zee->Close();

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================

  std::cout << "*" << std::endl;
  std::cout << "* SUMMARY" << std::endl;
  std::cout << "*--------------------------------------------------" << std::endl;
  std::cout << "Z -> e e" << std::endl;
  std::cout << " pT > " << PT_CUT << std::endl;
  std::cout << " |eta| < " << ETA_CUT << std::endl;
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "  <> Output saved in " << outFilename_Zee << "/" << std::endl;
  std::cout << std::endl;
}

// ------------ method called when starting to processes a run  ------------

void 
selectZee::beginRun(edm::Run const&, edm::EventSetup const&)
{
  outFile_Zee = new TFile(outFilename_Zee,"RECREATE");
  outTree_Zee = new TTree("Events","Events");
}


// ------------ method called when ending the processing of a run  ------------

void 
selectZee::endRun(edm::Run const&, edm::EventSetup const&)
{
  outFile_Zee->Write();
  outFile_Zee->Close();
}


// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
selectZee::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
selectZee::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
selectZee::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(selectZee);
