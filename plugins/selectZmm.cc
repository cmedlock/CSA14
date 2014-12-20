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

class selectZmm : public edm::EDAnalyzer {
   public:
      explicit selectZmm(const edm::ParameterSet&);
      ~selectZmm();

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
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<pat::METCollection> metToken_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
};

//
// constants, enums and typedefs
//

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

// Convert int to binary value
bool Convert_Zmm(unsigned int val,bool print=kFALSE)
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
const Double_t ETA_CUT   = 2.4;
const Double_t MUON_MASS  = 0.105658369;

Double_t nsel_Zmm=0, nsel_Zmmvar_Zmm=0;

TString outFilename_Zmm = TString("Zmm_select.root");
TFile *outFile_Zmm = new TFile();
TTree *outTree_Zmm = new TTree();

//
// Declare output ntuple variables
//
// Initialize everything to -10 so that it is obvious if an object (such
// as a Z or W generator level mother particle) doesn't exist
//
Int_t   nVtx_Zmm;
Float_t genVPt1_Zmm=-10, genVPhi1_Zmm=-10, genVy1_Zmm=-10, genVMass1_Zmm=-10;
Float_t genVPt2_Zmm=-10, genVPhi2_Zmm=-10, genVy2_Zmm=-10, genVMass2_Zmm=-10;
Float_t scale1fb_Zmm;
Float_t rawpfMETpx_Zmm, rawpfMETpy_Zmm; // Will be used in vrawpfMET_Zmm
TVector2 vtype1pfMET_Zmm, vrawpfMET_Zmm, vgenMET_Zmm;
Float_t u1_Zmm=-10, u2_Zmm=-10;
Int_t   q1_Zmm=-10, q2_Zmm=-10;
Int_t   dummynEvents_Zmm=0, nEvents_Zmm=0;
Float_t Mu1pt_Zmm=-10, Mu1eta_Zmm=-10, Mu1phi_Zmm=-10; // Will be used in lep1_Zmm
Float_t Mu2pt_Zmm=-10, Mu2eta_Zmm=-10, Mu2phi_Zmm=-10; // Will be used in lep2_Zmm
Float_t Musta1pt_Zmm=-10, Musta1eta_Zmm=-10, Musta1phi_Zmm=-10; // Will be used in sta1_Zmm
Float_t Musta2pt_Zmm=-10, Musta2eta_Zmm=-10, Musta2phi_Zmm=-10; // Will be used in sta2_Zmm
LorentzVector *dilep_Zmm=0;
LorentzVector *lep1_Zmm=0, *lep2_Zmm=0;
LorentzVector *sta1_Zmm=0, *sta2_Zmm=0;
///// muon specific /////
Float_t pfChIso1_Zmm=-10, pfGamIso1_Zmm=-10, pfNeuIso1_Zmm=-10;
Float_t pfChIso2_Zmm=-10, pfGamIso2_Zmm=-10, pfNeuIso2_Zmm=-10;
Float_t isLooseMuon1_Zmm=-10, isSoftMuon1_Zmm=-10, isTightMuon1_Zmm=-10;
Float_t isLooseMuon2_Zmm=-10, isSoftMuon2_Zmm=-10, isTightMuon2_Zmm=-10;

// Compute MC event weight_Zmm_sel per 1/fb
Double_t weight_Zmm = 1;
//const Double_t xsec = 1;
//if(xsec>0) weight_Zmm = 1000.*xsec/(Double_t)eventTree->GetEntries();

//
// constructors and destructor
//
selectZmm::selectZmm(const edm::ParameterSet& iConfig):
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
   muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
   metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
   pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands")))
{
   //now do what ever initialization is needed
}


selectZmm::~selectZmm()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
selectZmm::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   isLooseMuon1_Zmm=-10; isSoftMuon1_Zmm=-10; isTightMuon1_Zmm=-10;
   isLooseMuon2_Zmm=-10; isSoftMuon2_Zmm=-10; isTightMuon2_Zmm=-10;

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   // good vertex requirement
   if (vertices->empty()) return; // skip the event if no PV found
   const reco::Vertex &PV = vertices->front();
   nVtx_Zmm = vertices->size();

   Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);
   if(muons->size()==0) return;   // skip the event if no muons found

   dummynEvents_Zmm++;

   unsigned int idxMu1=-1, idxMu2=-1;
   float maxMupt=-1, secondmaxMupt=-1;
   // find the muon with the highest pT
   for (unsigned int jMuon=0;jMuon < muons->size();jMuon++) {
     const pat::Muon &mu = (*muons)[jMuon];
     if     (maxMupt==-1)                    { maxMupt=mu.pt(); idxMu1=jMuon; }
     else if(maxMupt!=-1 && mu.pt()>maxMupt) { maxMupt=mu.pt(); idxMu1=jMuon; }
   }
   const pat::Muon &Mu1 = (*muons)[idxMu1];
   q1_Zmm           = Mu1.charge();
   Mu1pt_Zmm        = Mu1.pt();
   Mu1eta_Zmm       = Mu1.eta();
   Mu1phi_Zmm       = Mu1.phi();
   Musta1pt_Zmm     = Mu1.superCluster()->energy()*(Mu1.pt()/Mu1.p());
   Musta1eta_Zmm    = Mu1.superCluster()->eta();
   Musta1phi_Zmm    = Mu1.superCluster()->phi();
   pfChIso1_Zmm     = Mu1.chargedHadronIso();
   pfGamIso1_Zmm    = Mu1.photonIso();
   pfNeuIso1_Zmm    = Mu1.neutralHadronIso();
   isLooseMuon1_Zmm = Mu1.isLooseMuon();
   isSoftMuon1_Zmm  = Mu1.isSoftMuon(PV);
   isTightMuon1_Zmm = Mu1.isTightMuon(PV);
   // get mother particle information
   const reco::GenParticle* gen = Mu1.genLepton();
   Bool_t foundGenV0 = kFALSE; // some events do not have a generated Z or W, this protects against a seg fault
   if(gen!=NULL) {
     const reco::Candidate* genCand = gen;
     while(genCand!=NULL && genCand->numberOfMothers()==1) {
       genCand = genCand->mother(0);
       if( (  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && !foundGenV0 ) foundGenV0=kTRUE;
       genVPt1_Zmm   = genCand->pt();
       genVPhi1_Zmm  = genCand->phi();
       genVy1_Zmm    = genCand->y();
       genVMass1_Zmm = genCand->mass();
     }
     // mother particle should always be either a Z or a W
     // if it isn't, trace back to the first daughter particle that is either a Z or a W
     if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && foundGenV0) {
       while( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) ) {
         genCand = genCand->daughter(0);
         genVPt1_Zmm   = genCand->pt();
         genVPhi1_Zmm  = genCand->phi();
         genVy1_Zmm    = genCand->y();
         genVMass1_Zmm = genCand->mass();
       }
     }
     // if there is no mother particle that is a Z or a W, mark this by saving the mother particle pT and phi both as -10
     else if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && !foundGenV0) {
       genVPt1_Zmm   = -10;
       genVPhi1_Zmm  = -10;
       genVy1_Zmm    = -10;
       genVMass1_Zmm = -10;
     }
   }
   // if there is > 1 muon, find the one with the second highest pT
   if(muons->size() > 1) {
     for (unsigned int kMuon=0; kMuon < muons->size();kMuon++) {
       const pat::Muon &mu = (*muons)[kMuon];
       if(kMuon==idxMu1)                                                    continue;
       else if(kMuon!=idxMu1 && secondmaxMupt==-1)                          { secondmaxMupt=mu.pt(); idxMu2=kMuon; }
       else if(kMuon!=idxMu1 && secondmaxMupt!=-1 && mu.pt()>secondmaxMupt) { secondmaxMupt=mu.pt(); idxMu2=kMuon; }
     }
     const pat::Muon &Mu2 = (*muons)[idxMu2];
     q2_Zmm           = Mu2.charge();
     Mu2pt_Zmm        = Mu2.pt();
     Mu2eta_Zmm       = Mu2.eta();
     Mu2phi_Zmm       = Mu2.phi();
     Musta2pt_Zmm     = Mu2.superCluster()->energy()*(Mu2.pt()/Mu2.p());
     Musta2eta_Zmm    = Mu2.superCluster()->eta();
     Musta2phi_Zmm    = Mu2.superCluster()->phi();
     pfChIso2_Zmm     = Mu2.chargedHadronIso();
     pfGamIso2_Zmm    = Mu2.photonIso();
     pfNeuIso2_Zmm    = Mu2.neutralHadronIso();
     isLooseMuon2_Zmm = Mu2.isLooseMuon();
     isSoftMuon2_Zmm  = Mu2.isSoftMuon(PV);
     isTightMuon2_Zmm = Mu2.isTightMuon(PV);
     // get mother particle information
     const reco::GenParticle* gen = Mu2.genLepton();
     Bool_t foundGenV1 = kFALSE; // some events do not have a generated Z or W, this protects against a seg fault
     if(gen!=NULL) {
       const reco::Candidate* genCand = gen;
       while(genCand!=NULL && genCand->numberOfMothers()==1) {
         genCand = genCand->mother(0);
         if( (  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && !foundGenV1 ) foundGenV1=kTRUE;
         genVPt2_Zmm   = genCand->pt();
         genVPhi2_Zmm  = genCand->phi();
         genVy2_Zmm    = genCand->y();
         genVMass2_Zmm = genCand->mass();
       }
       // mother particle should always be either a Z or a W
       // if it isn't, trace back to the first daughter particle that is either a Z or a W
       if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && foundGenV1) {
         while( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) ) {
           genCand = genCand->daughter(0);
           genVPt2_Zmm   = genCand->pt();
           genVPhi2_Zmm  = genCand->phi();
           genVy2_Zmm    = genCand->y();
           genVMass2_Zmm = genCand->mass();
         }
       }
       // if there is no mother particle that is a Z or a W, mark this by saving the mother particle pT and phi both as -10
       else if( !(  fabs(genCand->pdgId())==23 || fabs(genCand->pdgId())==24  ) && !foundGenV1) {
         genVPt2_Zmm   = -10;
         genVPhi2_Zmm  = -10;
         genVy2_Zmm    = -10;
         genVMass2_Zmm = -10;
       }
     }
   }

   LorentzVector vlep1_Zmm(Mu1pt_Zmm, Mu1eta_Zmm, Mu1phi_Zmm, MUON_MASS);
   LorentzVector vlep2_Zmm(Mu2pt_Zmm, Mu2eta_Zmm, Mu2phi_Zmm, MUON_MASS);
   LorentzVector vsta1_Zmm(Musta1pt_Zmm, Musta1eta_Zmm, Musta1phi_Zmm, MUON_MASS);
   LorentzVector vsta2_Zmm(Musta2pt_Zmm, Musta2eta_Zmm, Musta2phi_Zmm, MUON_MASS);

   // Save MET information

   edm::Handle<pat::METCollection> mets;
   iEvent.getByToken(metToken_, mets);
   const pat::MET &met = mets->front();

   edm::Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);
   rawpfMETpx_Zmm=0;
   rawpfMETpy_Zmm=0;
   // loop on PF candidates to calculate sum of Et's
   for(unsigned int jcand=0; jcand < pfs->size(); jcand++) {
     const pat::PackedCandidate &pf = (*pfs)[jcand];
     rawpfMETpx_Zmm -= pf.px();
     rawpfMETpy_Zmm -= pf.py();
   }

   vtype1pfMET_Zmm.Set(met.px(),met.py());
   vrawpfMET_Zmm.Set(rawpfMETpx_Zmm,rawpfMETpy_Zmm);
   vgenMET_Zmm.Set(met.genMET()->px(),met.genMET()->py());

   //
   // SELECTION PROCEDURE:
   //  (1) Find a "tag" muon that passes the tight selection
   //  (2) Find a stand-alone "probe" muon which gives a dilepton mass along with the tag inside the Z-mass window
   //
   Bool_t foundTag=kFALSE, foundProbe=kFALSE;
   // Find tag muon
   if(  fabs(vsta1_Zmm.Eta())<=ETA_CUT && vsta1_Zmm.Pt()>=PT_CUT && Convert_Zmm(isLooseMuon1_Zmm)  ) foundTag=kTRUE;
   if(  foundTag==kFALSE  ) return; // event not interesting if there is no tag muon
   LorentzVector vTag_Zmm(vlep1_Zmm.Pt(),vlep1_Zmm.Eta(),vlep1_Zmm.Phi(),MUON_MASS); lep1_Zmm = &vTag_Zmm;
   LorentzVector vTagSTA_Zmm(vsta1_Zmm.Pt(),vsta1_Zmm.Eta(),vsta1_Zmm.Phi(),MUON_MASS); sta1_Zmm  = &vTagSTA_Zmm;
   // Find probe muon
   if(  fabs(vsta2_Zmm.Eta())<=ETA_CUT && vsta2_Zmm.Pt()>=PT_CUT  ) foundProbe=kTRUE;
   if(  foundProbe==kFALSE  ) return; // event not interesting if there is no probe muon
   LorentzVector vProbe_Zmm(vlep2_Zmm.Pt(),vlep2_Zmm.Eta(),vlep2_Zmm.Phi(),MUON_MASS); lep2_Zmm = &vProbe_Zmm;
   LorentzVector vProbeSTA_Zmm(vsta2_Zmm.Pt(),vsta2_Zmm.Eta(),vsta2_Zmm.Phi(),MUON_MASS); sta2_Zmm  = &vProbeSTA_Zmm;
   // Mass window
   LorentzVector vdilep_Zmm = vTag_Zmm + vProbe_Zmm; dilep_Zmm = &vdilep_Zmm;
   TVector2 vdilepPt_Zmm(vTag_Zmm.px()+vProbe_Zmm.px(),vTag_Zmm.py()+vProbe_Zmm.py());
   if((vdilep_Zmm.M()<MASS_LOW) || (vdilep_Zmm.M()>MASS_HIGH)) return;

   // Calculate hadronic recoil
   TVector2 uT  = -1*(vdilepPt_Zmm+vtype1pfMET_Zmm);
   TVector2 vu1 = uT.Proj(vdilepPt_Zmm); u1_Zmm = vu1.Mod();
   TVector2 vu2 = uT-u1_Zmm; u2_Zmm = vu2.Mod();

   nsel_Zmm    += weight_Zmm;
   nsel_Zmmvar_Zmm += weight_Zmm*weight_Zmm;
   scale1fb_Zmm = weight_Zmm;

   //
   // Fill tree
   //
   outTree_Zmm->Fill();

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
selectZmm::beginJob()
{
   // Create output directory
//   gSystem->mkdir(outputDir,kTRUE);
//   gSystem->mkdir(ntupDir,kTRUE);

  //
  // Set up output ntuple
  //

  outTree_Zmm->Branch("nVtx",          &nVtx_Zmm,          "nVtx/I");         // number of vertices
  outTree_Zmm->Branch("genVPt1",       &genVPt1_Zmm,       "genVPt1/F");      // GEN boson pT  (signal MC)  (mother of tag lepton )
  outTree_Zmm->Branch("genVPhi1",      &genVPhi1_Zmm,      "genVPhi1/F");     // GEN boson phi (signal MC)  (mother of tag lepton )
  outTree_Zmm->Branch("genVy1",        &genVy1_Zmm,        "genVy1/F");       // GEN boson rapidity (signal MC)  (mother of tag lepton )
  outTree_Zmm->Branch("genVMass1",     &genVMass1_Zmm,     "genVMass1/F");    // GEN boson mass (signal MC) (mother of tag lepton )
  outTree_Zmm->Branch("genVPt2",       &genVPt2_Zmm,       "genVPt2/F");      // GEN boson pT  (signal MC)  (mother of probe lepton)
  outTree_Zmm->Branch("genVPhi2",      &genVPhi2_Zmm,      "genVPhi2/F");     // GEN boson phi (signal MC)  (mother of probe lepton)
  outTree_Zmm->Branch("genVy2",        &genVy2_Zmm,        "genVy2/F");       // GEN boson rapidity (signal MC)  (mother of probe lepton)
  outTree_Zmm->Branch("genVMass2",     &genVMass2_Zmm,     "genVMass2/F");    // GEN boson mass (signa MC)  (mother of probe lepton)
  outTree_Zmm->Branch("scale1fb",      &scale1fb_Zmm,      "scale1fb/F");     // event weight_Zmm per 1/fb (MC)
  outTree_Zmm->Branch("vtype1pfMET",   "TVector2",         &vtype1pfMET_Zmm);     // type-1 corrected PF MET
  outTree_Zmm->Branch("vrawpfMET",     "TVector2",         &vrawpfMET_Zmm);       // raw PF MET
  outTree_Zmm->Branch("vgenMET",       "TVector2",         &vgenMET_Zmm);         // generated MET
  outTree_Zmm->Branch("u1",            &u1_Zmm,            "u1/F");           // parallel component of recoil
  outTree_Zmm->Branch("u2",            &u2_Zmm,            "u2/F");           // perpendicular component of recoil
  outTree_Zmm->Branch("q1",            &q1_Zmm,            "q1/I");           // lepton charge (tag lepton )
  outTree_Zmm->Branch("q2",            &q2_Zmm,            "q2/I");           // lepton charge (probe lepton)
  outTree_Zmm->Branch("nEvents",       &nEvents_Zmm,       "nEvents/I");      // events in MC file
  outTree_Zmm->Branch("dilep","ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &dilep_Zmm);  // dilepton 4-vector
  outTree_Zmm->Branch("lep1", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep1_Zmm);   // lepton 4-vector (tag lepton )
  outTree_Zmm->Branch("lep2", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep2_Zmm);   // lepton 4-vector (probe lepton)
  outTree_Zmm->Branch("sta1", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sta1_Zmm);    // stand-alone 4-vector (tag lepton )
  outTree_Zmm->Branch("sta2", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sta2_Zmm);    // stand-alone 4-vector (probe lepton)
  ///// electron specific /////
  outTree_Zmm->Branch("pfChIso1",     &pfChIso1_Zmm,    "pfChIso1/F");      // PF charged hadron isolation of lepton (tag lepton )
  outTree_Zmm->Branch("pfGamIso1",    &pfGamIso1_Zmm,   "pfGamIso1/F");     // PF photon isolation of lepton         (tag lepton )
  outTree_Zmm->Branch("pfNeuIso1",    &pfNeuIso1_Zmm,   "pfNeuIso1/F");     // PF neutral hadron isolation of lepton (tag lepton )
  outTree_Zmm->Branch("pfChIso2",     &pfChIso2_Zmm,    "pfChIso2/F");      // PF charged hadron isolation of lepton (probe lepton)
  outTree_Zmm->Branch("pfGamIso2",    &pfGamIso2_Zmm,   "pfGamIso1/F");     // PF photon isolation of lepton         (probe lepton)
  outTree_Zmm->Branch("pfNeuIso2",    &pfNeuIso2_Zmm,   "pfNeuIso2/F");     // PF neutral hadron isolation of lepton (probe lepton)
  outTree_Zmm->Branch("isLooseMuon1", &isLooseMuon1_Zmm, "isLooseMuon1/F");   // loose muon ID (tag lepton )
  outTree_Zmm->Branch("isSoftMuon1",  &isSoftMuon1_Zmm,  "isSoftMuon1/F");    // soft muon ID (tag lepton )
  outTree_Zmm->Branch("isTightMuon1", &isTightMuon1_Zmm, "isTightMuon1/F");   // tight muon ID (tag lepton )
  outTree_Zmm->Branch("isLooseMuon2", &isLooseMuon2_Zmm, "isLooseMuon2/F");   // loose muon ID (probe lepton)
  outTree_Zmm->Branch("isSoftMuon2",  &isSoftMuon2_Zmm,  "isSoftMuon2/F");    // soft muon ID (probe lepton )
  outTree_Zmm->Branch("isTightMuon2", &isTightMuon2_Zmm, "isTightMuon2/F");   // tight muon ID (probe lepton)
}

// ------------ method called once each job just after ending the event loop  ------------
void 
selectZmm::endJob() 
{
   // The only information in the last entry of the tree is the total number of events
   // that were processed (not the same as the total number of events selected).
   nVtx_Zmm = -10;
   genVPt1_Zmm = -10; genVPhi1_Zmm = -10;
   genVPt2_Zmm = -10; genVPhi2_Zmm = -10;
   scale1fb_Zmm = -10;
   vtype1pfMET_Zmm.Set(-10.0,-10.0);
   vrawpfMET_Zmm.Set(-10.0,-10.0);
   vgenMET_Zmm.Set(-10.0,-10.0);
   q1_Zmm = -10; q2_Zmm = -10;
   nEvents_Zmm = dummynEvents_Zmm;
   LorentzVector dummyLep(-10, -10, -10, -10);
   lep1_Zmm = &dummyLep;
   lep2_Zmm = &dummyLep;
   sta1_Zmm = &dummyLep;
   sta2_Zmm = &dummyLep;
   pfChIso1_Zmm = -10; pfGamIso1_Zmm = -10; pfNeuIso1_Zmm = -10;
   pfChIso2_Zmm = -10; pfGamIso2_Zmm = -10; pfNeuIso2_Zmm = -10;
   isLooseMuon1_Zmm = kFALSE; isTightMuon1_Zmm = kFALSE;
   isLooseMuon2_Zmm = kFALSE; isTightMuon2_Zmm = kFALSE;
   outTree_Zmm->Fill();

   std::cout << nsel_Zmm << " +/- " << sqrt(nsel_Zmmvar_Zmm) << " per 1/fb" << std::endl;
   std::cout << "endJob: nEvents_Zmm is " << nEvents_Zmm << std::endl;
//   outFile_Zmm->Write();
//   outFile_Zmm->Close();

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================

  std::cout << "*" << std::endl;
  std::cout << "* SUMMARY" << std::endl;
  std::cout << "*--------------------------------------------------" << std::endl;
  std::cout << "Z -> m m" << std::endl;
  std::cout << " pT > " << PT_CUT << std::endl;
  std::cout << " |eta| < " << ETA_CUT << std::endl;
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "  <> Output saved in " << outFilename_Zmm << "/" << std::endl;
  std::cout << std::endl;
}

// ------------ method called when starting to processes a run  ------------

void 
selectZmm::beginRun(edm::Run const&, edm::EventSetup const&)
{
  outFile_Zmm = new TFile(outFilename_Zmm,"RECREATE");
  outTree_Zmm = new TTree("Events","Events");
}


// ------------ method called when ending the processing of a run  ------------

void 
selectZmm::endRun(edm::Run const&, edm::EventSetup const&)
{
  outFile_Zmm->Write();
  outFile_Zmm->Close();
}


// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
selectZmm::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
selectZmm::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
selectZmm::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(selectZmm);
