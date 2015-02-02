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
// This code is just an adaptation of Kevin Sung's original code used for the 8 TeV analysis:
// https://github.com/jaylawhorn/mitewk/blob/master/Selection/selectZmm.C

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

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

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<pat::METCollection> metToken_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
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

const Double_t MUON_MASS = 0.105658369;

Double_t nsel_Zmm=0;

TString outFilename_Zmm = TString("selectZmm.root");
TFile *outFile_Zmm = new TFile();
TTree *outTree_Zmm = new TTree();

//
// Declare output ntuple variables
//
Int_t  matchGen_Zmm=0, npv_Zmm=0;
Float_t genVPdgID_Zmm=0, genVPt_Zmm=0, genVPhi_Zmm=0, genVy_Zmm=0, genVMass_Zmm=0;
Float_t rawpfmet_Zmm=0, rawpfmetPhi_Zmm=0;
Float_t type1pfmet_Zmm=0, type1pfmetPhi_Zmm=0;
Float_t genmet_Zmm=0, genmetPhi_Zmm=0;
Float_t u1_Zmm=0, u2_Zmm=0;
Int_t q1_Zmm=0, q2_Zmm=0;
Float_t pfChIso1_Zmm=0, pfChIso2_Zmm=0;
LorentzVector *dilep_Zmm=0, *lep1_Zmm=0, *lep2_Zmm=0;
Int_t isLooseMuon1_Zmm=0, isSoftMuon1_Zmm=0, isTightMuon1_Zmm=0;
Int_t isLooseMuon2_Zmm=0, isSoftMuon2_Zmm=0, isTightMuon2_Zmm=0;
Int_t passSingleMuTrigger_Zmm=0, matchTrigObj1_Zmm=0, matchTrigObj2_Zmm=0;

//
// constructors and destructor
//
selectZmm::selectZmm(const edm::ParameterSet& iConfig):
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
   muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
   metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
   pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
   triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
   triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects")))
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

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   // Good vertex requirement
   if (vertices->empty()) return; // Skip the event if no PV found
   const reco::Vertex &PV = vertices->front();
   npv_Zmm = vertices->size();

   Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);
   if(muons->size()<2) return; // Skip the event if there is no possibility of both tag and probe leptons

   edm::Handle<edm::TriggerResults> triggerBits;
   edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   iEvent.getByToken(triggerBits_, triggerBits);
   iEvent.getByToken(triggerObjects_, triggerObjects);
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
   const TString singleMu("HLT_IsoMu20_eta2p1_IterTrk02_v1");
   for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
       if((const TString)names.triggerName(i)!=singleMu) continue;
       passSingleMuTrigger_Zmm = triggerBits->accept(i) ? 1 : 0;
       break;
   }

   //
   // SELECTION PROCEDURE:
   //  (1) Find a good muon that passes the soft ID -> this will be the "tag"
   //  (2) Pair the tag with a probe muon so that that tag+probe mass is
   //      inside the Z window
   //
   Bool_t foundTag = kFALSE;
   for (unsigned int i1=0;i1 < muons->size();i1++) {
     const pat::Muon &tag = (*muons)[i1];

     Int_t isLooseMuon = tag.isLooseMuon() ? 1 : 0;
     if(!isLooseMuon) continue; // lepton selection

     foundTag = kTRUE;

     // Tag lepton information
     LorentzVector vTag(tag.pt(),tag.eta(),tag.phi(),MUON_MASS);

     lep1_Zmm = &vTag;

     q1_Zmm = tag.charge();
     pfChIso1_Zmm = tag.chargedHadronIso();
     isLooseMuon1_Zmm = tag.isLooseMuon() ? 1 : 0;
     isSoftMuon1_Zmm  = tag.isSoftMuon(PV) ? 1 : 0;
     isTightMuon1_Zmm = tag.isTightMuon(PV) ? 1 : 0;

     // Match to a trigger object
     matchTrigObj1_Zmm = 0;
     Int_t nObj = 0, tagTrigObjIdx = triggerObjects->size();
     if(passSingleMuTrigger_Zmm) {
       for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
           matchTrigObj1_Zmm = (sqrt((tag.eta()-obj.eta())*(tag.eta()-obj.eta())+(tag.phi()-obj.phi())*(tag.phi()-obj.phi())) <= 0.2) ? 1 : 0;
           if(matchTrigObj1_Zmm) tagTrigObjIdx = nObj;
           nObj++;
       }
     }

     for (unsigned int i2=0;i2 < muons->size();i2++) {
       if(i2==i1) continue;
       const pat::Muon &probe = (*muons)[i2];

       // Probe lepton information
       q2_Zmm = probe.charge();
       pfChIso2_Zmm = probe.chargedHadronIso();
       isLooseMuon2_Zmm = probe.isLooseMuon() ? 1 : 0;
       isSoftMuon2_Zmm  = probe.isSoftMuon(PV) ? 1 : 0;
       isTightMuon2_Zmm = probe.isTightMuon(PV) ? 1 : 0;

       LorentzVector vProbe(probe.pt(),probe.eta(),probe.phi(),MUON_MASS);

       lep2_Zmm = &vProbe;

       // Match to a trigger object
       matchTrigObj2_Zmm = 0;
       if(passSingleMuTrigger_Zmm) {
         Int_t nProbeObj = 0;
         for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
             if(matchTrigObj1_Zmm && nProbeObj==tagTrigObjIdx) continue; // Do not want to match the tag and the probe to the same object
             matchTrigObj2_Zmm = (sqrt((probe.eta()-obj.eta())*(probe.eta()-obj.eta())+(probe.phi()-obj.phi())*(probe.phi()-obj.phi())) <= 0.2) ? 1 : 0;
             nProbeObj++;
         }
       }

       // Dilepton 4-vector
       LorentzVector vDilep = vTag + vProbe;
       dilep_Zmm = &vDilep;

       // Perform matching of leptons to generator level leptons
       const reco::GenParticle* gen1 = tag.genParticle();
       const reco::GenParticle* gen2 = probe.genParticle();
       Bool_t hasGenMatch = kFALSE;
       if(gen1!=NULL && gen2!=NULL) {
         Int_t id1 = gen1->pdgId(), id2 = gen2->pdgId();
         Float_t eta1 = gen1->eta(), eta2 = gen2->eta();
         Float_t phi1 = gen1->phi(), phi2 = gen2->phi();
         Bool_t match1 = ( fabs(id1)==13 && sqrt((tag.eta()-eta1)*(tag.eta()-eta1)+(tag.phi()-phi1)*(tag.phi()-phi1)) < 0.5 );
         Bool_t match2 = ( fabs(id2)==13 && sqrt((probe.eta()-eta2)*(probe.eta()-eta2)+(probe.phi()-phi2)*(probe.phi()-phi2)) < 0.5 );
         if(match1 && match2) hasGenMatch = kTRUE;
       }
       matchGen_Zmm  = hasGenMatch ? 1 : 0;

       // Type-1 corrected PF MET
       edm::Handle<pat::METCollection> mets;
       iEvent.getByToken(metToken_, mets);

       const pat::MET &met = mets->front();
       TVector2 vtype1pfMET_Zmm = TVector2(met.px(),met.py());
       type1pfmet_Zmm    = vtype1pfMET_Zmm.Mod();
       type1pfmetPhi_Zmm = vtype1pfMET_Zmm.Phi();

       // Generator level MET
       TVector2 vgenMET_Zmm = TVector2(met.genMET()->px(),met.genMET()->py());
       genmet_Zmm    = vgenMET_Zmm.Mod();
       genmetPhi_Zmm = vgenMET_Zmm.Phi();

       // Raw PF MET
       edm::Handle<pat::PackedCandidateCollection> pfs;
       iEvent.getByToken(pfToken_, pfs);
       Float_t rawpfMETpx=0, rawpfMETpy=0;

       for(unsigned int jcand=0; jcand < pfs->size(); jcand++) {
         const pat::PackedCandidate &pf = (*pfs)[jcand];
         rawpfMETpx -= pf.px();
         rawpfMETpy -= pf.py();
       }

       TVector2 vrawpfMET_Zmm = TVector2(rawpfMETpx,rawpfMETpy);
       rawpfmet_Zmm    = vrawpfMET_Zmm.Mod();
        rawpfmetPhi_Zmm = vrawpfMET_Zmm.Phi();

       // Generator level lepton information and hadronic recoil
       genVPdgID_Zmm = 0;
       genVPt_Zmm    = 0;
       genVPhi_Zmm   = 0;
       genVy_Zmm     = 0;
       genVMass_Zmm  = 0;
       u1_Zmm        = 0;
       u2_Zmm        = 0;
       if(gen1!=NULL) {
         const reco::Candidate* mother = gen1->mother(0);
         genVPdgID_Zmm = mother->pdgId();
         genVPt_Zmm    = mother->pt();
         genVPhi_Zmm   = mother->phi();
         genVy_Zmm     = mother->y();
         genVMass_Zmm  = mother->mass();
         TVector2 vZPt((vDilep.Pt())*cos(vDilep.Phi()),(vDilep.Pt())*sin(vDilep.Phi()));
         TVector2 vMet(type1pfmet_Zmm*cos(type1pfmetPhi_Zmm), type1pfmet_Zmm*sin(type1pfmetPhi_Zmm));
         TVector2 vU = -1.0*(vMet+vZPt);
         u1_Zmm = ((vDilep.Px())*(vU.Px()) + (vDilep.Py())*(vU.Py()))/(vDilep.Pt()); // u1 = (pT . u)/|pT|
         u2_Zmm = ((vDilep.Px())*(vU.Py()) - (vDilep.Py())*(vU.Px()))/(vDilep.Pt()); // u2 = (pT x u)/|pT|
       }
       // Count total number of events selected
       nsel_Zmm++;
       // Fill tree
       outTree_Zmm->Fill();
       return;
     } // End of i2 loop
   } // End of i1 loop
   if(!foundTag) return; // Skip event if there is no tag lepton

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

  //
  // Set up output ntuple
  //

  outFile_Zmm = new TFile(outFilename_Zmm,"RECREATE");
  outTree_Zmm = new TTree("Events","Events");

  outTree_Zmm->Branch("matchGen",      &matchGen_Zmm,      "matchGen/I");      // event has both leptons matched to MC Z->ll
  outTree_Zmm->Branch("npv",           &npv_Zmm,           "npv/I");           // number of vertices
  outTree_Zmm->Branch("genVPt",        &genVPt_Zmm,        "genVPt/F");        // GEN boson pT (signal MC)
  outTree_Zmm->Branch("genVPhi",       &genVPhi_Zmm,       "genVPhi/F");       // GEN boson phi (signal MC)
  outTree_Zmm->Branch("genVy",         &genVy_Zmm,         "genVy/F");         // GEN boson rapidity (signal MC)
  outTree_Zmm->Branch("genVMass",      &genVMass_Zmm,      "genVMass/F");      // GEN boson mass (signal MC)
  outTree_Zmm->Branch("rawpfmet",      &rawpfmet_Zmm,      "rawpfmet/F");      // Raw PF MET
  outTree_Zmm->Branch("rawpfmetPhi",   &rawpfmetPhi_Zmm,   "rawpfmetPhi/F");   // Raw PF MET phi
  outTree_Zmm->Branch("type1pfmet",    &type1pfmet_Zmm,    "type1pfmet/F");    // Type-1 corrected PF MET
  outTree_Zmm->Branch("type1pfmetPhi", &type1pfmetPhi_Zmm, "type1pfmetPhi/F"); // Type-1 corrected PF MET phi
  outTree_Zmm->Branch("genmet",        &genmet_Zmm,        "genmet/F");        // Generator level MET
  outTree_Zmm->Branch("genmetPhi",     &genmetPhi_Zmm,     "genmetPhi/F");     // Generator level MET phi
  outTree_Zmm->Branch("u1",            &u1_Zmm,            "u1/F");            // parallel component of recoil
  outTree_Zmm->Branch("u2",            &u2_Zmm,            "u2/F");            // perpendicular component of recoil
  outTree_Zmm->Branch("q1",            &q1_Zmm,            "q1/I");            // tag lepton charge
  outTree_Zmm->Branch("q2",            &q2_Zmm,            "q2/I");            // probe lepton charge
  outTree_Zmm->Branch("pfChIso1",      &pfChIso1_Zmm,      "pfChIso1/F");      // tag lepton charged hadron isolation
  outTree_Zmm->Branch("pfChIso2",      &pfChIso2_Zmm,      "pfChIso2/F");      // probe lepton charged hadron isolation
  outTree_Zmm->Branch("dilep","ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &dilep_Zmm); // dilepton 4-vector
  outTree_Zmm->Branch("lep1", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep1_Zmm);  // lepton 4-vector (tag lepton)
  outTree_Zmm->Branch("lep2", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep2_Zmm);  // lepton 4-vector (probe lepton)
  outTree_Zmm->Branch("isLooseMuon1",  &isLooseMuon1_Zmm, "isLooseMuon1/I");  // tag lepton loose muon ID
  outTree_Zmm->Branch("isSoftMuon1",   &isSoftMuon1_Zmm,  "isSoftMuon1/I");   // tag lepton soft muon ID
  outTree_Zmm->Branch("isTightMuon1",  &isTightMuon1_Zmm, "isTightMuon1/I");  // tag lepton tight muon ID
  outTree_Zmm->Branch("isLooseMuon2",  &isLooseMuon2_Zmm, "isLooseMuon2/I");  // probe lepton loose muon ID
  outTree_Zmm->Branch("isSoftMuon2",   &isSoftMuon2_Zmm,  "isSoftMuon2/I");   // probe lepton soft muon ID
  outTree_Zmm->Branch("isTightMuon2",  &isTightMuon2_Zmm, "isTightMuon2/I");  // probe lepton tight muon ID
  outTree_Zmm->Branch("passSingleMuTrigger", &passSingleMuTrigger_Zmm, "passSingleMuTrigger/I"); // single muon trigger
  outTree_Zmm->Branch("matchTrigObj1", &matchTrigObj1_Zmm,  "matchTrigObj1/I"); // tag lepton match to trigger object
  outTree_Zmm->Branch("matchTrigObj2", &matchTrigObj2_Zmm,  "matchTrigObj2/I"); // probe lepton match to trigger object
}

// ------------ method called once each job just after ending the event loop  ------------
void 
selectZmm::endJob() 
{
   // Save tree in output file
   outFile_Zmm->Write();
   outFile_Zmm->Close();

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================

  std::cout << "*" << std::endl;
  std::cout << "* SUMMARY" << std::endl;
  std::cout << "*--------------------------------------------------" << std::endl;
  std::cout << "Z -> m m" << std::endl;
  std::cout << nsel_Zmm << " +/- " << sqrt(nsel_Zmm) << " per 1/fb" << std::endl;
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "  <> Output saved in " << outFilename_Zmm << "/" << std::endl;
  std::cout << std::endl;

}

// ------------ method called when starting to processes a run  ------------
/*
void 
selectZmm::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
selectZmm::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

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
