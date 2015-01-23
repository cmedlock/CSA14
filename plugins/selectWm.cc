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
// https://github.com/jaylawhorn/mitewk/blob/master/Selection/selectWm.C

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

Double_t nsel_Wm=0, eventCounter_Wm=0;

TString outFilename_Wm = TString("selectWm.root");
TFile *outFile_Wm = new TFile();
TTree *outTree_Wm = new TTree();

//
// Declare output ntuple variables
//
Int_t  npv_Wm=0, nEvents_Wm=0;
Float_t genVPdgID_Wm=0, genVPt_Wm=0, genVPhi_Wm=0, genVy_Wm=0, genVMass_Wm=0;
Float_t genLepPdgID_Wm=0, genLepPt_Wm=0, genLepPhi_Wm=0;
Float_t rawpfmet_Wm=0, rawpfmetPhi_Wm=0;
Float_t type1pfmet_Wm=0, type1pfmetPhi_Wm=0;
Float_t genmet_Wm=0, genmetPhi_Wm=0;
Float_t mt_Wm=0, u1_Wm=0, u2_Wm=0;
Int_t q_Wm=0;
LorentzVector *lep_Wm=0;
///// muon specific /////
Float_t pfChIso_Wm=0, pfGamIso_Wm=0, pfNeuIso_Wm=0;
Int_t isLooseMuon_Wm=0, isSoftMuon_Wm=0, isTightMuon_Wm=0;

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

   // Count total number of MC events processed
   eventCounter_Wm++;

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   // good vertex requirement
   if (vertices->empty()) return; // skip the event if no PV found
   const reco::Vertex &PV = vertices->front();
   npv_Wm = vertices->size();

   Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);
   if(muons->size()==0) return;   // skip the event if no muons found

   //
   // SELECTION PROCEDURE:
   //  (1) Look for 1 good muon matched to trigger
   //  (2) Reject event if another muon is present passing looser cuts
   //
   Int_t  nLooseLep=0;
   Bool_t passSel=kFALSE;
   Int_t  goodMuonIdx=0;
   for (unsigned int jMuon=0;jMuon < muons->size();jMuon++) {
     const pat::Muon &mu = (*muons)[jMuon];

     isLooseMuon_Wm = mu.isLooseMuon();
     isTightMuon_Wm = mu.isTightMuon(PV);

     if(fabs(mu.eta()) > 2.4) continue;    // lepton |eta| cut
     if(mu.pt()        < 10 ) continue;    // lepton pT cut
     if(isLooseMuon_Wm)   nLooseLep++; // loose lepton selection
     if(nLooseLep>1) { // extra lepton veto
       passSel=kFALSE;
       break;
     }

     Double_t iso = mu.chargedHadronIso() + TMath::Max(mu.neutralHadronIso() + mu.photonIso() - 0.5*mu.puChargedHadronIso(),Double_t(0));

     if(fabs(mu.eta()) > ETA_CUT) continue; // lepton |eta| cut
     if(mu.pt()        < PT_CUT ) continue; // lepton pT cut
     if(!(isTightMuon_Wm && iso <= 0.12*mu.pt())) continue; // lepton selection

     passSel = kTRUE;
     goodMuonIdx = jMuon;
   }
   if(passSel==kFALSE) return;

   // Count total number of events selected
   nsel_Wm++;

   const pat::Muon &goodMuon = (*muons)[goodMuonIdx];
   // Lepton information
   q_Wm           = goodMuon.charge();
   pfChIso_Wm     = goodMuon.chargedHadronIso();
   pfGamIso_Wm    = goodMuon.photonIso();
   pfNeuIso_Wm    = goodMuon.neutralHadronIso();
   isLooseMuon_Wm = goodMuon.isLooseMuon() ? 1 : 0;
   isSoftMuon_Wm  = goodMuon.isSoftMuon(PV) ? 1 : 0;
   isTightMuon_Wm = goodMuon.isTightMuon(PV) ? 1 : 0;

   // Lepton and supercluster 4-vectors
   LorentzVector vLep(goodMuon.pt(),goodMuon.eta(),goodMuon.phi(),MUON_MASS);

   lep_Wm = &vLep;

   // Type-1 corrected PF MET (default)
   edm::Handle<pat::METCollection> mets;
   iEvent.getByToken(metToken_, mets);

   const pat::MET &met = mets->front();
   TVector2 vtype1pfMET_Wm = TVector2(met.px(),met.py());
   type1pfmet_Wm    = vtype1pfMET_Wm.Mod();
   type1pfmetPhi_Wm = vtype1pfMET_Wm.Phi();

   // Generator level MET
   TVector2 vgenMET_Wm = TVector2(met.genMET()->px(),met.genMET()->py());
   genmet_Wm    = vgenMET_Wm.Mod();
   genmetPhi_Wm = vgenMET_Wm.Phi();

   // Raw PF MET
   edm::Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);

   Float_t rawpfMETpx=0, rawpfMETpy=0;
   for(unsigned int jcand=0; jcand < pfs->size(); jcand++) {
     const pat::PackedCandidate &pf = (*pfs)[jcand];
     rawpfMETpx -= pf.px();
     rawpfMETpy -= pf.py();
   }

   TVector2 vrawpfMET_Wm = TVector2(rawpfMETpx,rawpfMETpy);
   rawpfmet_Wm    = vrawpfMET_Wm.Mod();
   rawpfmetPhi_Wm = vrawpfMET_Wm.Phi();

   // Lepton transverse mass
   mt_Wm       = sqrt(2.0*vLep.Pt()*vtype1pfMET_Wm.Mod()*(1.0-cos(vLep.Phi()-vtype1pfMET_Wm.Phi())));

   // Generator level lepton information and hadronic recoil
   const reco::GenParticle* genLep = goodMuon.genLepton();
   if(genLep!=NULL) {
     genLepPdgID_Wm = genLep->pdgId();
     genLepPt_Wm    = genLep->pt();
     genLepPhi_Wm   = genLep->phi();
     const reco::Candidate* mother = genLep->mother(0);
     genVPdgID_Wm = mother->pdgId();
     genVPt_Wm    = mother->pt();
     genVPhi_Wm   = mother->phi();
     genVy_Wm     = mother->y();
     genVMass_Wm  = mother->mass();
     TVector2 vWPt(genVPt_Wm*cos(genVPhi_Wm),genVPt_Wm*sin(genVPhi_Wm));
     TVector2 vLepPt(vLep.Px(),vLep.Py());
     TVector2 vU = -1.0*(vtype1pfMET_Wm+vLepPt);
     u1_Wm = (vWPt.Px()*vU.Px()+vWPt.Py()*vU.Py())/genVPt_Wm; // u1 = (pT . u)/|pT|
     u2_Wm = (vWPt.Px()*vU.Px()-vWPt.Py()*vU.Py())/genVPt_Wm; // u1 = (pT x u)/|pT|
   }

   // Fill tree
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

  //
  // Set up output ntuple
  //

  outFile_Wm = new TFile(outFilename_Wm,"RECREATE");
  outTree_Wm = new TTree("Events","Events");

  outTree_Wm->Branch("npv",           &npv_Wm,           "npv/I");          // number of primary vertices
  outTree_Wm->Branch("nEvents",       &nEvents_Wm,       "nEvents/I");      // total number of MC events processed
                                                                            // (will be used for event weighting in signal extraction)
  outTree_Wm->Branch("genVPt",        &genVPt_Wm,        "genVPt/F");       // GEN boson pT (signal MC)
  outTree_Wm->Branch("genVPhi",       &genVPhi_Wm,       "genVPhi/F");      // GEN boson phi (signal MC)
  outTree_Wm->Branch("genVy",         &genVy_Wm,         "genVy/F");        // GEN boson rapidity (signal MC)
  outTree_Wm->Branch("genVMass",      &genVMass_Wm,      "genVMass/F");     // GEN boson mass (signal MC)
  outTree_Wm->Branch("genLepPt",      &genLepPt_Wm,      "genLepPt/F");     // GEN lepton pT (signal MC)
  outTree_Wm->Branch("genLepPhi",     &genLepPhi_Wm,     "genLepPhi/F");    // GEN lepton phi (signal MC)
  outTree_Wm->Branch("rawpfmet",      &rawpfmet_Wm,      "rawpfmet/F");     // Raw PF MET
  outTree_Wm->Branch("rawpfmetPhi",   &rawpfmetPhi_Wm,   "rawpfmetPhi/F");  // Raw PF MET phi
  outTree_Wm->Branch("type1pfmet",    &type1pfmet_Wm,    "type1pfmet/F");   // Type-1 corrected PF MET
  outTree_Wm->Branch("type1pfmetPhi", &type1pfmetPhi_Wm, "type1pfmetPhi/F");// Type-1 corrected PF MET phi
  outTree_Wm->Branch("genmet",        &genmet_Wm,        "genmet/F");       // Generator level MET
  outTree_Wm->Branch("genmetPhi",     &genmetPhi_Wm,     "genmetPhi/F");    // Generator level MET phi
  outTree_Wm->Branch("mt",            &mt_Wm,            "mt/F");           // transverse mass
  outTree_Wm->Branch("u1",            &u1_Wm,            "u1/F");           // parallel component of recoil
  outTree_Wm->Branch("u2",            &u2_Wm,            "u2/F");           // perpendicular component of recoil 
  outTree_Wm->Branch("q",             &q_Wm,             "q/I");            // lepton charge
  outTree_Wm->Branch("lep", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep_Wm);   // lepton 4-vector
  ///// muon specific /////
  outTree_Wm->Branch("pfChIso",       &pfChIso_Wm,       "pfChIso/F");      // PF charged hadron isolation of muon
  outTree_Wm->Branch("pfGamIso",      &pfGamIso_Wm,      "pfGamIso/F");     // PF photon isolation of muon
  outTree_Wm->Branch("pfNeuIso",      &pfNeuIso_Wm,      "pfNeuIso/F");     // PF neutral hadron isolation of muon
  outTree_Wm->Branch("isLooseMuon",   &isLooseMuon_Wm,  "isLooseMuon/I");   // loose muon ID
  outTree_Wm->Branch("isSoftMuon",    &isSoftMuon_Wm,   "isSoftMuon/I");    // loose muon ID
  outTree_Wm->Branch("isTightMuon",   &isTightMuon_Wm,  "isTightMuon/I");   // tight muon ID

}

// ------------ method called once each job just after ending the event loop  ------------
void 
selectWm::endJob() 
{
   // The only information in the last entry of the output tree
   // is the total number of MC events processed
   nEvents_Wm = eventCounter_Wm;
   // Set everything else to zero
   npv_Wm=0;
   genVPdgID_Wm=0; genVPt_Wm=0; genVPhi_Wm=0; genVy_Wm=0; genVMass_Wm=0;
   genLepPdgID_Wm=0; genLepPt_Wm=0; genLepPhi_Wm=0;
   rawpfmet_Wm=0; rawpfmetPhi_Wm=0;
   type1pfmet_Wm=0; type1pfmetPhi_Wm=0;
   genmet_Wm=0; genmetPhi_Wm=0;
   mt_Wm=0; u1_Wm=0; u2_Wm=0;
   q_Wm=0;
   lep_Wm=0;
   ///// muon specific /////
   pfChIso_Wm=0; pfGamIso_Wm=0; pfNeuIso_Wm=0;
   isLooseMuon_Wm=0; isSoftMuon_Wm=0; isTightMuon_Wm=0;

   outTree_Wm->Fill();

   // Save tree in output file
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
  std::cout << nsel_Wm << " +/- " << sqrt(nsel_Wm) << " per 1/fb" << std::endl;
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
