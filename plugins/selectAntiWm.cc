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
// https://github.com/jaylawhorn/mitewk/blob/master/Selection/selectAntiWm.C

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

class selectAntiWm : public edm::EDAnalyzer {
   public:
      explicit selectAntiWm(const edm::ParameterSet&);
      ~selectAntiWm();

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

// Convert int to binary value
bool Convert_AWm(unsigned int val,bool print=kFALSE)
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

const Double_t PT_CUT    = 20;
const Double_t ETA_CUT   = 2.4;
const Double_t MUON_MASS = 0.105658369;

Double_t nsel_AWm=0, nselvar_AWm=0;

TString outFilename_AWm = TString("Wmunu_select.root");
TFile *outFile_AWm = new TFile();
TTree *outTree_AWm = new TTree();

//
// Declare output ntuple variables
//
UInt_t  npv_AWm;
Float_t genVPdgID_AWm, genVPt_AWm, genVPhi_AWm, genVy_AWm, genVMass_AWm;
Float_t genLepPdgID_AWm, genLepPt_AWm, genLepPhi_AWm;
Float_t scale1fb_AWm;
Float_t rawpfmet_AWm, rawpfmetPhi_AWm;
Float_t type1pfmet_AWm, type1pfmetPhi_AWm;
Float_t genmet_AWm, genmetPhi_AWm;
Float_t mt_AWm, u1_AWm, u2_AWm;
Int_t q_AWm;
LorentzVector *lep_AWm=0;
///// muon specific /////
Float_t pfChIso_AWm, pfGamIso_AWm, pfNeuIso_AWm;
Float_t isLooseMuon_AWm, isSoftMuon_AWm, isTightMuon_AWm;
Bool_t passMuonLooseID_AWm, passMuonSoftID_AWm, passMuonTightID_AWm;

// Compute MC event weightWm_sel per 1/fb
Double_t weight_AWm = 1;
//const Double_t xsec = 1;
//if(xsec>0) weightWm = 1000.*xsec/(Double_t)eventTree->GetEntries();

//
// constructors and destructor
//
selectAntiWm::selectAntiWm(const edm::ParameterSet& iConfig):
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
   muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
   metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
   pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands")))
{
   //now do what ever initialization is needed

}


selectAntiWm::~selectAntiWm()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
selectAntiWm::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   // good vertex requirement
   if (vertices->empty()) return; // skip the event if no PV found
   const reco::Vertex &PV = vertices->front();
   npv_AWm = vertices->size();

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

     isLooseMuon_AWm = mu.isLooseMuon(); passMuonLooseID_AWm = Convert_AWm(isLooseMuon_AWm);
     isTightMuon_AWm = mu.isTightMuon(PV); passMuonTightID_AWm = Convert_AWm(isTightMuon_AWm);

     if(fabs(mu.eta()) > 2.4) continue;    // lepton |eta| cut
     if(mu.pt()        < 10 ) continue;    // lepton pT cut
     if(passMuonLooseID_AWm)   nLooseLep++; // loose lepton selection
     if(nLooseLep>1) { // extra lepton veto
       passSel=kFALSE;
       break;
     }

     Double_t iso = mu.chargedHadronIso() + TMath::Max(mu.neutralHadronIso() + mu.photonIso() - 0.5*mu.puChargedHadronIso(),Double_t(0));

     if(fabs(mu.eta()) > ETA_CUT) continue;                          // lepton |eta| cut
     if(mu.pt()        < PT_CUT ) continue;                          // lepton pT cut
     if(!(passMuonTightID_AWm && iso >= 0.3*mu.pt()))      continue; // lepton selection

     passSel = kTRUE;
     goodMuonIdx = jMuon;
   }
   if(passSel==kFALSE) return;

   const pat::Muon &goodMuon = (*muons)[goodMuonIdx];
   // Lepton information
   q_AWm           = goodMuon.charge();
   pfChIso_AWm     = goodMuon.chargedHadronIso();
   pfGamIso_AWm    = goodMuon.photonIso();
   pfNeuIso_AWm    = goodMuon.neutralHadronIso();
   isLooseMuon_AWm = goodMuon.isLooseMuon(); passMuonLooseID_AWm = Convert_AWm(isLooseMuon_AWm);
   isSoftMuon_AWm  = goodMuon.isSoftMuon(PV); passMuonSoftID_AWm = Convert_AWm(isSoftMuon_AWm);
   isTightMuon_AWm = goodMuon.isTightMuon(PV); passMuonTightID_AWm = Convert_AWm(isTightMuon_AWm);

   // Event counting and scale factors
   nsel_AWm    += weight_AWm;
   nselvar_AWm += weight_AWm*weight_AWm;
   scale1fb_AWm = weight_AWm;

   // Lepton and supercluster 4-vectors
   LorentzVector vLep(goodMuon.pt(),goodMuon.eta(),goodMuon.phi(),MUON_MASS);

   lep_AWm = &vLep;

   // Type-1 corrected PF MET (default)
   edm::Handle<pat::METCollection> mets;
   iEvent.getByToken(metToken_, mets);

   const pat::MET &met = mets->front();
   TVector2 vtype1pfMET_AWm = TVector2(met.px(),met.py());
   type1pfmet_AWm    = vtype1pfMET_AWm.Mod();
   type1pfmetPhi_AWm = vtype1pfMET_AWm.Phi();

   // Generator level MET
   TVector2 vgenMET_AWm = TVector2(met.genMET()->px(),met.genMET()->py());
   genmet_AWm    = vgenMET_AWm.Mod();
   genmetPhi_AWm = vgenMET_AWm.Phi();

   // Raw PF MET
   edm::Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);

   Float_t rawpfMETpx=0, rawpfMETpy=0;
   for(unsigned int jcand=0; jcand < pfs->size(); jcand++) {
     const pat::PackedCandidate &pf = (*pfs)[jcand];
     rawpfMETpx -= pf.px();
     rawpfMETpy -= pf.py();
   }

   TVector2 vrawpfMET_AWm = TVector2(rawpfMETpx,rawpfMETpy);
   rawpfmet_AWm    = vrawpfMET_AWm.Mod();
   rawpfmetPhi_AWm = vrawpfMET_AWm.Phi();

   // Lepton transverse mass
   mt_AWm       = sqrt(2.0*vLep.Pt()*vtype1pfMET_AWm.Mod()*(1.0-cos(vLep.Phi()-vtype1pfMET_AWm.Phi())));

   // Generator level lepton information and hadronic recoil
   const reco::GenParticle* genLep = goodMuon.genLepton();
   if(genLep!=NULL) {
     genLepPdgID_AWm = genLep->pdgId();
     genLepPt_AWm    = genLep->pt();
     genLepPhi_AWm   = genLep->phi();
     const reco::Candidate* mother = genLep->mother(0);
     genVPdgID_AWm = mother->pdgId();
     genVPt_AWm    = mother->pt();
     genVPhi_AWm   = mother->phi();
     genVy_AWm     = mother->y();
     genVMass_AWm  = mother->mass();
     TVector2 vWPt(genVPt_AWm*cos(genVPhi_AWm),genVPt_AWm*sin(genVPhi_AWm));
     TVector2 vLepPt(vLep.Px(),vLep.Py());
     TVector2 vU = -1.0*(vtype1pfMET_AWm+vLepPt);
     u1_AWm = (vWPt.Px()*vU.Px()+vWPt.Py()*vU.Py())/genVPt_AWm; // u1 = (pT . u)/|pT|
     u2_AWm = (vWPt.Px()*vU.Px()-vWPt.Py()*vU.Py())/genVPt_AWm; // u1 = (pT x u)/|pT|
   }

   // Fill tree
   outTree_AWm->Fill();

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
selectAntiWm::beginJob()
{

  //
  // Set up output ntuple
  //

  outFile_AWm = new TFile(outFilename_AWm,"RECREATE");
  outTree_AWm = new TTree("Events","Events");

  outTree_AWm->Branch("npv",           &npv_AWm,           "npv/I");          // number of primary vertices
  outTree_AWm->Branch("genVPt",        &genVPt_AWm,        "genVPt/F");       // GEN boson pT (signal MC)
  outTree_AWm->Branch("genVPhi",       &genVPhi_AWm,       "genVPhi/F");      // GEN boson phi (signal MC)
  outTree_AWm->Branch("genVy",         &genVy_AWm,         "genVy/F");        // GEN boson rapidity (signal MC)
  outTree_AWm->Branch("genVMass",      &genVMass_AWm,      "genVMass/F");     // GEN boson mass (signal MC)
  outTree_AWm->Branch("genLepPt",      &genLepPt_AWm,      "genLepPt/F");     // GEN lepton pT (signal MC)
  outTree_AWm->Branch("genLepPhi",     &genLepPhi_AWm,     "genLepPhi/F");    // GEN lepton phi (signal MC)
  outTree_AWm->Branch("scale1fb",      &scale1fb_AWm,      "scale1fb/F");     // event weight per 1/fb (MC)
  outTree_AWm->Branch("rawpfmet",      &rawpfmet_AWm,      "rawpfmet/F");     // Raw PF MET
  outTree_AWm->Branch("rawpfmetPhi",   &rawpfmetPhi_AWm,   "rawpfmetPhi/F");  // Raw PF MET phi
  outTree_AWm->Branch("type1pfmet",    &type1pfmet_AWm,    "type1pfmet/F");   // Type-1 corrected PF MET
  outTree_AWm->Branch("type1pfmetPhi", &type1pfmetPhi_AWm, "type1pfmetPhi/F");// Type-1 corrected PF MET phi
  outTree_AWm->Branch("genmet",        &genmet_AWm,        "genmet/F");       // Generator level MET
  outTree_AWm->Branch("genmetPhi",     &genmetPhi_AWm,     "genmetPhi/F");    // Generator level MET phi
  outTree_AWm->Branch("mt",            &mt_AWm,            "mt/F");           // transverse mass
  outTree_AWm->Branch("u1",            &u1_AWm,            "u1/F");           // parallel component of recoil
  outTree_AWm->Branch("u2",            &u2_AWm,            "u2/F");           // perpendicular component of recoil 
  outTree_AWm->Branch("q",             &q_AWm,             "q/I");            // lepton charge
  outTree_AWm->Branch("lep", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep_AWm);   // lepton 4-vector
  ///// muon specific /////
  outTree_AWm->Branch("pfChIso",       &pfChIso_AWm,       "pfChIso/F");      // PF charged hadron isolation of muon
  outTree_AWm->Branch("pfGamIso",      &pfGamIso_AWm,      "pfGamIso/F");     // PF photon isolation of muon
  outTree_AWm->Branch("pfNeuIso",      &pfNeuIso_AWm,      "pfNeuIso/F");     // PF neutral hadron isolation of muon
  outTree_AWm->Branch("isLooseMuon",   &isLooseMuon_AWm,  "isLooseMuon/F");  // loose muon ID
  outTree_AWm->Branch("isSoftMuon",    &isSoftMuon_AWm,    "isSoftMuon/F");   // loose muon ID
  outTree_AWm->Branch("isTightMuon",   &isTightMuon_AWm,  "isTightMuon/F");  // tight muon ID

}

// ------------ method called once each job just after ending the event loop  ------------
void 
selectAntiWm::endJob() 
{

   // Save tree in output file
   outFile_AWm->Write();
   outFile_AWm->Close();

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================

  std::cout << "*" << std::endl;
  std::cout << "* SUMMARY" << std::endl;
  std::cout << "*--------------------------------------------------" << std::endl;
  std::cout << "W -> mu nu anti-selection" << std::endl;
  std::cout << " pT > " << PT_CUT << std::endl;
  std::cout << " |eta| < " << ETA_CUT << std::endl;
  std::cout << nsel_AWm << " +/- " << sqrt(nselvar_AWm) << " per 1/fb" << std::endl;
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "  <> Output saved in " << outFilename_AWm << "/" << std::endl;
  std::cout << std::endl;
}

// ------------ method called when starting to processes a run  ------------
/*
void 
selectAntiWm::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
selectAntiWm::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
selectAntiWm::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
selectAntiWm::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
selectAntiWm::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(selectAntiWm);
