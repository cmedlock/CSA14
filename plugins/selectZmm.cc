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
const Double_t MUON_MASS = 0.105658369;

Double_t nsel_Zmm=0, nselvar_Zmm=0;

TString outFilename_Zmm = TString("selectZmm.root");
TFile *outFile_Zmm = new TFile();
TTree *outTree_Zmm = new TTree();

//
// Declare output ntuple variables
//
UInt_t  matchGen_Zmm, npv_Zmm;
Float_t genVPdgID_Zmm, genVPt_Zmm, genVPhi_Zmm, genVy_Zmm, genVMass_Zmm;
Float_t scale1fb_Zmm;
Float_t rawpfmet_Zmm, rawpfmetPhi_Zmm;
Float_t type1pfmet_Zmm, type1pfmetPhi_Zmm;
Float_t genmet_Zmm, genmetPhi_Zmm;
Float_t u1_Zmm, u2_Zmm;
Int_t q1_Zmm, q2_Zmm;
LorentzVector *dilep_Zmm=0, *lep1_Zmm=0, *lep2_Zmm=0;
Bool_t passMuonLooseID1_Zmm, passMuonSoftID1_Zmm, passMuonTightID1_Zmm;
Bool_t passMuonLooseID2_Zmm, passMuonSoftID2_Zmm, passMuonTightID2_Zmm;

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

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   // Good vertex requirement
   if (vertices->empty()) return; // Skip the event if no PV found
   const reco::Vertex &PV = vertices->front();
   npv_Zmm       = vertices->size();

   Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);
   if(muons->size()<2) return; // Skip the event if there is no possibility of both tag and probe leptons

   // Look for tag lepton
   unsigned int tagIdx=muons->size();
   Float_t tagPt=0;
   Float_t isTightMuon=0;
   Bool_t passMuonTightID=kFALSE;
   for (unsigned int jMuon=0;jMuon < muons->size();jMuon++) {
     const pat::Muon &mu = (*muons)[jMuon];

     isTightMuon = mu.isTightMuon(PV); passMuonTightID = Convert_Zmm(isTightMuon);

     if(mu.pt()>tagPt && passMuonTightID) {
       tagPt  = mu.pt();
       tagIdx = jMuon;
     } else {
       tagPt  = tagPt;
       tagIdx = tagIdx;
     }
   }
   if(tagIdx==muons->size()) return; // Skip event if there is no tag lepton

   const pat::Muon &tag = (*muons)[tagIdx];
   // Tag lepton information
   q1_Zmm = tag.charge();
   Float_t isLooseMuon1 = tag.isLooseMuon();   passMuonLooseID1_Zmm = Convert_Zmm(isLooseMuon1);
   Float_t isSoftMuon1  = tag.isSoftMuon(PV);  passMuonSoftID1_Zmm  = Convert_Zmm(isSoftMuon1);
   Float_t isTightMuon1 = tag.isTightMuon(PV); passMuonTightID1_Zmm = Convert_Zmm(isTightMuon1);

   // Look for probe lepton
   unsigned int probeIdx=0;
   Float_t probePt=0;
   for (unsigned int kMuon=0;kMuon < muons->size();kMuon++) {
     if(kMuon==tagIdx) continue;
     const pat::Muon &mu = (*muons)[kMuon];
     if(mu.pt()>probePt) {
       probePt  = mu.pt();
       probeIdx = kMuon;
     }
   }

   const pat::Muon &probe = (*muons)[probeIdx];
   // Probe lepton information
   q2_Zmm = probe.charge();
   Float_t isLooseMuon2 = probe.isLooseMuon();   passMuonLooseID2_Zmm = Convert_Zmm(isLooseMuon2);
   Float_t isSoftMuon2  = probe.isSoftMuon(PV);  passMuonSoftID2_Zmm  = Convert_Zmm(isSoftMuon2);
   Float_t isTightMuon2 = probe.isTightMuon(PV); passMuonTightID2_Zmm = Convert_Zmm(isTightMuon2);

   // Lepton and supercluster 4-vectors
   LorentzVector vLep1(tag.pt(),tag.eta(),tag.phi(),MUON_MASS);
   LorentzVector vLep2(probe.pt(),probe.eta(),probe.phi(),MUON_MASS);
   LorentzVector vDilep = vLep1 + vLep2;

   lep1_Zmm     = &vLep1;
   lep2_Zmm     = &vLep2;
   dilep_Zmm    = &vDilep;

   // Mass window requirement
   if((vDilep.M()<MASS_LOW) || (vDilep.M()>MASS_HIGH)) return;

   // Perform matching of dileptons to generator level leptons
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

   // Event counting and scale factors
   nsel_Zmm    += weight_Zmm;
   nselvar_Zmm += weight_Zmm*weight_Zmm;
   scale1fb_Zmm = weight_Zmm;

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
   const reco::GenParticle* genLep = tag.genLepton();
   if(genLep!=NULL) {
     const reco::Candidate* mother = genLep->mother(0);
     genVPdgID_Zmm = mother->pdgId();
     genVPt_Zmm    = mother->pt();
     genVPhi_Zmm   = mother->phi();
     genVy_Zmm     = mother->y();
     genVMass_Zmm  = mother->mass();
     TVector2 vVPt(genVPt_Zmm*cos(genVPhi_Zmm),genVPt_Zmm*sin(genVPhi_Zmm));
     TVector2 vLepPt(vLep1.Px(),vLep1.Py());
     TVector2 vU = -1.0*(vtype1pfMET_Zmm+vLepPt);
     u1_Zmm = (vVPt.Px()*vU.Px()+vVPt.Py()*vU.Py())/genVPt_Zmm; // u1 = (pT . u)/|pT|
     u2_Zmm = (vVPt.Px()*vU.Px()-vVPt.Py()*vU.Py())/genVPt_Zmm; // u2 = (pT x u)/|pT|
   }

   // Fill tree
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

  //
  // Set up output ntuple
  //

  outFile_Zmm = new TFile(outFilename_Zmm,"RECREATE");
  outTree_Zmm = new TTree("Events","Events");

  outTree_Zmm->Branch("matchGen",      &matchGen_Zmm,      "matchGen/I");     // event has both leptons matched to MC Z->ll
  outTree_Zmm->Branch("npv",           &npv_Zmm,           "npv/I");          // number of vertices
  outTree_Zmm->Branch("genVPt",        &genVPt_Zmm,        "genVPt/F");       // GEN boson pT (signal MC)
  outTree_Zmm->Branch("genVPhi",       &genVPhi_Zmm,       "genVPhi/F");      // GEN boson phi (signal MC)
  outTree_Zmm->Branch("genVy",         &genVy_Zmm,         "genVy/F");        // GEN boson rapidity (signal MC)
  outTree_Zmm->Branch("genVMass",      &genVMass_Zmm,      "genVMass/F");     // GEN boson mass (signal MC)
  outTree_Zmm->Branch("scale1fb",      &scale1fb_Zmm,      "scale1fb/F");     // event weight_Zmm per 1/fb (MC)
  outTree_Zmm->Branch("rawpfmet",      &rawpfmet_Zmm,      "rawpfmet/F");         // Raw PF MET
  outTree_Zmm->Branch("rawpfmetPhi",   &rawpfmetPhi_Zmm,   "rawpfmetPhi/F");      // Raw PF MET phi
  outTree_Zmm->Branch("type1pfmet",    &type1pfmet_Zmm,    "type1pfmet/F");       // Type-1 corrected PF MET
  outTree_Zmm->Branch("type1pfmetPhi", &type1pfmetPhi_Zmm, "type1pfmetPhi/F");    // Type-1 corrected PF MET phi
  outTree_Zmm->Branch("genmet",        &genmet_Zmm,        "genmet/F");           // Generator level MET
  outTree_Zmm->Branch("genmetPhi",     &genmetPhi_Zmm,     "genmetPhi/F");        // Generator level MET phi
  outTree_Zmm->Branch("u1",            &u1_Zmm,            "u1/F");           // parallel component of recoil
  outTree_Zmm->Branch("u2",            &u2_Zmm,            "u2/F");           // perpendicular component of recoil
  outTree_Zmm->Branch("q1",            &q1_Zmm,            "q1/I");           // lepton charge (tag lepton )
  outTree_Zmm->Branch("q2",            &q2_Zmm,            "q2/I");           // lepton charge (probe lepton)
  outTree_Zmm->Branch("dilep","ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &dilep_Zmm); // dilepton 4-vector
  outTree_Zmm->Branch("lep1", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep1_Zmm);  // lepton 4-vector (tag lepton)
  outTree_Zmm->Branch("lep2", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep2_Zmm);  // lepton 4-vector (probe lepton)
  outTree_Zmm->Branch("passMuonLooseID1",  &passMuonLooseID1_Zmm, "passMuonLooseID1/O");  // tag lepton loose muon ID
  outTree_Zmm->Branch("passMuonSoftID1",   &passMuonSoftID1_Zmm,  "passMuonSoftID1/O");   // tag lepton loose muon ID
  outTree_Zmm->Branch("passMuonTightID1",  &passMuonTightID1_Zmm, "passMuonTightID1/O");  // tag lepton tight muon ID
  outTree_Zmm->Branch("passMuonLooseID2",  &passMuonLooseID2_Zmm, "passMuonLooseID2/O");  // probe lepton loose muon ID
  outTree_Zmm->Branch("passMuonSoftID2",   &passMuonSoftID2_Zmm,  "passMuonSoftID2/O");   // tag lepton loose muon ID
  outTree_Zmm->Branch("passMuonTightID2",  &passMuonTightID2_Zmm, "passMuonTightID2/O");  // probe lepton tight muon ID
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
  std::cout << nsel_Zmm << " +/- " << sqrt(nselvar_Zmm) << " per 1/fb" << std::endl;
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
