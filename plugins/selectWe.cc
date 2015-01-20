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
//         Created:  Tue, 22 Jul 2014 11:26:17 Gmt_We
// This code is just an adaptation of Kevin Sung's original code used for the 8 TeV analysis:
// https://github.com/jaylawhorn/mitewk/blob/master/Selection/selectWe.C

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
      edm::EDGetTokenT<double> rhoToken_;
//      edm::InputTag rhoFixedGrid;
};

//
// constants, enums and typedefs
//

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

// Convert int to binary value
bool Convert_We(unsigned int val,bool print=kFALSE)
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
const Double_t ETA_CUT   = 2.5;
const Double_t ELE_MASS  = 0.000511;

const Double_t ECAL_GAP_LOW  = 1.4442;
const Double_t ECAL_GAP_HIGH = 1.566;

Double_t nsel_We=0, nselvar_We=0;

TString outFilename = TString("Wenu_select.root");
TFile *outFile = new TFile();
TTree *outTree_We = new TTree();

//
// Declare output ntuple variables
//
UInt_t  npv_We;
Float_t genVPdgID_We, genVPt_We, genVPhi_We, genVy_We, genVMass_We;
Float_t genLepPdgID_We, genLepPt_We, genLepPhi_We;
Float_t scale1fb_We;
Float_t rawpfmet_We, rawpfmetPhi_We;
Float_t type1pfmet_We, type1pfmetPhi_We;
Float_t genmet_We, genmetPhi_We;
Float_t mt_We, u1_We, u2_We;
Int_t q_We;
LorentzVector *lep_We=0;
///// electron specific /////
Float_t pfChIso_We, pfGamIso_We, pfNeuIso_We;
Float_t isVetoEle_We, isLooseEle_We, isMediumEle_We, isTightEle_We;
LorentzVector *sc_We=0;

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
   pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
   rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhos")))
{
   //now do what ever initialization is needed
//   rhoFixedGrid = iConfig.getParameter<edm::InputTag>("fixedGridRhoAll");
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

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   // good vertex requirement
   if (vertices->empty()) return; // skip the event if no PV found
   npv_We = vertices->size();
//   const reco::Vertex &PV = vertices->front();

   Handle<double> rhoHandle;
//   iEvent.getByLabel(rhoFixedGrid, rhoHandle);
   iEvent.getByToken(rhoToken_, rhoHandle);
   Double_t rho = *rhoHandle.product();

   Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);
   if(electrons->size()==0) return; // skip the event if there are no electrons

   //
   // SELECTION PROCEDURE:
   //  (1) Look for 1 good electron matched to trigger
   //  (2) Reject event if another electron is present passing looser cuts
   //
   Int_t  nLooseLep=0;
   Bool_t passSel=kFALSE;
   Int_t  goodEleIdx=0;
   for (unsigned int jElectron=0;jElectron < electrons->size();jElectron++) {
     const pat::Electron &ele = (*electrons)[jElectron];

     // check ECAL gap
     if(fabs(ele.superCluster()->eta())>=ECAL_GAP_LOW && fabs(ele.superCluster()->eta())<=ECAL_GAP_HIGH) continue;

     Double_t ea = 0;
     if (fabs(ele.eta()) < 1.0) ea = 0.100;
     else if(fabs(ele.eta()) < 1.479) ea = 0.120;
     else if(fabs(ele.eta()) < 2.0) ea = 0.085;
     else if(fabs(ele.eta()) < 2.2) ea = 0.110;
     else if(fabs(ele.eta()) < 2.3) ea = 0.120;
     else if(fabs(ele.eta()) < 2.4) ea = 0.120;
     else ea = 0.130;

     Double_t iso = ele.chargedHadronIso() + TMath::Max(ele.neutralHadronIso() + ele.photonIso() - rho*ea,0.);

     isLooseEle_We = ele.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-loose");
     isTightEle_We = ele.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-tight");

     if(fabs(ele.superCluster()->eta())> 2.5)  continue;    // lepton |eta| cut
     if(ele.superCluster()->energy()   < 20 )  continue;    // lepton pT cut
     if(isLooseEle_We && iso <= 0.15*ele.pt()) nLooseLep++; // loose lepton selection
     if(nLooseLep>1) { // extra lepton veto
       passSel=kFALSE;
       break;
     }

     if(fabs(ele.superCluster()->eta()) > ETA_CUT) continue; // lepton |eta| cut
     if(ele.superCluster()->energy()    < PT_CUT ) continue; // lepton pT cut
     if(!(isTightEle_We && iso <= 0.15*ele.pt()))  continue; // lepton selection

     passSel = kTRUE;
     goodEleIdx = jElectron;
   }
   if(passSel==kFALSE) return;

   const pat::Electron &goodEle = (*electrons)[goodEleIdx];
   // Lepton information
   q_We           = goodEle.charge();
   pfChIso_We     = goodEle.chargedHadronIso();
   pfGamIso_We    = goodEle.photonIso();
   pfNeuIso_We    = goodEle.neutralHadronIso();
   isVetoEle_We   = goodEle.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-veto");
   isLooseEle_We  = goodEle.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-loose");
   isMediumEle_We = goodEle.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-medium");
   isTightEle_We  = goodEle.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-tight");

   // Event counting and scale factors
   nsel_We    += weight;
   nselvar_We += weight*weight;
   scale1fb_We = weight;

   // Lepton and supercluster 4-vectors
   LorentzVector vLep(goodEle.pt(),goodEle.eta(),goodEle.phi(),ELE_MASS);
   LorentzVector vSC(goodEle.superCluster()->energy(),goodEle.superCluster()->eta(),goodEle.superCluster()->phi(),ELE_MASS);

   lep_We = &vLep;
   sc_We = &vSC;

   // Type-1 corrected PF MET (default)
   edm::Handle<pat::METCollection> mets;
   iEvent.getByToken(metToken_, mets);

   const pat::MET &met = mets->front();
   TVector2 vtype1pfMET_We = TVector2(met.px(),met.py());
   type1pfmet_We    = vtype1pfMET_We.Mod();
   type1pfmetPhi_We = vtype1pfMET_We.Phi();

   // Generator level MET
   TVector2 vgenMET_We = TVector2(met.genMET()->px(),met.genMET()->py());
   genmet_We    = vgenMET_We.Mod();
   genmetPhi_We = vgenMET_We.Phi();

   // Raw PF MET
   edm::Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);
   Float_t rawpfMETpx=0, rawpfMETpy=0;

   for(unsigned int jcand=0; jcand < pfs->size(); jcand++) {
     const pat::PackedCandidate &pf = (*pfs)[jcand];
     rawpfMETpx -= pf.px();
     rawpfMETpy -= pf.py();
   }

   TVector2 vrawpfMET_We = TVector2(rawpfMETpx,rawpfMETpy);
   rawpfmet_We    = vrawpfMET_We.Mod();
   rawpfmetPhi_We = vrawpfMET_We.Phi();

   // Lepton transverse mass
   mt_We = sqrt(2.0*vLep.Pt()*vtype1pfMET_We.Mod()*(1.0-cos(vLep.Phi()-vtype1pfMET_We.Phi())));

   // Generator level lepton information and hadronic recoil
   const reco::GenParticle* genLep = goodEle.genLepton();
   if(genLep!=NULL) {
     genLepPdgID_We = genLep->pdgId();
     genLepPt_We    = genLep->pt();
     genLepPhi_We   = genLep->phi();
     const reco::Candidate* mother = genLep->mother(0);
     genVPdgID_We = mother->pdgId();
     genVPt_We    = mother->pt();
     genVPhi_We   = mother->phi();
     genVy_We     = mother->y();
     genVMass_We  = mother->mass();
     TVector2 vWPt(genVPt_We*cos(genVPhi_We),genVPt_We*sin(genVPhi_We));
     TVector2 vLepPt(vLep.Px(),vLep.Py());
     TVector2 vU = -1.0*(vtype1pfMET_We+vLepPt);
     u1_We = (vWPt.Px()*vU.Px()+vWPt.Py()*vU.Py())/genVPt_We; // u1_We = (pT . u)/|pT|
     u2_We = (vWPt.Px()*vU.Px()-vWPt.Py()*vU.Py())/genVPt_We; // u1_We = (pT x u)/|pT|
   }

   // Fill tree
   outTree_We->Fill();

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

  //
  // Set up output ntuple
  //

  outFile = new TFile(outFilename,"RECREATE");
  outTree_We = new TTree("Events","Events");

  outTree_We->Branch("npv",           &npv_We,           "npv/I");          // number of primary vertices
  outTree_We->Branch("genVPt",        &genVPt_We,        "genVPt/F");       // GEN boson pT (signal MC)
  outTree_We->Branch("genVPhi",       &genVPhi_We,       "genVPhi/F");      // GEN boson phi (signal MC)
  outTree_We->Branch("genVy",         &genVy_We,         "genVy/F");        // GEN boson rapidity (signal MC)
  outTree_We->Branch("genVMass",      &genVMass_We,      "genVMass/F");     // GEN boson mass (signal MC)
  outTree_We->Branch("genLepPt",      &genLepPt_We,      "genLepPt/F");     // GEN lepton pT (signal MC)
  outTree_We->Branch("genLepPhi",     &genLepPhi_We,     "genLepPhi/F");    // GEN lepton phi (signal MC)
  outTree_We->Branch("scale1fb",      &scale1fb_We,      "scale1fb/F");     // event weight per 1/fb (MC)
  outTree_We->Branch("rawpfmet",      &rawpfmet_We,      "rawpfmet/F");     // Raw PF MET
  outTree_We->Branch("rawpfmetPhi",   &rawpfmetPhi_We,   "rawpfmetPhi/F");  // Raw PF MET phi
  outTree_We->Branch("type1pfmet",    &type1pfmet_We,    "type1pfmet/F");   // Type-1 corrected PF MET
  outTree_We->Branch("type1pfmetPhi", &type1pfmetPhi_We, "type1pfmetPhi/F");// Type-1 corrected PF MET phi
  outTree_We->Branch("genmet",        &genmet_We,        "genmet/F");       // Generator level MET
  outTree_We->Branch("genmetPhi",     &genmetPhi_We,     "genmetPhi/F");    // Generator level MET phi
  outTree_We->Branch("mt",            &mt_We,            "mt/F");           // transverse mass
  outTree_We->Branch("u1",            &u1_We,            "u1/F");           // parallel component of recoil
  outTree_We->Branch("u2",            &u2_We,            "u2/F");           // perpendicular component of recoil 
  outTree_We->Branch("q",             &q_We,             "q/I");            // lepton charge
  outTree_We->Branch("lep", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep_We);   // lepton 4-vector
  outTree_We->Branch("sc",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sc_We);    // supercluster 4-vector
  ///// electron specific /////
  outTree_We->Branch("pfChIso",       &pfChIso_We,       "pfChIso/F");      // PF charged hadron isolation of electron
  outTree_We->Branch("pfGamIso",      &pfGamIso_We,      "pfGamIso/F");     // PF photon isolation of electron
  outTree_We->Branch("pfNeuIso",      &pfNeuIso_We,      "pfNeuIso/F");     // PF neutral hadron isolation of electron
  outTree_We->Branch("isVetoEle",     &isVetoEle_We,     "isVetoEle/F");    // tag lepton veto electron ID
  outTree_We->Branch("isLooseEle",    &isLooseEle_We,    "isLooseEle/F");   // tag lepton loose electron ID
  outTree_We->Branch("isMediumEle",   &isMediumEle_We,   "isMediumEle/F");  // tag lepton medium electron ID
  outTree_We->Branch("isTightEle",    &isTightEle_We,    "isTightEle/F");   // tag lepton tight electron ID
}

// ------------ method called once each job just after ending the event loop  ------------
void 
selectWe::endJob() 
{
   // Save tree to output file
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
  std::cout << nsel_We << " +/- " << sqrt(nselvar_We) << " per 1/fb" << std::endl;
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "  <> Output saved in " << outFilename << "/" << std::endl;
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
