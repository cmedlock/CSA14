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

Double_t nsel=0, nselvar=0;

TString outFilename = TString("Wenu_select.root");
TFile *outFile = new TFile();
TTree *outTree = new TTree();

//
// Declare output ntuple variables
//
UInt_t   npv;
Float_t genVPdgID, genVPt, genVPhi, genVy, genVMass;
Float_t genLepPdgID, genLepPt, genLepPhi;
Float_t scale1fb;
Float_t rawpfMETpx, rawpfMETpy;
TVector2 vtype1pfMET, vrawpfMET, vgenMET;
Float_t mt, u1, u2;
Int_t q;
Int_t dummynEvents=0, nEvents=0;
LorentzVector *lep=0;
///// electron specific /////
Float_t pfChIso, pfGamIso, pfNeuIso;
Float_t isLooseEle, isTightEle;
Bool_t passEleLooseID, passEleTightID;
LorentzVector *sc=0;

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

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   // good vertex requirement
   if (vertices->empty()) return; // skip the event if no PV found
//   const reco::Vertex &PV = vertices->front();

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

     isLooseEle = ele.electronID("eidLoose"); passEleLooseID = Convert_We(isLooseEle);
     isTightEle = ele.electronID("eidTight"); passEleTightID = Convert_We(isTightEle);

     if(fabs(ele.superCluster()->eta())> 2.5) continue; // lepton |eta| cut
     if(ele.superCluster()->energy()   < 20 ) continue; // lepton pT cut
     if(passEleLooseID)                       nLooseLep++; // loose lepton selection
     if(nLooseLep>1) { // extra lepton veto
       passSel=kFALSE;
       break;
     }

     if(fabs(ele.superCluster()->eta()) > ETA_CUT) continue; // lepton |eta| cut
     if(ele.superCluster()->energy()    < PT_CUT ) continue; // lepton pT cut
     if(!passEleTightID)                           continue; // lepton selection

     passSel = kTRUE;
     goodEleIdx = jElectron;
   }

   const pat::Electron &goodEle = (*electrons)[goodEleIdx];

   if(passSel==kFALSE) return;

   nsel    += weight;
   nselvar += weight*weight;

   LorentzVector vLep(goodEle.pt(),goodEle.eta(),goodEle.phi(),ELE_MASS);
   LorentzVector vSC(goodEle.superCluster()->energy(),goodEle.superCluster()->eta(),goodEle.superCluster()->phi(),ELE_MASS);

   // Save MET information

   edm::Handle<pat::METCollection> mets;
   iEvent.getByToken(metToken_, mets);
   const pat::MET &met = mets->front();

   edm::Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);
   rawpfMETpx=0;
   rawpfMETpy=0;
   // loop on PF candidates to calculate sum of Et's
   for(unsigned int jcand=0; jcand < pfs->size(); jcand++) {
     const pat::PackedCandidate &pf = (*pfs)[jcand];
     rawpfMETpx -= pf.px();
     rawpfMETpy -= pf.py();
   }

   vtype1pfMET.Set(met.px(),met.py());
   vrawpfMET.Set(rawpfMETpx,rawpfMETpy);
   vgenMET.Set(met.genMET()->px(),met.genMET()->py());

   //
   // Fill tree
   //
   npv = vertices->size();
   genVPdgID   = 0;
   genVPt      = 0;
   genVPhi     = 0;
   genVy       = 0;
   genVMass    = 0;
   genLepPdgID = 0;
   genLepPt    = 0;
   genLepPhi   = 0;
   u1          = 0;
   u2          = 0;
   const reco::GenParticle* genLep = goodEle.genLepton();
   if(genLep!=NULL) {
     genLepPdgID = genLep->pdgId();
     genLepPt    = genLep->pt();
     genLepPhi   = genLep->phi();
     const reco::Candidate* mother = genLep->mother(0);
     genVPdgID = mother->pdgId();
     genVPt    = mother->pt();
     genVPhi   = mother->phi();
     genVy     = mother->y();
     genVMass  = mother->mass();
     TVector2 vWPt(genVPt*cos(genVPhi),genVPt*sin(genVPhi));
     TVector2 vLepPt(vLep.Px(),vLep.Py());
     TVector2 vU = -1.0*(vtype1pfMET+vLepPt);
     u1 = (vWPt.Px()*vU.Px()+vWPt.Py()*vU.Py())/genVPt; // u1 = (pT . u)/|pT|
     u2 = (vWPt.Px()*vU.Px()-vWPt.Py()*vU.Py())/genVPt; // u1 = (pT x u)/|pT|
   }

   scale1fb = weight;
   mt       = sqrt(2.0*vLep.Pt()*vtype1pfMET.Mod()*(1.0-cos(vLep.Phi()-vtype1pfMET.Phi())));
   q        = goodEle.charge();
   lep      = &vLep;

   ///// electron specific /////
   sc = &vSC;
   pfChIso = goodEle.chargedHadronIso();
   pfGamIso = goodEle.photonIso();
   pfNeuIso = goodEle.neutralHadronIso();
   isLooseEle = goodEle.electronID("eidLoose"); passEleLooseID = Convert_We(isLooseEle);
   isTightEle = goodEle.electronID("eidTight"); passEleTightID = Convert_We(isTightEle);

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
//   gSystem->mkdir(outputDir,kTRUE);
//   gSystem->mkdir(ntupDir,kTRUE);

  //
  // Set up output ntuple
  //

  outFile = new TFile(outFilename,"RECREATE");
  outTree = new TTree("Events","Events");

  outTree->Branch("npv",           &npv,           "npv/I");          // number of primary vertices
  outTree->Branch("genVPt",        &genVPt,        "genVPt/F");       // GEN boson pT (signal MC)
  outTree->Branch("genVPhi",       &genVPhi,       "genVPhi/F");      // GEN boson phi (signal MC)
  outTree->Branch("genVy",         &genVy,         "genVy/F");        // GEN boson rapidity (signal MC)
  outTree->Branch("genVMass",      &genVMass,      "genVMass/F");     // GEN boson mass (signal MC)
  outTree->Branch("genLepPt",      &genLepPt,      "genLepPt/F");     // GEN lepton pT (signal MC)
  outTree->Branch("genLepPhi",     &genLepPhi,     "genLepPhi/F");    // GEN lepton phi (signal MC)
  outTree->Branch("scale1fb",      &scale1fb,      "scale1fb/F");     // event weight per 1/fb (MC)
  outTree->Branch("vtype1pfmet",   "TVector2",     &vtype1pfMET);     // type-1 corrected pf MET
  outTree->Branch("vrawpfmet",     "TVector2",     &vrawpfMET);       // raw pf MET
  outTree->Branch("vgenmet",       "TVector2",     &vgenMET);         // generated MET
  outTree->Branch("mt",            &mt,            "mt/F");           // transverse mass
  outTree->Branch("u1",            &u1,            "u1/F");           // parallel component of recoil
  outTree->Branch("u2",            &u2,            "u2/F");           // perpendicular component of recoil 
  outTree->Branch("q",             &q,             "q/I");            // lepton charge
  outTree->Branch("nEvents",       &nEvents,       "nEvents/I");      // events in MC file
  outTree->Branch("lep", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep);   // lepton 4-vector
  ///// electron specific /////
  outTree->Branch("pfChIso",       &pfChIso,       "pfChIso/F");      // PF charged hadron isolation of electron
  outTree->Branch("pfGamIso",      &pfGamIso,      "pfGamIso/F");     // PF photon isolation of electron
  outTree->Branch("pfNeuIso",      &pfNeuIso,      "pfNeuIso/F");     // PF neutral hadron isolation of electron
  outTree->Branch("isLooseEle",    &isLooseEle,    "isLooseEle/F");   // loose electron ID
  outTree->Branch("isTightEle",    &isTightEle,    "isTightEle/F");   // tight electron ID
}

// ------------ method called once each job just after ending the event loop  ------------
void 
selectWe::endJob() 
{

   // The only information in the last entry of the tree is the total number of events
   // that were processed (not the same as the total number of events selected).
   npv=0;
   genVPdgID=0, genVPt=0, genVPhi=0, genVy=0, genVMass=0;
   genLepPt=0, genLepPhi=0;
   scale1fb=0;
   vtype1pfMET.Set(0.0,0.0);
   vrawpfMET.Set(0.0,0.0);
   vgenMET.Set(0.0,0.0);
   q=0;
   nEvents = dummynEvents;
   LorentzVector dummyLep(0, 0, 0, 0);
   lep = &dummyLep;
   ///// electron specific /////
   pfChIso=0, pfGamIso=0, pfNeuIso=0;
   isLooseEle=0; isTightEle=0;
   passEleLooseID=kFALSE; passEleTightID=kFALSE;
   sc=0;

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
