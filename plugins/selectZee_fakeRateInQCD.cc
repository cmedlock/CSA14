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
// https://github.com/jaylawhorn/mitewk/blob/master/Selection/selectZee_fakeRateInQCD.C

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
#include <TMath.h>

//
// class declaration
//

class selectZee_fakeRateInQCD : public edm::EDAnalyzer {
   public:
      explicit selectZee_fakeRateInQCD(const edm::ParameterSet&);
      ~selectZee_fakeRateInQCD();

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
bool Convert_fake(unsigned int val,bool print=kFALSE)
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

TString outFilename_fake = TString("selectZee_fakeRateInQCD.root");
TFile *outFile_fake = new TFile();
TTree *outTree_fake = new TTree();

//
// Declare output ntuple variables
//
std::vector<Float_t> elePts;
std::vector<Float_t> eleEtas;
std::vector<Float_t> isVetoEle;
std::vector<Float_t> isLooseEle;
std::vector<Float_t> isMediumEle;
std::vector<Float_t> isTightEle;

std::vector<Float_t> *v_elePts      = &elePts;
std::vector<Float_t> *v_eleEtas     = &eleEtas;
std::vector<Float_t> *v_isVetoEle   = &isVetoEle;
std::vector<Float_t> *v_isLooseEle  = &isLooseEle;
std::vector<Float_t> *v_isMediumEle = &isMediumEle;
std::vector<Float_t> *v_isTightEle  = &isTightEle;

//
// constructors and destructor
//
selectZee_fakeRateInQCD::selectZee_fakeRateInQCD(const edm::ParameterSet& iConfig):
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
   electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
   metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
   pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands")))
{
   //now do what ever initialization is needed

}


selectZee_fakeRateInQCD::~selectZee_fakeRateInQCD()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
selectZee_fakeRateInQCD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   // Good vertex requirement
   if (vertices->empty()) return; // Skip the event if no PV found

   Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);
   if(electrons->size()==0) return;

   // Save electron information
   for (unsigned int jElectron=0;jElectron < electrons->size();jElectron++) {
     const pat::Electron &ele = (*electrons)[jElectron];

     // Perform matching of lepton to GEN lepton
     const reco::GenParticle* gen = ele.genParticle();
     if(gen==NULL) continue;
     Int_t   id    = gen->pdgId();
     Float_t eta   = gen->eta();
     Float_t phi   = gen->phi();
     Bool_t  match = ( fabs(id)==11 && sqrt((ele.eta()-eta)*(ele.eta()-eta)+(ele.phi()-phi)*(ele.phi()-phi)) < 0.5 );
     if(!match) continue;

     elePts.push_back(ele.pt());
     eleEtas.push_back(ele.eta());
     isVetoEle.push_back(ele.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-veto"));
     isLooseEle.push_back(ele.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-loose"));
     isMediumEle.push_back(ele.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-medium"));
     isTightEle.push_back(ele.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-tight"));

   }

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
selectZee_fakeRateInQCD::beginJob()
{

  //
  // Set up output ntuple
  //

  outFile_fake = new TFile(outFilename_fake,"RECREATE");
  outTree_fake = new TTree("Events","Events");

  outTree_fake->Branch("elePts",      "std::vector<float>", &v_elePts);
  outTree_fake->Branch("eleEtas",     "std::vector<float>", &v_eleEtas);
  outTree_fake->Branch("isVetoEle",   "std::vector<float>", &v_isVetoEle);
  outTree_fake->Branch("isLooseEle",  "std::vector<float>", &v_isLooseEle);
  outTree_fake->Branch("isMediumEle", "std::vector<float>", &v_isMediumEle);
  outTree_fake->Branch("isTightEle",  "std::vector<float>", &v_isTightEle);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
selectZee_fakeRateInQCD::endJob() 
{
   outTree_fake->Fill();
   outFile_fake->Write();
   outFile_fake->Close();

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================

  std::cout << "*" << std::endl;
  std::cout << "* SUMMARY" << std::endl;
  std::cout << "*--------------------------------------------------" << std::endl;
  std::cout << "QCD" << std::endl;
  std::cout << std::endl;

}

// ------------ method called when starting to processes a run  ------------
/*
void 
selectZee_fakeRateInQCD::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
selectZee_fakeRateInQCD::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
selectZee_fakeRateInQCD::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
selectZee_fakeRateInQCD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
selectZee_fakeRateInQCD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(selectZee_fakeRateInQCD);
