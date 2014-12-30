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
// https://github.com/jaylawhorn/mitewk/blob/master/Efficiency/selectProbesEleEff.C

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

class selectProbesEleEff : public edm::EDAnalyzer {
   public:
      explicit selectProbesEleEff(const edm::ParameterSet&);
      ~selectProbesEleEff();

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
bool Convert_EleEff(unsigned int val,bool print=kFALSE)
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

const Double_t ELE_MASS   = 0.000511;
const Double_t TAG_PT_CUT = 25;

Double_t nsel_EleEff=0, nselvar_EleEff=0;

TString outFilename_EleEff = TString("selectEleProbes.root");
TFile *outFile_EleEff = new TFile();
TTree *outTree_EleEff = new TTree();

//
// Declare output ntuple variables
//
UInt_t  npv_EleEff;
Float_t scale1fb_EleEff;
Float_t rawpfMETpx_EleEff, rawpfMETpy_EleEff;
TVector2 vtype1pfMET_EleEff, vrawpfMET_EleEff, vgenMET_EleEff;
Float_t mt_EleEff, u1_EleEff, u2_EleEff;
Int_t q1_EleEff, q2_EleEff;
Int_t dummynEvents_EleEff_EleEff=0, nEvents_EleEff=0;
LorentzVector *dilep_EleEff=0, *lep1_EleEff=0, *lep2_EleEff=0;
LorentzVector *sc1_EleEff=0, *sc2_EleEff=0;
Bool_t passEleLooseID1, passEleTightID1;
Bool_t passEleLooseID2, passEleTightID2;

// Compute MC event weight_EleEff_sel per 1/fb
Double_t weight_EleEff = 1;
//const Double_t xsec = 1;
//if(xsec>0) weight_EleEff = 1000.*xsec/(Double_t)eventTree->GetEntries();

//
// constructors and destructor
//
selectProbesEleEff::selectProbesEleEff(const edm::ParameterSet& iConfig):
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
   electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
   metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
   pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands")))
{
   //now do what ever initialization is needed

}


selectProbesEleEff::~selectProbesEleEff()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
selectProbesEleEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;

   dummynEvents_EleEff_EleEff++;

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   // good vertex requirement
   if (vertices->empty()) return; // skip the event if no PV found
//   const reco::Vertex &PV = vertices->front();

   Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);
   if(electrons->size()<2) return;// skip the event if there is no possibility of both tag and probe leptons

   // look for tag lepton
   unsigned int tagIdx=electrons->size();
   Float_t tagPt=0;
   Float_t isTightEle=0;
   Bool_t passEleTightID=kFALSE;
   for (unsigned int jElectron=0;jElectron < electrons->size();jElectron++) {
     const pat::Electron &ele = (*electrons)[jElectron];

     isTightEle = ele.electronID("eidTight"); passEleTightID = Convert_EleEff(isTightEle);

     if(ele.pt()>tagPt && ele.pt()>TAG_PT_CUT && passEleTightID) {
       tagPt = ele.pt();
       tagIdx = jElectron;
     }
   }
   if(tagIdx==electrons->size()) return; // skip event if there is no tag lepton
   const pat::Electron &tag = (*electrons)[tagIdx];

   nsel_EleEff    += weight_EleEff;
   nselvar_EleEff += weight_EleEff*weight_EleEff;

   // look for probe lepton
   unsigned int probeIdx=0;
   Float_t probePt=0;
   for (unsigned int kElectron=0;kElectron < electrons->size();kElectron++) {
     if(kElectron==tagIdx) continue;
     const pat::Electron &ele = (*electrons)[kElectron];
     if(ele.pt()>probePt) {
       probePt = ele.pt();
       probeIdx = kElectron;
     }
   }
   const pat::Electron &probe = (*electrons)[probeIdx];

   LorentzVector vLep1(tag.pt(),tag.eta(),tag.phi(),ELE_MASS);
   LorentzVector vLep2(probe.pt(),probe.eta(),probe.phi(),ELE_MASS);
   LorentzVector vDilep = vLep1 + vLep2;
   LorentzVector vSC1(tag.superCluster()->energy(),tag.superCluster()->eta(),tag.superCluster()->phi(),ELE_MASS);
   LorentzVector vSC2(probe.superCluster()->energy(),probe.superCluster()->eta(),probe.superCluster()->phi(),ELE_MASS);

   // Save MET information

   edm::Handle<pat::METCollection> mets;
   iEvent.getByToken(metToken_, mets);
   const pat::MET &met = mets->front();

   edm::Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);
   rawpfMETpx_EleEff=0;
   rawpfMETpy_EleEff=0;
   // loop on PF candidates to calculate sum of Et's
   for(unsigned int jcand=0; jcand < pfs->size(); jcand++) {
     const pat::PackedCandidate &pf = (*pfs)[jcand];
     rawpfMETpx_EleEff -= pf.px();
     rawpfMETpy_EleEff -= pf.py();
   }

   vtype1pfMET_EleEff.Set(met.px(),met.py());
   vrawpfMET_EleEff.Set(rawpfMETpx_EleEff,rawpfMETpy_EleEff);
   vgenMET_EleEff.Set(met.genMET()->px(),met.genMET()->py());

   //
   // Fill tree
   //
   npv_EleEff = vertices->size();
   u1_EleEff          = 0;
   u2_EleEff          = 0;
   const reco::GenParticle* genLep = tag.genLepton();
   if(genLep!=NULL) {
     const reco::Candidate* mother = genLep->mother(0);
     Float_t genVPt    = mother->pt();
     Float_t genVPhi   = mother->phi();
     TVector2 vVPt(genVPt*cos(genVPhi),genVPt*sin(genVPhi));
     TVector2 vLepPt(vLep1.Px(),vLep1.Py());
     TVector2 vU = -1.0*(vtype1pfMET_EleEff+vLepPt);
     u1_EleEff = (vVPt.Px()*vU.Px()+vVPt.Py()*vU.Py())/genVPt; // u1_EleEff = (pT . u)/|pT|
     u2_EleEff = (vVPt.Px()*vU.Px()-vVPt.Py()*vU.Py())/genVPt; // u1_EleEff = (pT x u)/|pT|
   }

   scale1fb_EleEff = weight_EleEff;
   mt_EleEff       = sqrt(2.0*vLep1.Pt()*vtype1pfMET_EleEff.Mod()*(1.0-cos(vLep1.Phi()-vtype1pfMET_EleEff.Phi())));
   q1_EleEff       = tag.charge();
   q2_EleEff       = probe.charge();
   lep1_EleEff     = &vLep1;
   lep2_EleEff     = &vLep2;
   dilep_EleEff    = &vDilep;
   sc1_EleEff      = &vSC1;
   sc2_EleEff      = &vSC2;
   Float_t isLooseEle1 = tag.electronID("eidLoose"); passEleLooseID1 = Convert_EleEff(isLooseEle1);
   Float_t isTightEle1 = tag.electronID("eidTight"); passEleTightID1 = Convert_EleEff(isTightEle1);
   Float_t isLooseEle2 = probe.electronID("eidLoose"); passEleLooseID2 = Convert_EleEff(isLooseEle2);
   Float_t isTightEle2 = probe.electronID("eidTight"); passEleTightID2 = Convert_EleEff(isTightEle2);

   outTree_EleEff->Fill();

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
selectProbesEleEff::beginJob()
{
   // Create output directory
//   gSystem->mkdir(outputDir,kTRUE);
//   gSystem->mkdir(ntupDir,kTRUE);

  //
  // Set up output ntuple
  //

  outFile_EleEff = new TFile(outFilename_EleEff,"RECREATE");
  outTree_EleEff = new TTree("Events","Events");

  outTree_EleEff->Branch("npv",           &npv_EleEff,           "npv/I");              // number of primary vertices
  outTree_EleEff->Branch("scale1fb",      &scale1fb_EleEff,      "scale1fb/F");         // event weight_EleEff per 1/fb (MC)
  outTree_EleEff->Branch("vtype1pfmet",   "TVector2",             &vtype1pfMET_EleEff); // type-1 corrected pf MET
  outTree_EleEff->Branch("vrawpfmet",     "TVector2",             &vrawpfMET_EleEff);   // raw pf MET
  outTree_EleEff->Branch("vgenmet",       "TVector2",             &vgenMET_EleEff);     // generated MET
  outTree_EleEff->Branch("mt",            &mt_EleEff,            "mt/F");               // transverse mass
  outTree_EleEff->Branch("u1",            &u1_EleEff,            "u1/F");               // parallel component of recoil
  outTree_EleEff->Branch("u2",            &u2_EleEff,            "u2/F");               // perpendicular component of recoil 
  outTree_EleEff->Branch("q1",            &q1_EleEff,            "q1/I");               // tag lepton charge
  outTree_EleEff->Branch("q2",            &q2_EleEff,            "q2/I");               // probe lepton charge
  outTree_EleEff->Branch("nEvents",       &nEvents_EleEff,       "nEvents_EleEff/I");   // events in MC file
  outTree_EleEff->Branch("dilep", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &dilep_EleEff); // dilep_EleEffton 4-vector
  outTree_EleEff->Branch("lep1",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep1_EleEff);  // tag lepton 4-vector
  outTree_EleEff->Branch("lep2",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep2_EleEff);  // probe lepton 4-vector
  outTree_EleEff->Branch("sc1",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sc1_EleEff);   // tag supercluster 4-vector
  outTree_EleEff->Branch("sc2",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sc2_EleEff);   // probe supercluster 4-vector
  outTree_EleEff->Branch("passEleLooseID1",   &passEleLooseID1,   "passEleLooseID1/O"); // tag lepton loose electron ID
  outTree_EleEff->Branch("passEleTightID1",   &passEleTightID1,   "passEleTightID1/O"); // tag lepton tight electron ID
  outTree_EleEff->Branch("passEleLooseID2",   &passEleLooseID2,   "passEleLooseID2/O"); // probe lepton loose electron ID
  outTree_EleEff->Branch("passEleTightID2",   &passEleTightID2,   "passEleTightID2/O"); // probe lepton tight electron ID
}

// ------------ method called once each job just after ending the event loop  ------------
void 
selectProbesEleEff::endJob() 
{

   // The only information in the last entry of the tree is the total number of events
   // that were processed (not the same as the total number of events selected).
   npv_EleEff=0;
   scale1fb_EleEff=0;
   vtype1pfMET_EleEff.Set(0.0,0.0);
   vrawpfMET_EleEff.Set(0.0,0.0);
   vgenMET_EleEff.Set(0.0,0.0);
   u1_EleEff=0; u2_EleEff=0;
   q1_EleEff=0; q2_EleEff=0;
   nEvents_EleEff = dummynEvents_EleEff_EleEff;
   dilep_EleEff=0; lep1_EleEff=0; lep2_EleEff=0;
   sc1_EleEff=0; sc2_EleEff=0;
   passEleLooseID1=kFALSE, passEleTightID1=kFALSE;
   passEleLooseID2=kFALSE, passEleTightID2=kFALSE;

   outTree_EleEff->Fill();

   std::cout << nsel_EleEff << " +/- " << sqrt(nselvar_EleEff) << " per 1/fb" << std::endl;
   std::cout << "endJob: nEvents_EleEff is " << nEvents_EleEff << std::endl;
   outFile_EleEff->Write();
   outFile_EleEff->Close();

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================

  std::cout << "*" << std::endl;
  std::cout << "* SUMMARY" << std::endl;
  std::cout << "*--------------------------------------------------" << std::endl;
  std::cout << "W -> e nu" << std::endl;
  std::cout << " pT > " << TAG_PT_CUT << std::endl;
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "  <> Output saved in " << outFilename_EleEff << "/" << std::endl;
  std::cout << std::endl;
}

// ------------ method called when starting to processes a run  ------------
/*
void 
selectProbesEleEff::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
selectProbesEleEff::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
selectProbesEleEff::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
selectProbesEleEff::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
selectProbesEleEff::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(selectProbesEleEff);
