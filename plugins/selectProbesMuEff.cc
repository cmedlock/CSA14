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
// https://github.com/jaylawhorn/mitewk/blob/master/Efficiency/selectProbesMuEff.C

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

class selectProbesMuEff : public edm::EDAnalyzer {
   public:
      explicit selectProbesMuEff(const edm::ParameterSet&);
      ~selectProbesMuEff();

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
bool Convert_MuEff(unsigned int val,bool print=kFALSE)
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

const Double_t TAG_PT_CUT = 25;
const Double_t MUON_MASS = 0.105658369;

Double_t nsel_MuEff=0, nselvar_MuEff=0;

TString outFilename_MuEff = TString("selectMuonProbes.root");
TFile *outFile_MuEff = new TFile();
TTree *outTree_MuEff = new TTree();

//
// Declare output ntuple variables
//
UInt_t  npv_MuEff;
Float_t scale1fb_MuEff;
Float_t rawpfMETpx_MuEff, rawpfMETpy_MuEff;
TVector2 vtype1pfMET_MuEff, vrawpfMET_MuEff, vgenMET_MuEff;
Float_t mt_MuEff, u1_MuEff, u2_MuEff;
Int_t q1_MuEff, q2_MuEff;
Int_t dummynEvents_MuEff=0, nEvents_MuEff=0;
LorentzVector *dilep_MuEff=0, *lep1_MuEff=0, *lep2_MuEff=0;
//LorentzVector *sc1_MuEff=0, *sc2_MuEff=0;
Bool_t passMuonLooseID1_MuEff, passMuonSoftID1_MuEff, passMuonTightID1_MuEff;
Bool_t passMuonLooseID2_MuEff, passMuonSoftID2_MuEff, passMuonTightID2_MuEff;

// Compute MC event weightWm_sel per 1/fb
Double_t weight_MuEff = 1;
//const Double_t xsec = 1;
//if(xsec>0) weightWm = 1000.*xsec/(Double_t)eventTree->GetEntries();

//
// constructors and destructor
//
selectProbesMuEff::selectProbesMuEff(const edm::ParameterSet& iConfig):
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
   muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
   metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
   pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands")))
{
   //now do what ever initialization is needed

}


selectProbesMuEff::~selectProbesMuEff()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
selectProbesMuEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   dummynEvents_MuEff++;

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   // good vertex requirement
   if (vertices->empty()) return; // skip the event if no PV found
   const reco::Vertex &PV = vertices->front();
   npv_MuEff = vertices->size();

   Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);
   if(muons->size()<2) return; // skip the event if there is no possibility of both tag and probe leptons

   // look for tag lepton
   unsigned int tagIdx=muons->size();
   Float_t tagPt=0;
   Float_t isTightMuon=0;
   Bool_t passMuonTightID=kFALSE;
   for (unsigned int jMuon=0;jMuon < muons->size();jMuon++) {
     const pat::Muon &mu = (*muons)[jMuon];

     isTightMuon = mu.isTightMuon(PV); passMuonTightID = Convert_MuEff(isTightMuon);

     if(mu.pt()>tagPt && mu.pt()>TAG_PT_CUT && passMuonTightID) {
       tagPt = mu.pt();
       tagIdx = jMuon;
     }
   }
   if(tagIdx==muons->size()) return; // skip event if there is no tag lepton
   const pat::Muon &tag = (*muons)[tagIdx];

   nsel_MuEff    += weight_MuEff;
   nselvar_MuEff += weight_MuEff*weight_MuEff;

   // look for probe lepton
   unsigned int probeIdx=0;
   Float_t probePt=0;
   for (unsigned int kMuon=0;kMuon < muons->size();kMuon++) {
     if(kMuon==tagIdx) continue;
     const pat::Muon &mu = (*muons)[kMuon];
     if(mu.pt()>probePt) {
       probePt = mu.pt();
       probeIdx = kMuon;
     }
   }
   const pat::Muon &probe = (*muons)[probeIdx];

   LorentzVector vLep1(tag.pt(),tag.eta(),tag.phi(),MUON_MASS);
   LorentzVector vLep2(probe.pt(),probe.eta(),probe.phi(),MUON_MASS);
   LorentzVector vDilep = vLep1 + vLep2;
//   LorentzVector vSC1(tag.superCluster()->energy(),tag.superCluster()->eta(),tag.superCluster()->phi(),MUON_MASS);
//   LorentzVector vSC2(probe.superCluster()->energy(),probe.superCluster()->eta(),probe.superCluster()->phi(),MUON_MASS);

   // Save MET information

   edm::Handle<pat::METCollection> mets;
   iEvent.getByToken(metToken_, mets);
   const pat::MET &met = mets->front();

   edm::Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);
   rawpfMETpx_MuEff=0;
   rawpfMETpy_MuEff=0;
   // loop on pf candidates to calculate sum of Et's
   for(unsigned int jcand=0; jcand < pfs->size(); jcand++) {
     const pat::PackedCandidate &pf = (*pfs)[jcand];
     rawpfMETpx_MuEff -= pf.px();
     rawpfMETpy_MuEff -= pf.py();
   }

   vtype1pfMET_MuEff.Set(met.px(),met.py());
   vrawpfMET_MuEff.Set(rawpfMETpx_MuEff,rawpfMETpy_MuEff);
   vgenMET_MuEff.Set(met.genMET()->px(),met.genMET()->py());

   //
   // Fill tree
   //
   npv_MuEff = vertices->size();
   u1_MuEff  = 0;
   u2_MuEff  = 0;
   const reco::GenParticle* genLep = tag.genLepton();
   if(genLep!=NULL) {
     const reco::Candidate* mother = genLep->mother(0);
     Float_t genVPt    = mother->pt();
     Float_t genVPhi   = mother->phi();
     TVector2 vVPt(genVPt*cos(genVPhi),genVPt*sin(genVPhi));
     TVector2 vLepPt(vLep1.Px(),vLep1.Py());
     TVector2 vU = -1.0*(vtype1pfMET_MuEff+vLepPt);
     u1_MuEff = (vVPt.Px()*vU.Px()+vVPt.Py()*vU.Py())/genVPt; // u1_MuEff = (pT . u)/|pT|
     u2_MuEff = (vVPt.Px()*vU.Px()-vVPt.Py()*vU.Py())/genVPt; // u1_MuEff = (pT x u)/|pT|
   }

   scale1fb_MuEff = weight_MuEff;
   mt_MuEff       = sqrt(2.0*vLep1.Pt()*vtype1pfMET_MuEff.Mod()*(1.0-cos(vLep1.Phi()-vtype1pfMET_MuEff.Phi())));
   q1_MuEff       = tag.charge();
   q2_MuEff       = probe.charge();
   lep1_MuEff     = &vLep1;
   lep2_MuEff     = &vLep2;
   dilep_MuEff    = &vDilep;
//   sc1_MuEff      = &vSC1;
//   sc2_MuEff      = &vSC2;
   Float_t isLooseMuon1 = tag.isLooseMuon();     passMuonLooseID1_MuEff = Convert_MuEff(isLooseMuon1);
   Float_t isSoftMuon1  = tag.isSoftMuon(PV);    passMuonSoftID1_MuEff  = Convert_MuEff(isSoftMuon1);
   Float_t isTightMuon1 = tag.isTightMuon(PV);   passMuonTightID1_MuEff = Convert_MuEff(isTightMuon1);
   Float_t isLooseMuon2 = probe.isLooseMuon();   passMuonLooseID2_MuEff = Convert_MuEff(isLooseMuon2);
   Float_t isSoftMuon2  = probe.isSoftMuon(PV);  passMuonSoftID2_MuEff  = Convert_MuEff(isSoftMuon2);
   Float_t isTightMuon2 = probe.isTightMuon(PV); passMuonTightID2_MuEff = Convert_MuEff(isTightMuon2);

   outTree_MuEff->Fill();

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
selectProbesMuEff::beginJob()
{
   // Create output directory
//   gSystem->mkdir(outputDir,kTRUE);
//   gSystem->mkdir(ntupDir,kTRUE);

  //
  // Set up output ntuple
  //

  outFile_MuEff = new TFile(outFilename_MuEff,"RECREATE");
  outTree_MuEff = new TTree("Events","Events");

  outTree_MuEff->Branch("npv",           &npv_MuEff,           "npv/I");                  // number of primary vertices
  outTree_MuEff->Branch("scale1fb",      &scale1fb_MuEff,      "scale1fb/F");             // event weight_MuEff per 1/fb (MC)
  outTree_MuEff->Branch("vtype1pfmet",   "TVector2",             &vtype1pfMET_MuEff);     // type-1 corrected pf MET
  outTree_MuEff->Branch("vrawpfmet",     "TVector2",             &vrawpfMET_MuEff);       // raw pf MET
  outTree_MuEff->Branch("vgenmet",       "TVector2",             &vgenMET_MuEff);         // generated MET
  outTree_MuEff->Branch("mt",            &mt_MuEff,            "mt/F");                   // transverse mass
  outTree_MuEff->Branch("u1",            &u1_MuEff,            "u1/F");                   // parallel component of recoil
  outTree_MuEff->Branch("u2",            &u2_MuEff,            "u2/F");                   // perpendicular component of recoil 
  outTree_MuEff->Branch("q1",            &q1_MuEff,            "q1/I");                   // tag lepton charge
  outTree_MuEff->Branch("q2",            &q2_MuEff,            "q2/I");                   // probe lepton charge
  outTree_MuEff->Branch("nEvents",       &nEvents_MuEff,       "nEvents_MuEff/I");        // events in MC file
  outTree_MuEff->Branch("dilep", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &dilep_MuEff); // dilepton 4-vector
  outTree_MuEff->Branch("lep1",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep1_MuEff);  // tag lepton 4-vector
  outTree_MuEff->Branch("lep2",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep2_MuEff);  // probe lepton 4-vector
//  outTree_MuEff->Branch("sc1",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sc1_MuEff);   // tag supercluster 4-vector
//  outTree_MuEff->Branch("sc2",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sc2_MuEff);   // probe supercluster 4-vector
  outTree_MuEff->Branch("passMuonLooseID1",   &passMuonLooseID1_MuEff,   "passMuonLooseID1/O");  // tag lepton loose muon ID
  outTree_MuEff->Branch("passMuonSoftID1",    &passMuonSoftID1_MuEff,    "passMuonSoftID1/O");   // tag lepton soft muon ID
  outTree_MuEff->Branch("passMuonTightID1",   &passMuonTightID1_MuEff,   "passMuonTightID/O");   // tag lepton tight muon ID
  outTree_MuEff->Branch("passMuonLooseID2",   &passMuonLooseID2_MuEff,   "passMuonLooseID2/O");  // probe lepton loose muon ID
  outTree_MuEff->Branch("passMuonSoftID2",    &passMuonSoftID2_MuEff,    "passMuonSoftID2/O");   // probe lepton soft muon ID
  outTree_MuEff->Branch("passMuonTightID2",   &passMuonTightID2_MuEff,   "passMuonTightID2/O");  // probe lepton tight muon ID

}

// ------------ method called once each job just after ending the event loop  ------------
void 
selectProbesMuEff::endJob() 
{

   // The only information in the last entry of the tree is the total number of events
   // that were processed (not the same as the total number of events selected).
   npv_MuEff=0;
   scale1fb_MuEff=0;
   vtype1pfMET_MuEff.Set(0.0,0.0);
   vrawpfMET_MuEff.Set(0.0,0.0);
   vgenMET_MuEff.Set(0.0,0.0);
   u1_MuEff=0; u2_MuEff=0;
   q1_MuEff=0; q2_MuEff=0;
   nEvents_MuEff = dummynEvents_MuEff;
   dilep_MuEff=0; lep1_MuEff=0; lep2_MuEff=0;
//   sc1_MuEff=0; sc2_MuEff=0;
   passMuonLooseID1_MuEff=kFALSE, passMuonSoftID1_MuEff=kFALSE; passMuonTightID1_MuEff=kFALSE;
   passMuonLooseID2_MuEff=kFALSE, passMuonSoftID2_MuEff=kFALSE; passMuonTightID2_MuEff=kFALSE;

   outTree_MuEff->Fill();

   std::cout << nsel_MuEff << " +/- " << sqrt(nselvar_MuEff) << " per 1/fb" << std::endl;
   std::cout << "endJob: nEvents_MuEff is " << nEvents_MuEff << std::endl;
   outFile_MuEff->Write();
   outFile_MuEff->Close();

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================

  std::cout << "*" << std::endl;
  std::cout << "* SUMMARY" << std::endl;
  std::cout << "*--------------------------------------------------" << std::endl;
  std::cout << "W -> mu nu" << std::endl;
  std::cout << " pT > " << TAG_PT_CUT << std::endl;
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "  <> Output saved in " << outFilename_MuEff << "/" << std::endl;
  std::cout << std::endl;
}

// ------------ method called when starting to processes a run  ------------
/*
void 
selectProbesMuEff::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
selectProbesMuEff::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
selectProbesMuEff::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
selectProbesMuEff::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
selectProbesMuEff::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(selectProbesMuEff);
