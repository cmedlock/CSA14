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
// https://github.com/jaylawhorn/mitewk/blob/master/Selection/selectZee.C

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
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TVector2.h>               // 2D vector class
#include <TMath.h>

//
// class declaration
//

class selectZee : public edm::EDAnalyzer {
   public:
      explicit selectZee(const edm::ParameterSet&);
      ~selectZee();

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
      edm::EDGetTokenT<edm::View<reco::GsfElectron>> gsfelectronToken_;
};

//
// constants, enums and typedefs
//

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

// Convert int to binary value
bool Convert_Zee(unsigned int val,bool print=kFALSE)
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

Double_t nsel_Zee=0, nselvar_Zee=0;

TString outFilename_Zee = TString("selectZee.root");
TFile *outFile_Zee = new TFile();
TTree *outTree_Zee = new TTree();

//
// Declare output ntuple variables
//
UInt_t  matchGen_Zee, npv_Zee;
Float_t genVPdgID_Zee, genVPt_Zee, genVPhi_Zee, genVy_Zee, genVMass_Zee;
Float_t scale1fb_Zee;
Float_t rawpfMETpx_Zee, rawpfMETpy_Zee;
TVector2 vtype1pfMET_Zee, vrawpfMET_Zee, vgenMET_Zee;
Float_t mt_Zee, u1_Zee, u2_Zee;
Int_t q1_Zee, q2_Zee;
Int_t dummynEvents_Zee_Zee=0, nEvents_Zee=0;
LorentzVector *dilep_Zee=0, *lep1_Zee=0, *lep2_Zee=0;
LorentzVector *sc1_Zee=0, *sc2_Zee=0;
Bool_t passEleLooseID1, passEleTightID1;
Bool_t passEleLooseID2, passEleTightID2;
Float_t dEtaIn1_Zee, dEtaIn2_Zee, dPhiIn1_Zee, dPhiIn2_Zee;
Float_t full5x5_sIeIe1_Zee, full5x5_sIeIe2_Zee;
Float_t hOverE1_Zee, hOverE2_Zee;
Float_t d0vtx1_Zee, d0vtx2_Zee, dzvtx1_Zee, dzvtx2_Zee;
Float_t EInvOverPInv1_Zee, EInvOverPInv2_Zee;
Float_t relIsoWithDBeta1_Zee, relIsoWithDBeta2_Zee;
Float_t vtxFitProb_Zee;
UInt_t expectedMissingInnerHits1_Zee, expectedMissingInnerHits2_Zee;

// Compute MC event weight_Zee_sel per 1/fb
Double_t weight_Zee = 1;
//const Double_t xsec = 1;
//if(xsec>0) weight_Zee = 1000.*xsec/(Double_t)eventTree->GetEntries();

//
// constructors and destructor
//
selectZee::selectZee(const edm::ParameterSet& iConfig):
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
   electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
   metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
   pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands")))
{
   //now do what ever initialization is needed

}


selectZee::~selectZee()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
selectZee::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;

   dummynEvents_Zee_Zee++;

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   // Good vertex requirement
   if (vertices->empty()) return; // Skip the event if no PV found
   const reco::Vertex &PV = vertices->front();
   vtxFitProb_Zee = TMath::Prob(PV.chi2(),PV.ndof());
   if(vtxFitProb_Zee <= 1e-6) return; // From CSA14 selection

   Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);
   if(electrons->size()<2) return;// Skip the event if there is no possibility of both tag and probe leptons

   // Look for tag lepton
   unsigned int tagIdx=electrons->size();
   Float_t tagPt=0;
   // Float_t isTightEle=0;
   // Bool_t passEleTightID=kFALSE;
   for (unsigned int jElectron=0;jElectron < electrons->size();jElectron++) {
     const pat::Electron &ele = (*electrons)[jElectron];

     // isTightEle = ele.electronID("eidTight"); passEleTightID = Convert_Zee(isTightEle);

     //
     // The calculations below were taken from:
     // https://github.com/lgray/cmssw/blob/common_isolation_selection_70X/TestElectronID/ElectronIDAnalyzer/plugins/ElectronIDAnalyzer.cc
     //
     // ID and matching
     Float_t dEtaIn = ele.deltaEtaSuperClusterTrackAtVtx();
     Float_t dPhiIn = ele.deltaPhiSuperClusterTrackAtVtx();
     Float_t hOverE = ele.hcalOverEcal();
     Float_t full5x5_sIeIe = ele.full5x5_sigmaIetaIeta();
     // Impact parameters
     Float_t d0vtx = (-1)*ele.gsfTrack()->dxy(PV.position());
     Float_t dzvtx = (-1)*ele.gsfTrack()->dz(PV.position());
     // |1/E - 1/p|
     Float_t EInvOverPInv = 0;
     if(ele.ecalEnergy()==0 || !std::isfinite(ele.ecalEnergy())) EInvOverPInv = 1e30;
     else EInvOverPInv = fabs(1.0/ele.ecalEnergy()-ele.eSuperClusterOverP()/ele.ecalEnergy());
     // Isolation
     reco::GsfElectron::PflowIsolationVariables pfIso = ele.pfIsolationVariables();
     Float_t absiso = pfIso.sumChargedHadronPt + std::max(0.0,pfIso.sumNeutralHadronEt+pfIso.sumPhotonEt-0.5*pfIso.sumPUPt);
     Float_t relIsoWithDBeta = absiso/ele.pt();
     // Conversion rejection - vertex fit probability and missing hits
     Float_t expectedMissingInnerHits = ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

    // CSA14 selection - barrel cuts
    if( fabs(ele.superCluster()->eta()) <= 1.479 ) {
      if( fabs(dEtaIn)  >= 0.02   ) continue;
      if( fabs(dPhiIn)  >= 0.2579 ) continue;
      if( full5x5_sIeIe >= 0.0125 ) continue;
      if( hOverE        >= 0.2564 ) continue;
      if( fabs(d0vtx)        >= 0.025  ) continue;
      if( fabs(dzvtx)        >= 0.5863 ) continue;
      if( fabs(EInvOverPInv) >= 0.1508 ) continue;
      if( relIsoWithDBeta    >= 0.3313 ) continue;
      if( expectedMissingInnerHits > 2 ) continue;
      tagPt  = ele.pt();
      tagIdx = jElectron;
    }
    // CSA14 selection - endcap cuts
    else if( 1.479 < fabs(ele.superCluster()->eta()) && fabs(ele.superCluster()->eta()) < 2.5 ) {
      if( fabs(dEtaIn)  >= 0.0141 ) continue;
      if( fabs(dPhiIn)  >= 0.2591 ) continue;
      if( full5x5_sIeIe >= 0.0371 ) continue;
      if( hOverE        >= 0.1335 ) continue;
      if( fabs(d0vtx)        >= 0.2232 ) continue;
      if( fabs(dzvtx)        >= 0.9513 ) continue;
      if( fabs(EInvOverPInv) >= 0.1542 ) continue;
      if( relIsoWithDBeta    >= 0.3816 ) continue;
      if( expectedMissingInnerHits > 3 ) continue;
      tagPt  = ele.pt();
      tagIdx = jElectron;
    }
    else {
      tagPt  = tagPt;
      tagIdx = tagIdx;
    }

   }
   if(tagIdx==electrons->size()) return; // Skip event if there is no tag lepton
   const pat::Electron &tag = (*electrons)[tagIdx];

   nsel_Zee    += weight_Zee;
   nselvar_Zee += weight_Zee*weight_Zee;

   // Look for probe lepton
   unsigned int probeIdx=0;
   Float_t probePt=0;
   for (unsigned int kElectron=0;kElectron < electrons->size();kElectron++) {
     if(kElectron==tagIdx) continue;
     const pat::Electron &ele = (*electrons)[kElectron];
     if(ele.pt()>probePt) {
       probePt  = ele.pt();
       probeIdx = kElectron;
     }
   }
   const pat::Electron &probe = (*electrons)[probeIdx];

   LorentzVector vLep1(tag.pt(),tag.eta(),tag.phi(),ELE_MASS);
   LorentzVector vLep2(probe.pt(),probe.eta(),probe.phi(),ELE_MASS);
   LorentzVector vDilep = vLep1 + vLep2;
   LorentzVector vSC1(tag.superCluster()->energy(),tag.superCluster()->eta(),tag.superCluster()->phi(),ELE_MASS);
   LorentzVector vSC2(probe.superCluster()->energy(),probe.superCluster()->eta(),probe.superCluster()->phi(),ELE_MASS);

   // Perform matching of dileptons to GEN leptons from Z decay
   const reco::GenParticle* gen1 = tag.genParticle();
   const reco::GenParticle* gen2 = probe.genParticle();
   Bool_t hasGenMatch = kFALSE;
   if(gen1!=NULL && gen2!=NULL) {
     Int_t id1 = gen1->pdgId(), id2 = gen2->pdgId();
     Float_t eta1 = gen1->eta(), eta2 = gen2->eta();
     Float_t phi1 = gen1->phi(), phi2 = gen2->phi();
     Bool_t match1 = ( fabs(id1)==11 && sqrt((tag.eta()-eta1)*(tag.eta()-eta1)+(tag.phi()-phi1)*(tag.phi()-phi1)) < 0.5 );
     Bool_t match2 = ( fabs(id2)==11 && sqrt((probe.eta()-eta2)*(probe.eta()-eta2)+(probe.phi()-phi2)*(probe.phi()-phi2)) < 0.5 );
     if(match1 && match2) hasGenMatch = kTRUE;
   }

   // Save MET information

   edm::Handle<pat::METCollection> mets;
   iEvent.getByToken(metToken_, mets);
   const pat::MET &met = mets->front();

   edm::Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);
   rawpfMETpx_Zee=0;
   rawpfMETpy_Zee=0;
   // Loop on PF candidates to calculate sum of Et's
   for(unsigned int jcand=0; jcand < pfs->size(); jcand++) {
     const pat::PackedCandidate &pf = (*pfs)[jcand];
     rawpfMETpx_Zee -= pf.px();
     rawpfMETpy_Zee -= pf.py();
   }

   vtype1pfMET_Zee.Set(met.px(),met.py());
   vrawpfMET_Zee.Set(rawpfMETpx_Zee,rawpfMETpy_Zee);
   vgenMET_Zee.Set(met.genMET()->px(),met.genMET()->py());

   //
   // Fill tree
   //
   matchGen_Zee  = hasGenMatch ? 1 : 0;
   npv_Zee       = vertices->size();
   genVPdgID_Zee = 0;
   genVPt_Zee    = 0;
   genVPhi_Zee   = 0;
   genVy_Zee     = 0;
   genVMass_Zee  = 0;
   u1_Zee        = 0;
   u2_Zee        = 0;
   // Generator level information
   const reco::GenParticle* genLep = tag.genLepton();
   if(genLep!=NULL) {
     const reco::Candidate* mother = genLep->mother(0);
     genVPdgID_Zee = mother->pdgId();
     genVPt_Zee    = mother->pt();
     genVPhi_Zee   = mother->phi();
     genVy_Zee     = mother->y();
     genVMass_Zee  = mother->mass();
     TVector2 vVPt(genVPt_Zee*cos(genVPhi_Zee),genVPt_Zee*sin(genVPhi_Zee));
     TVector2 vLepPt(vLep1.Px(),vLep1.Py());
     TVector2 vU = -1.0*(vtype1pfMET_Zee+vLepPt);
     u1_Zee = (vVPt.Px()*vU.Px()+vVPt.Py()*vU.Py())/genVPt_Zee; // u1_Zee = (pT . u)/|pT|
     u2_Zee = (vVPt.Px()*vU.Px()-vVPt.Py()*vU.Py())/genVPt_Zee; // u1_Zee = (pT x u)/|pT|
   }
   // 4-vectors and pre-implemented ID's
   scale1fb_Zee = weight_Zee;
   mt_Zee       = sqrt(2.0*vLep1.Pt()*vtype1pfMET_Zee.Mod()*(1.0-cos(vLep1.Phi()-vtype1pfMET_Zee.Phi())));
   q1_Zee       = tag.charge();
   q2_Zee       = probe.charge();
   lep1_Zee     = &vLep1;
   lep2_Zee     = &vLep2;
   dilep_Zee    = &vDilep;
   sc1_Zee      = &vSC1;
   sc2_Zee      = &vSC2;
   Float_t isLooseEle1 = tag.electronID("eidLoose");   passEleLooseID1 = Convert_Zee(isLooseEle1);
   Float_t isTightEle1 = tag.electronID("eidTight");   passEleTightID1 = Convert_Zee(isTightEle1);
   Float_t isLooseEle2 = probe.electronID("eidLoose"); passEleLooseID2 = Convert_Zee(isLooseEle2);
   Float_t isTightEle2 = probe.electronID("eidTight"); passEleTightID2 = Convert_Zee(isTightEle2);
   // ID and matching
   dEtaIn1_Zee   = tag.deltaEtaSuperClusterTrackAtVtx();
   dEtaIn2_Zee   = probe.deltaEtaSuperClusterTrackAtVtx();
   dPhiIn1_Zee   = tag.deltaPhiSuperClusterTrackAtVtx();
   dPhiIn2_Zee   = probe.deltaPhiSuperClusterTrackAtVtx();
   hOverE1_Zee   = tag.hcalOverEcal();
   hOverE2_Zee   = probe.hcalOverEcal();
   full5x5_sIeIe1_Zee = tag.full5x5_sigmaIetaIeta();
   full5x5_sIeIe2_Zee = probe.full5x5_sigmaIetaIeta();
   // Impact parameters
   d0vtx1_Zee    = (-1)*tag.gsfTrack()->dxy(PV.position());
   d0vtx2_Zee    = (-1)*probe.gsfTrack()->dxy(PV.position());
   dzvtx1_Zee    = (-1)*tag.gsfTrack()->dz(PV.position());
   dzvtx2_Zee    = (-1)*probe.gsfTrack()->dz(PV.position());
   // |1/E - 1/p|
   if(tag.ecalEnergy()==0 || !std::isfinite(tag.ecalEnergy())) EInvOverPInv1_Zee = 1e30;
   else EInvOverPInv1_Zee = fabs(1.0/tag.ecalEnergy()-tag.eSuperClusterOverP()/tag.ecalEnergy());
   if(probe.ecalEnergy()==0 || !std::isfinite(probe.ecalEnergy())) EInvOverPInv2_Zee = 1e30;
   else EInvOverPInv2_Zee = fabs(1.0/probe.ecalEnergy()-probe.eSuperClusterOverP()/probe.ecalEnergy());
   // Isolation
   reco::GsfElectron::PflowIsolationVariables pfIso1 = tag.pfIsolationVariables();
   Float_t absiso1 = pfIso1.sumChargedHadronPt + std::max(0.0,pfIso1.sumNeutralHadronEt+pfIso1.sumPhotonEt-0.5*pfIso1.sumPUPt);
   relIsoWithDBeta1_Zee = absiso1/tag.pt();
   reco::GsfElectron::PflowIsolationVariables pfIso2 = probe.pfIsolationVariables();
   Float_t absiso2 = pfIso2.sumChargedHadronPt + std::max(0.0,pfIso2.sumNeutralHadronEt+pfIso2.sumPhotonEt-0.5*pfIso2.sumPUPt);
   relIsoWithDBeta2_Zee = absiso2/probe.pt();
   // Conversion rejection - vertex fit probability and missing hits
   vtxFitProb_Zee = TMath::Prob(PV.chi2(),PV.ndof());
   expectedMissingInnerHits1_Zee = tag.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
   expectedMissingInnerHits2_Zee = probe.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

   outTree_Zee->Fill();

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
selectZee::beginJob()
{
   // Create output directory
//   gSystem->mkdir(outputDir,kTRUE);
//   gSystem->mkdir(ntupDir,kTRUE);

  //
  // Set up output ntuple
  //

  outFile_Zee = new TFile(outFilename_Zee,"RECREATE");
  outTree_Zee = new TTree("Events","Events");

  outTree_Zee->Branch("matchGen",      &matchGen_Zee,      "matchGen/i");         // event has both leptons matched to MC Z->ll
  outTree_Zee->Branch("npv",           &npv_Zee,           "npv/I");              // number of primary vertices
  outTree_Zee->Branch("genVPt",        &genVPt_Zee,        "genVPt/F");           // GEN boson pT (signal MC)
  outTree_Zee->Branch("genVPhi",       &genVPhi_Zee,       "genVPhi/F");          // GEN boson phi (signal MC)
  outTree_Zee->Branch("genVy",         &genVy_Zee,         "genVy/F");            // GEN boson rapidity (signal MC)
  outTree_Zee->Branch("genVMass",      &genVMass_Zee,      "genVMass/F");         // GEN boson mass (signal MC)
  outTree_Zee->Branch("scale1fb",      &scale1fb_Zee,      "scale1fb/F");         // event weight_Zee per 1/fb (MC)
  outTree_Zee->Branch("vtype1pfmet",   "TVector2",         &vtype1pfMET_Zee);     // type-1 corrected pf MET
  outTree_Zee->Branch("vrawpfmet",     "TVector2",         &vrawpfMET_Zee);       // raw pf MET
  outTree_Zee->Branch("vgenmet",       "TVector2",         &vgenMET_Zee);         // generated MET
  outTree_Zee->Branch("mt",            &mt_Zee,            "mt/F");               // transverse mass
  outTree_Zee->Branch("u1",            &u1_Zee,            "u1/F");               // parallel component of recoil
  outTree_Zee->Branch("u2",            &u2_Zee,            "u2/F");               // perpendicular component of recoil 
  outTree_Zee->Branch("q1",            &q1_Zee,            "q1/I");               // tag lepton charge
  outTree_Zee->Branch("q2",            &q2_Zee,            "q2/I");               // probe lepton charge
  outTree_Zee->Branch("nEvents",       &nEvents_Zee,       "nEvents_Zee/I");   // events in MC file
  outTree_Zee->Branch("dilep", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &dilep_Zee); // dilep_Zeeton 4-vector
  outTree_Zee->Branch("lep1",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep1_Zee);  // tag lepton 4-vector
  outTree_Zee->Branch("lep2",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep2_Zee);  // probe lepton 4-vector
  outTree_Zee->Branch("sc1",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sc1_Zee);   // tag supercluster 4-vector
  outTree_Zee->Branch("sc2",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sc2_Zee);   // probe supercluster 4-vector
  outTree_Zee->Branch("passEleLooseID1",  &passEleLooseID1,      "passEleLooseID1/O");  // tag lepton loose electron ID
  outTree_Zee->Branch("passEleTightID1",  &passEleTightID1,      "passEleTightID1/O");  // tag lepton tight electron ID
  outTree_Zee->Branch("passEleLooseID2",  &passEleLooseID2,      "passEleLooseID2/O");  // probe lepton loose electron ID
  outTree_Zee->Branch("passEleTightID2",  &passEleTightID2,      "passEleTightID2/O");  // probe lepton tight electron ID
  outTree_Zee->Branch("dEtaIn1",          &dEtaIn1_Zee,          "dEtaIn1/F");
  outTree_Zee->Branch("dEtaIn2",          &dEtaIn2_Zee,          "dEtaIn2/F");
  outTree_Zee->Branch("dPhiIn1",          &dPhiIn1_Zee,          "dPhiIn1/F");
  outTree_Zee->Branch("dPhiIn2",          &dPhiIn2_Zee,          "dPhiIn2/F");
  outTree_Zee->Branch("full5x5_sIeIe1",   &full5x5_sIeIe1_Zee,   "full5x5_sIeIe1/F");
  outTree_Zee->Branch("full5x5_sIeIe2",   &full5x5_sIeIe2_Zee,   "full5x5_sIeIe2/F");
  outTree_Zee->Branch("hOverE1",          &hOverE1_Zee,          "hOverE1/F");
  outTree_Zee->Branch("hOverE2",          &hOverE2_Zee,          "hOverE2/F");
  outTree_Zee->Branch("d0vtx1",           &d0vtx1_Zee,           "d0vtx1/F");
  outTree_Zee->Branch("d0vtx2",           &d0vtx2_Zee,           "d0vtx2/F");
  outTree_Zee->Branch("dzvtx1",           &dzvtx1_Zee,           "dzvtx1/F");
  outTree_Zee->Branch("dzvtx2",           &dzvtx2_Zee,           "dzvtx2/F");
  outTree_Zee->Branch("EInvOverPInv1",    &EInvOverPInv1_Zee,    "EInvOverPInv1/F");
  outTree_Zee->Branch("EInvOverPInv2",    &EInvOverPInv2_Zee,    "EInvOverPInv2/F");
  outTree_Zee->Branch("relIsoWithDBeta1", &relIsoWithDBeta1_Zee, "relIsoWithDBeta1/F");
  outTree_Zee->Branch("relIsoWithDBeta2", &relIsoWithDBeta2_Zee, "relIsoWithDBeta2/F");
  outTree_Zee->Branch("vtxFitProb",       &vtxFitProb_Zee,       "vtxFitProb/F");
  outTree_Zee->Branch("expectedMissingInnerHits1", &expectedMissingInnerHits1_Zee, "expectedMissingInnerHits1/I");
  outTree_Zee->Branch("expectedMissingInnerHits2", &expectedMissingInnerHits2_Zee, "expectedMissingInnerHits2/I");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
selectZee::endJob() 
{

   // The only information in the last entry of the tree is the total number of events
   // that were processed (not the same as the total number of events selected).
   npv_Zee=0;
   scale1fb_Zee=0;
   vtype1pfMET_Zee.Set(0.0,0.0);
   vrawpfMET_Zee.Set(0.0,0.0);
   vgenMET_Zee.Set(0.0,0.0);
   u1_Zee=0; u2_Zee=0;
   q1_Zee=0; q2_Zee=0;
   nEvents_Zee = dummynEvents_Zee_Zee;
   dilep_Zee=0; lep1_Zee=0; lep2_Zee=0;
   sc1_Zee=0; sc2_Zee=0;
   passEleLooseID1=kFALSE, passEleTightID1=kFALSE;
   passEleLooseID2=kFALSE, passEleTightID2=kFALSE;

   outTree_Zee->Fill();

   std::cout << nsel_Zee << " +/- " << sqrt(nselvar_Zee) << " per 1/fb" << std::endl;
   std::cout << "endJob: nEvents_Zee is " << nEvents_Zee << std::endl;
   outFile_Zee->Write();
   outFile_Zee->Close();

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================

  std::cout << "*" << std::endl;
  std::cout << "* SUMMARY" << std::endl;
  std::cout << "*--------------------------------------------------" << std::endl;
  std::cout << "Z -> e e" << std::endl;
  std::cout << " pT > " << TAG_PT_CUT << std::endl;
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "  <> Output saved in " << outFilename_Zee << "/" << std::endl;
  std::cout << std::endl;
}

// ------------ method called when starting to processes a run  ------------
/*
void 
selectZee::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
selectZee::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
selectZee::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
selectZee::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
selectZee::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(selectZee);
