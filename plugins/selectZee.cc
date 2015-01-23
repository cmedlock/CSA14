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
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

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

const Double_t ELE_MASS = 0.000511;

Double_t nsel_Zee=0, eventCounter_Zee=0;

TString outFilename_Zee = TString("selectZee.root");
TFile *outFile_Zee = new TFile();
TTree *outTree_Zee = new TTree();

//
// Declare output ntuple variables
//
Int_t  matchGen_Zee=0, npv_Zee=0, nEvents_Zee=0;
Float_t genVPdgID_Zee=0, genVPt_Zee=0, genVPhi_Zee=0, genVy_Zee=0, genVMass_Zee=0;
Float_t rawpfmet_Zee=0, rawpfmetPhi_Zee=0;
Float_t type1pfmet_Zee=0, type1pfmetPhi_Zee=0;
Float_t genmet_Zee=0, genmetPhi_Zee=0;
Float_t u1_Zee=0, u2_Zee=0;
Int_t q1_Zee=0, q2_Zee=0;
Float_t pfChIso1_Zee=0, pfChIso2_Zee=0;
LorentzVector *dilep_Zee=0, *lep1_Zee=0, *lep2_Zee=0;
LorentzVector *sc1_Zee=0, *sc2_Zee=0;
Int_t isVetoEle1=0, isLooseEle1=0, isMediumEle1=0, isTightEle1=0;
Int_t isVetoEle2=0, isLooseEle2=0, isMediumEle2=0, isTightEle2=0;

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

   // Count the total number of MC events processed
   eventCounter_Zee++;

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   // Good vertex requirement
   if (vertices->empty()) return; // Skip the event if no PV found
   npv_Zee = vertices->size();

   Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);
   if(electrons->size()<2) return; // Skip the event if there is no possibility of both tag and probe leptons

   //
   // SELECTION PROCEDURE:
   //  (1) Find a good electron that passes the CSA14 veto ID -> this will be the "tag"
   //  (2) Pair the tag with a probe electron
   //
   Bool_t foundTag = kFALSE;
   for (unsigned int i1=0;i1 < electrons->size();i1++) {
     const pat::Electron &tag = (*electrons)[i1];

     Int_t isVetoEle = tag.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-veto");
     if(!isVetoEle) continue; // lepton selection

     foundTag = kTRUE;

     // Tag lepton information
     LorentzVector vTag(tag.pt(),tag.eta(),tag.phi(),ELE_MASS);
     LorentzVector vTagSC(tag.superCluster()->energy()*(tag.pt()/tag.p()),tag.superCluster()->eta(),tag.superCluster()->phi(),ELE_MASS);

     lep1_Zee = &vTag;
     sc1_Zee  = &vTagSC;

     q1_Zee       = tag.charge();
     pfChIso1_Zee = tag.chargedHadronIso();
     isVetoEle1   = tag.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-veto");
     isLooseEle1  = tag.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-loose");
     isMediumEle1 = tag.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-medium");
     isTightEle1  = tag.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-tight");

     for (unsigned int i2=0;i2 < electrons->size();i2++) {
       if(i2==i1) continue;
       const pat::Electron &probe = (*electrons)[i2];

       // Probe lepton information
       LorentzVector vProbe(probe.pt(),probe.eta(),probe.phi(),ELE_MASS);
       LorentzVector vProbeSC(probe.superCluster()->energy(),probe.superCluster()->eta(),probe.superCluster()->phi(),ELE_MASS);

       lep2_Zee = &vProbe;
       sc2_Zee  = &vProbeSC;

       q2_Zee       = probe.charge();
       pfChIso2_Zee = probe.chargedHadronIso();
       isVetoEle2   = probe.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-veto");
       isLooseEle2  = probe.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-loose");
       isMediumEle2 = probe.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-medium");
       isTightEle2  = probe.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-tight");

       // Dilepton 4-vector
       LorentzVector vDilep = vTag + vProbe;
       dilep_Zee = &vDilep;

       // Perform matching of dileptons to geneverator level leptons
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
       matchGen_Zee  = hasGenMatch ? 1 : 0;

       // Type-1 corrected PF MET
       edm::Handle<pat::METCollection> mets;
       iEvent.getByToken(metToken_, mets);

       const pat::MET &met = mets->front();
       TVector2 vtype1pfMET_Zee = TVector2(met.px(),met.py());
       type1pfmet_Zee    = vtype1pfMET_Zee.Mod();
       type1pfmetPhi_Zee = vtype1pfMET_Zee.Phi();

       // Generator level MET
       TVector2 vgenMET_Zee = TVector2(met.genMET()->px(),met.genMET()->py());
       genmet_Zee    = vgenMET_Zee.Mod();
       genmetPhi_Zee = vgenMET_Zee.Phi();

       // Raw PF MET
       edm::Handle<pat::PackedCandidateCollection> pfs;
       iEvent.getByToken(pfToken_, pfs);
       Float_t rawpfMETpx=0, rawpfMETpy=0;

       for(unsigned int jcand=0; jcand < pfs->size(); jcand++) {
         const pat::PackedCandidate &pf = (*pfs)[jcand];
         rawpfMETpx -= pf.px();
         rawpfMETpy -= pf.py();
       }

       TVector2 vrawpfMET_Zee = TVector2(rawpfMETpx,rawpfMETpy);
       rawpfmet_Zee    = vrawpfMET_Zee.Mod();
       rawpfmetPhi_Zee = vrawpfMET_Zee.Phi();

       // Generator level lepton information and hadronic recoil
       genVPdgID_Zee = 0;
       genVPt_Zee    = 0;
       genVPhi_Zee   = 0;
       genVy_Zee     = 0;
       genVMass_Zee  = 0;
       u1_Zee        = 0;
       u2_Zee        = 0;
       if(gen1!=NULL) {
         const reco::Candidate* mother = gen1->mother(0);
         genVPdgID_Zee = mother->pdgId();
         genVPt_Zee    = mother->pt();
         genVPhi_Zee   = mother->phi();
         genVy_Zee     = mother->y();
         genVMass_Zee  = mother->mass();	
         TVector2 vZPt((vDilep.Pt())*cos(vDilep.Phi()),(vDilep.Pt())*sin(vDilep.Phi()));
         TVector2 vMet(type1pfmet_Zee*cos(type1pfmetPhi_Zee), type1pfmet_Zee*sin(type1pfmetPhi_Zee));
         TVector2 vU = -1.0*(vMet+vZPt);
         u1_Zee = ((vDilep.Px())*(vU.Px()) + (vDilep.Py())*(vU.Py()))/(vDilep.Pt()); // u1 = (pT . u)/|pT|
         u2_Zee = ((vDilep.Px())*(vU.Py()) - (vDilep.Py())*(vU.Px()))/(vDilep.Pt()); // u2 = (pT x u)/|pT|
       }
       // Count total number of events selected
       nsel_Zee++;
       // Fill tree
       outTree_Zee->Fill();
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
selectZee::beginJob()
{

  //
  // Set up output ntuple
  //

  outFile_Zee = new TFile(outFilename_Zee,"RECREATE");
  outTree_Zee = new TTree("Events","Events");

  outTree_Zee->Branch("matchGen",      &matchGen_Zee,      "matchGen/I");         // event has both leptons matched to MC Z->ll
  outTree_Zee->Branch("npv",           &npv_Zee,           "npv/I");              // number of primary vertices
  outTree_Zee->Branch("nEvents",       &nEvents_Zee,       "nEvents/I");          // total number of MC events processed
                                                                                  // (will be used for event weighting in signal extraction)
  outTree_Zee->Branch("genVPt",        &genVPt_Zee,        "genVPt/F");           // GEN boson pT (signal MC)
  outTree_Zee->Branch("genVPhi",       &genVPhi_Zee,       "genVPhi/F");          // GEN boson phi (signal MC)
  outTree_Zee->Branch("genVy",         &genVy_Zee,         "genVy/F");            // GEN boson rapidity (signal MC)
  outTree_Zee->Branch("genVMass",      &genVMass_Zee,      "genVMass/F");         // GEN boson mass (signal MC)
  outTree_Zee->Branch("rawpfmet",      &rawpfmet_Zee,      "rawpfmet/F");         // Raw PF MET
  outTree_Zee->Branch("rawpfmetPhi",   &rawpfmetPhi_Zee,   "rawpfmetPhi/F");      // Raw PF MET phi
  outTree_Zee->Branch("type1pfmet",    &type1pfmet_Zee,    "type1pfmet/F");       // Type-1 corrected PF MET
  outTree_Zee->Branch("type1pfmetPhi", &type1pfmetPhi_Zee, "type1pfmetPhi/F");    // Type-1 corrected PF MET phi
  outTree_Zee->Branch("genmet",        &genmet_Zee,        "genmet/F");           // Generator level MET
  outTree_Zee->Branch("genmetPhi",     &genmetPhi_Zee,     "genmetPhi/F");        // Generator level MET phi
  outTree_Zee->Branch("u1",            &u1_Zee,            "u1/F");               // parallel component of recoil
  outTree_Zee->Branch("u2",            &u2_Zee,            "u2/F");               // perpendicular component of recoil 
  outTree_Zee->Branch("q1",            &q1_Zee,            "q1/I");               // tag lepton charge
  outTree_Zee->Branch("q2",            &q2_Zee,            "q2/I");               // probe lepton charge
  outTree_Zee->Branch("pfChIso1",      &pfChIso1_Zee,      "pfChIso1/F");         // tag lepton charged hadron isolation
  outTree_Zee->Branch("pfChIso2",      &pfChIso2_Zee,      "pfChIso2/F");         // probe lepton charged hadron isolation
  outTree_Zee->Branch("dilep", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &dilep_Zee); // dilepton 4-vector
  outTree_Zee->Branch("lep1",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep1_Zee);  // tag lepton 4-vector
  outTree_Zee->Branch("lep2",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep2_Zee);  // probe lepton 4-vector
  outTree_Zee->Branch("sc1",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sc1_Zee);   // tag supercluster 4-vector
  outTree_Zee->Branch("sc2",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sc2_Zee);   // probe supercluster 4-vector
  outTree_Zee->Branch("isVetoEle1",   &isVetoEle1,      "isVetoEle1/I");   // tag lepton veto electron ID
  outTree_Zee->Branch("isLooseEle1",  &isLooseEle1,     "isLooseEle1/I");  // tag lepton loose electron ID
  outTree_Zee->Branch("isMediumEle1", &isMediumEle1,    "isMediumEle1/I"); // tag lepton medium electron ID
  outTree_Zee->Branch("isTightEle1",  &isTightEle1,     "isTightEle1/I");  // tag lepton tight electron ID
  outTree_Zee->Branch("isVetoEle2",   &isVetoEle2,      "isVetoEle2/I");   // probe lepton veto electron ID
  outTree_Zee->Branch("isLooseEle2",  &isLooseEle2,     "isLooseEle2/I");  // probe lepton loose electron ID
  outTree_Zee->Branch("isMediumEle2", &isMediumEle2,    "isMediumEle2/I"); // probe lepton medium electron ID
  outTree_Zee->Branch("isTightEle2",  &isTightEle2,     "isTightEle2/I");  // probe lepton tight electron ID
}

// ------------ method called once each job just after ending the event loop  ------------
void 
selectZee::endJob() 
{
   // The only information in the last entry of the output tree
   // is the total number of MC events processed
   nEvents_Zee = eventCounter_Zee;
   // Set everything else to zero
   matchGen_Zee=0; npv_Zee=0;
   genVPdgID_Zee=0; genVPt_Zee=0; genVPhi_Zee=0; genVy_Zee=0; genVMass_Zee=0;
   rawpfmet_Zee=0; rawpfmetPhi_Zee=0;
   type1pfmet_Zee=0; type1pfmetPhi_Zee=0;
   genmet_Zee=0; genmetPhi_Zee=0;
   u1_Zee=0; u2_Zee=0;
   q1_Zee=0; q2_Zee=0;
   pfChIso1_Zee=0; pfChIso2_Zee=0;
   dilep_Zee=0; lep1_Zee=0; lep2_Zee=0;
   sc1_Zee=0; sc2_Zee=0;
   isVetoEle1=0; isLooseEle1=0; isMediumEle1=0; isTightEle1=0;
   isVetoEle2=0; isLooseEle2=0; isMediumEle2=0; isTightEle2=0;

   outTree_Zee->Fill();

   // Save tree in output file
   outFile_Zee->Write();
   outFile_Zee->Close();

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================

  std::cout << "*" << std::endl;
  std::cout << "* SUMMARY" << std::endl;
  std::cout << "*--------------------------------------------------" << std::endl;
  std::cout << "Z -> e e" << std::endl;
  std::cout << nsel_Zee << " +/- " << sqrt(nsel_Zee) << " per 1/fb" << std::endl;
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
