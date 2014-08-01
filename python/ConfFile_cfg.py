import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:/afs/cern.ch/user/c/cmedlock/CMSSW_7_0_6_patch1/src/Test/MiniAnalyzer/TT_Tune4C_13TeV-pythia8-tauola_PU_S14_PAT.root'
        'file:/afs/cern.ch/work/a/arapyan/public/forCatherine/miniAOD-prod_PAT.root'
    )
)

process.demo = cms.EDAnalyzer('MiniAnalyzer',
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("slimmedJets"),
    fatjets = cms.InputTag("slimmedJetsAKS"),
    mets = cms.InputTag("slimmedMETs"),
)


process.p = cms.Path(process.demo)
