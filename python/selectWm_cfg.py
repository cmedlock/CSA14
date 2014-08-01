# To run
# >> cmsenv
# >> cmsRun selectWm_cfg.py inputFiles_load=WplusToMuNu_CT10_13TeV-powheg-pythia8.txt
# The output should be a file called
# WplusToMuNu_CT10_13TeV-powheg-pythia8_SELECT.root
# in your working directory.

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

options = VarParsing.VarParsing("analysis")
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )

process.source = cms.Source("PoolSource",
    # Uncomment the next 2 lines to run over a single file
    # For local files, use "  fileNames = cms.untracked.vstring('file:file_name.root')  "
#    fileNames = cms.untracked.vstring(
#        'file:/afs/cern.ch/work/a/arapyan/public/forCatherine/miniAOD-prod_PAT.root'
    fileNames = cms.untracked.vstring(options.inputFiles
    )
)

process.demo = cms.EDAnalyzer('selectWm',
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    mets = cms.InputTag("slimmedMETs"),
    pfCands = cms.InputTag("packedPFCandidates"),
)

process.p = cms.Path(process.demo)
