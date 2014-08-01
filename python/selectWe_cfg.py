# To run
# >> cmsenv
# >> cmsRun selectWm_cfg.py inputFiles_load=WplusToENu_CT10_13TeV-powheg-pythia8.txt
# The output should be a file called
# WplusToENu_CT10_13TeV-powheg-pythia8_SELECT.root
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
#        '/store/user/arapyan/WplusToENu_CT10_13TeV-powheg-pythia8/Spring14dr-PU-S14-POSTLS170_V6-v1/140723_170241/0000/miniAOD-prod_PAT_19.root'
    fileNames = cms.untracked.vstring(options.inputFiles
    )
)

process.demo = cms.EDAnalyzer('selectWe',
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    electrons = cms.InputTag("slimmedElectrons"),
    mets = cms.InputTag("slimmedMETs"),
    pfCands = cms.InputTag("packedPFCandidates"),
)

process.p = cms.Path(process.demo)
