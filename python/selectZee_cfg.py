import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

options = VarParsing.VarParsing("analysis")
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

# load the input files
process.source = cms.Source("PoolSource",
    # for a single input file
#    fileNames = cms.untracked.vstring(
#        '/store/user/arapyan/WplusToENu_CT10_13TeV-powheg-pythia8/Spring14dr-PU-S14-POSTLS170_V6-v1/140723_170241/0000/miniAOD-prod_PAT_19.root'
    fileNames = cms.untracked.vstring(options.inputFiles
    )
)

# load the selection plugin
process.selectZee = cms.EDAnalyzer('selectZee',
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    electrons = cms.InputTag("slimmedElectrons"),
    mets = cms.InputTag("slimmedMETs"),
    pfCands = cms.InputTag("packedPFCandidates"),
)

# define path of execution
process.p = cms.Path(process.selectZee)
