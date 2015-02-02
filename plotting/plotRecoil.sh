# This script makes plots to study the hadronic recoil resolution in Z->ee and Z->mm samples

SAMPLE_DIR=/afs/cern.ch/user/c/cmedlock/CMSSW_7_2_0/src/Test/MiniAnalyzer
ZEE_SAMPLE=DYJetsToLL_M-50_13TeV-madgraph-pythia8_PU20bx25_PHYS14_V1-v1_00000_Zee_select.root
ZMM_SAMPLE=DYJetsToLL_M-50_13TeV-madgraph-pythia8_PU20bx25_PHYS14_V1-v1_00000_Zmm_select.root

root -l plotZll_recoil.C+\(\"$SAMPLE_DIR/$ZEE_SAMPLE\"\)

#root -l plotZll_recoil.C+\(\"$SAMPLE_DIR/$ZMM_SAMPLE\"\)
