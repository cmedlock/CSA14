# This script makes plots of
#	1) Different types of MET
#	2) Basic lepton characteristics (pT, eta, phi)
#	3) Hadronic recoil resolution for Z->ee and Z->mm samples

SAMPLE_DIR=/scratch/cmedlock/PHYS14
ZLL_SAMPLES=DYJetsToLL_M-50_13TeV-madgraph-pythia8/PU20bx25_PHYS14_V1-v1/00000
WJETS_SAMPLES=WJetsToLNu_13TeV-madgraph-pythia8-tauola/PU20bx25_PHYS14_25_V1-v1/00000

#root -l -q plotWe.C+\(\"$SAMPLE_DIR/$WJETS_SAMPLES/selectWe.root\"\)
root -l -q plotWm.C+\(\"$SAMPLE_DIR/$WJETS_SAMPLES/selectWm.root\"\)

#root -l -q plotZll_recoil.C+\(\"$SAMPLE_DIR/$ZLL_SAMPLES/selectZee.root\"\)
#root -l -q plotZll_recoil.C+\(\"$SAMPLE_DIR/$ZLL_SAMPLES/selectZmm.root\"\)

rm *~ *.d *.so
