# Lxplus Batch Job Script
CMSSW_PROJECT_SRC="CMSSW_7_0_6_patch1/src"
CFG_FILE="Test/MiniAnalyzer/python/selectWm_cfg.py"
OUTPUT_FILE="Analyzer_Output.root"
TOP="$PWD"

cd /afs/cern.ch/user/c/cmedlock/$CMSSW_PROJECT_SRC
echo "$PWD"
eval `scramv1 runtime -sh`
cd $TOP
cmsRun /afs/cern.ch/user/c/cmedlock/$CMSSW_PROJECT_SRC/$CFG_FILE
rfcp Analyzer_Output.root /afs/cern.ch/user/c/cmedlock/public/$OUTPUT_FILE

