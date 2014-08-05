This folder contains the python configuration files needed to run the selection macros in the plugins/ folder.

The relevant files are:

python/selectWe_cfg.py + WplusToENu_CT10_13TeV-powheg-pythia8.txt
python/selectWm_cfg.py + WplusToMuNu_CT10_13TeV-powheg-pythia8.txt

To run:

cmsRun python/selectWe_cfg.py inputFiles_load=WplusToENu_CT10_13TeV-powheg-pythia8.txt  --> The output should be a flat ntuple called Wenu_p_select.root in your working directory.
cmsRun python/selectWm_cfg.py inputFiles_load=WplusToMuNu_CT10_13TeV-powheg-pythia8.txt --> The output should be a flat ntuple called Wmunu_p_select.root in your working directory.

