CSA14
=====

This repository contains a setup for analyzing the miniAOD samples used in the CSA14 computing exercise in the spring and summer of 2014. The code was written in the interest of preparing for the W/Z inclusive cross section measurement at 13 TeV. Almost all of it is just an adaptation of Kevin Sung's original code used for the 8 TeV analysis, which can be found here:

https://github.com/ksung25/UserCode/tree/master/EWKAna 

Here are the main components, organized by folder, and some instructions on how to use them:

1.plugins folder

The selection code used to make the flat ntuples are here. The cuts (pT, eta) applied to leptons are slightly looser than the final cuts used for signal extraction.

Also included is code to make ntuples used for tag-and-probe efficiency studies. These ntuples differ from the ones used in signal extraction in that they save information about two leptons rather than one, as well as information about their associated superclusters. They are meant to be inputs for the plotEff.C macro in the Efficiency/TagAndProbe folder.

2.python folder

This folder contains the configuration files needed to run over the samples and create the ntuples. The general syntax is

cmsenv
cmsRun **configuration file** [ inputFiles_load=**text file with list of sample files** ]

To run over just one file at a time, enter the file name in the configuration file and don't use the "inputFiles_load" option.

The number of events per sample can be changed in the configuration files in the line that looks like

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(**nEvents**) )

3.plotting folder

The macros in this folder make plots relevant for both W and Z signal extraction, and compare different pile-up scenarios:

1.Dilepton mass in the Z signal samples
2.Type-1 corrected PF MET (the default in CSA14 miniAOD samples) in the W signal samples
3.Type-1 corrected PF MET resolution in the W signal samples
4.Hadronic recoil parallel and perpendicular component resolution

All plots can be generated through the script makePlots.sh.
