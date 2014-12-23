CSA14
=====

This repository contains a setup for analyzing the miniAOD samples used in the CSA14 computing exercise in the spring and summer of 2014. The code was written in the interest of preparing for the W/Z inclusive cross section measurement at 13 TeV.

Here are the main components, organized by folder, and some instructions on how to use them:

1.plugins folder

The selection code used to make the flat ntuples are here. The cuts (pT, eta) applied to leptons are slightly looser than the final cuts used for signal extraction.

2.python folder

This folder contains the configuration files needed to run over the samples and create the ntuples. The general syntax is

cmsenv
cmsRun **configuration file** [ inputFiles_load=**text file with list of sample files** ]

To run over just one file at a time, enter the file name in the configuration file and don't use the "inputFiles_load" option.

The number of events per sample can be changed in the configuration files in the line that looks like

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(**nEvents**) )

3.plotting folder

The macros in this folder apply the final kinematic cuts to leptons and make plots of both MET and lepton characteristics. The MET plots compare 4 different types of MET:

1. Generated MET
2. Type-1 corrected PF MET - the default in CSA14 miniAOD samples
3. Type-1 + xy-shift corrected PF MET - the xy-shift correction mitigates the modulation of the phi component of the MET vector
4. Raw PF MET - the negative sum of the ET vectors of all PF candidates in an event

Also included is a macro to study the resolution of the hadronic recoil in Z->ee and Z->mm samples.

All plots can be generated through the script makePlots.sh.
