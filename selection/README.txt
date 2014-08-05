This folder contains the code used to:

1.Implement the xy-shift corrections to the type-1 corrected PF MET (the default MET for CSA14 miniAOD samples)
2.Make plots of interesting variables in order to verify that the samples are usable

The relevant files are:

plotWe.C
plotWm.C

To run:

root -l -q plotWe.C+ --> The output should be *.png files of electron pT, eta, and phi and MET pT and phi in your working directory.
root -l -q plotWm.C+ --> The output should be *.png files of muon pT, eta, and phi and MET pT and phi in your working directory.

