SAMPLE_DIR=/scratch3/cmedlock/PHYS14/Selection

#
# Z->ee and Z->mm dilepton mass plots: dependency on PU
#
root -l plotDileptonMass.C+\(\"$SAMPLE_DIR/Zee/zee_select_PU4bx50.root\",\"$SAMPLE_DIR/Zee/zee_select_PU20bx25.root\"\)
root -l plotDileptonMass.C+\(\"$SAMPLE_DIR/Zmm/zmm_select_PU4bx50.root\",\"$SAMPLE_DIR/Zmm/zmm_select_PU20bx25.root\"\)

#
# Z->ee MET resolution and hadronic recoil plots
#
#root -l -q plotMET_res_recoil.C+\(\"$SAMPLE_DIR/Zee/zee_select_PU4bx50.root\",\"./zee_metplots_PU4bx50.root\"\)
#root -l -q plotMET_res_recoil.C+\(\"$SAMPLE_DIR/Zee/zee_select_PU20bx25.root\",\"./zee_metplots_PU20bx25.root\"\)
#root -l mergeMetPlots_res_recoil.C+\(\"./zee_metplots_PU4bx50.root\",\"./zee_metplots_PU20bx25.root\"\)

#
# Z->mm MET resolution and hadronic recoil plots
#
#root -l -q plotMET_res_recoil.C+\(\"$SAMPLE_DIR/Zmm/zmm_select_PU4bx50.root\",\"./zmm_metplots_PU4bx50.root\"\)
#root -l -q plotMET_res_recoil.C+\(\"$SAMPLE_DIR/Zmm/zmm_select_PU20bx25.root\",\"./zmm_metplots_PU20bx25.root\"\)
#root -l mergeMetPlots_res_recoil.C+\(\"./zmm_metplots_PU4bx50.root\",\"./zmm_metplots_PU20bx25.root\"\)

#
# W->ev MET distributions: dependency on PU
#
#root -l -q plotMET.C+\(\"$SAMPLE_DIR/Wenu/we_select_PU4bx50.root\",\"./wenu_metplots_PU4bx50.root\"\)
#root -l -q plotMET.C+\(\"$SAMPLE_DIR/Wenu/we_select_PU20bx25.root\",\"./wenu_metplots_PU20bx25.root\"\)
#root -l mergeMetPlots.C+\(\"./wenu_metplots_PU4bx50.root\",\"./wenu_metplots_PU20bx25.root\"\)

#
# W->mv MET distributions: dependency on PU
#
#root -l -q plotMET.C+\(\"$SAMPLE_DIR/Wmunu/wm_select_PU4bx50.root\",\"./wmunu_metplots_PU4bx50.root\"\)
#root -l -q plotMET.C+\(\"$SAMPLE_DIR/Wmunu/wm_select_PU20bx25.root\",\"./wmunu_metplots_PU20bx25.root\"\)
#root -l mergeMetPlots.C+\(\"./wmunu_metplots_PU4bx50.root\",\"./wmunu_metplots_PU20bx25.root\"\)

rm *~ *.d *.so
