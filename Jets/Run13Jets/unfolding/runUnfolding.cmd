#!/bin/csh
cd /direct/phenix+u/arbint/Jets/Analysis/makePlots/plotPostQM15/R3/unfolding

root.exe -l <<EOF
gSystem->Load("/direct/phenix+u/arbint/Jets/Analysis/makePlots/plotJetFinal/3_Unfolding/RooUnfold/libRooUnfold.so") 
.include /direct/phenix+u/arbint/Jets/Analysis/makePlots/plotJetFinal/3_Unfolding/RooUnfold/src
.x doUnfoldingR3.C+
EOF

purgeAll
