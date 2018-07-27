#!/bin/csh
cd /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pcDeadArea
rm -rf *.dat tempdead;

echo "\nMaking Files\n"
root -l -q 'getBadRoc.C(0)'
root -l -q 'getBadRoc.C(1)'
root -l -q 'getDeadCh.C(0)'
root -l -q 'getDeadCh.C(1)'
purgeAll

rm -rf tempdead;

