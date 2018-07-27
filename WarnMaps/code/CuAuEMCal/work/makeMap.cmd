#!/bin/csh

echo "\nDeleting plots \n"
cd /direct/phenix+hhj/arbint/plots/WarnMaps;
rm -rf *;

echo "\nDeleting old map and info\n"
cd /direct/phenix+u/arbint/Jets/Analysis/warnmap;
rm -rf warnmapCuAu.txt infoMapCuAu;

cd /direct/phenix+u/arbint/WarnMaps/code/CuAuEMCal/work; 

echo "\nDeleting old mean and sigma\n"
rm -rf meanSigma.txt;

echo "\nGetting mean and sigma\n"
root -l -q 'getSigma.C+'
purgeAll

echo "\nMaking map\n"
root -l -q 'makeMap.C+'
purgeAll
