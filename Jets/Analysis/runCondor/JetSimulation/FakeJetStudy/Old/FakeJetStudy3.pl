#!/usr/local/bin/perl

use strict;

my $runDir = "/direct/phenix+hhj/arbint/JetSimulation/FakeJetStudy/FakeJetStudy3";

my $segments = 500;

my $ii = 0;
my $nevents = 500;
for($ii=0;$ii<$segments;$ii++) {
    my $dir = "FakeJetStudy_$ii";
    print "\n***********************************\n";
    print "Making directory $dir\n";
    system("cd $runDir; rm -rf $dir; mkdir -p $dir");

#First run sHIJING
    print "Making the xml file for sHijing\n";

    my $SEED= $ii + 333333;
    print "Seed set to ${SEED}\n";

    system("cd $runDir/$dir; echo \"<?xml version=\"1.0\"?>\" > sHijing.xml");
    system("cd $runDir/$dir; echo \"<HIJING>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \" <EFRM>200.0</EFRM>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \" <FRAME>CMS</FRAME>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \" <PROJ>A</PROJ>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \" <TARG>A</TARG>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \" <IAP>197</IAP>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \" <IZP>79</IZP>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \" <IAT>64</IAT>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \" <IZT>29</IZT>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \" <N>${nevents}</N>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \" <BMIN>0.0</BMIN>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \" <BMAX>20.0</BMAX>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \" <KEEP_SPECTATORS>1</KEEP_SPECTATORS>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \" <OUTPUT>sHijing.dat</OUTPUT>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \" <RANDOM>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \"  <SEED>${SEED}</SEED>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \" </RANDOM>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \" <IHPR2>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \"  <4>0</4>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \"  <12>0</12>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \"  <21>1</21>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \" </IHPR2>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \" <FastJet>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \"  <Algorithm>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \"   <Name>AntikT</Name>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \"   <R>0.2</R>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \"   <PID>2000000</PID>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \"  </Algorithm>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \" </FastJet>\" >> sHijing.xml");
    system("cd $runDir/$dir; echo \"</HIJING>\" >> sHijing.xml");

#Create FakeJetStudy.cmd script and run sHIJING
    system("cd $runDir/$dir; echo \'#!/bin/csh\' >> FakeJetStudy.cmd");
    system("cd $runDir/$dir; echo \'cd $runDir/$dir\' >> FakeJetStudy.cmd");
    system("cd $runDir/$dir; echo \'sHijing' >> FakeJetStudy.cmd");

#Convert HepMC output to PHPYTHIA output
    print "Setting up HepMCPyNodeReader \n";
    system("cd $runDir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/runHepMCPyNodeReader.C");
    system("cd $runDir/$dir; echo \'root -l -q -b runHepMCPyNodeReader.C\' >> FakeJetStudy.cmd");

#Create glogon.kumac
    print "Creating glogon.kumac\n";
    system("cd $runDir/$dir; echo \'phpythia 10 phpythia.root\nptrig $nevents\nexit\' > glogon.kumac");

#Set PISA
    print "Setting PISA\n";
    system("cd $runDir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pisaLinker/* .");
    system("cd $runDir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pisaToDSTLinker/* .");
    system("cd $runDir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pcDeadArea/pad_deadroc_372524.dat pad_deadroc.dat");
    system("cd $runDir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pcDeadArea/pad_deadch_372524.dat pad_deadch.dat");
    system("cd $runDir/$dir; echo \'pisa <pisa.input\' >> FakeJetStudy.cmd");

#Set PISA to DST
    print "Setting PISA to DST\n";
    system("cd $runDir/$dir; echo \'root -l -q -b \"pisaToDST.C($nevents, 372524)\"\' >> FakeJetStudy.cmd");

#Set the code
    print "Setting the code\n";
    system("cd $runDir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/runFakeJetStudy.C");
    system("cd $runDir/$dir; echo \'root -l -q -b runFakeJetStudy.C\' >> FakeJetStudy.cmd");

#Want to delete all files
    system("cd $runDir/$dir; echo \'find . -type f -not -name 'FakeJetStudy.root' -delete\' >> FakeJetStudy.cmd");
    system("cd $runDir/$dir; echo \'find . -type l -delete\' >> FakeJetStudy.cmd");

#Submit the job
    system("cd $runDir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/runCondor/JetSimulation/FakeJetStudy/FakeJetStudy.job");
    system("cd $runDir/$dir; chmod a+x FakeJetStudy.cmd");
    system("cd $runDir/$dir; condor_submit FakeJetStudy.job");

    print "Done!!\n";

#Sleep for 15 seconds
    sleep(15);
}


