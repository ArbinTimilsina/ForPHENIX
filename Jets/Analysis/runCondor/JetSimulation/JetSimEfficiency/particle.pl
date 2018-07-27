#!/usr/local/bin/perl                                                                                                                             
use strict;

my $segments = 20000;
print "\nRunning over $segments segments.\n\n";
sleep(1);

my $nevents = 500;
print "\nNumber of events in each segment is: $nevents\n\n";
sleep(1);

my $R = 0.2;
print "R-parameter is $R\n\n";
sleep(1);

my $minPt = 7.9;
print "Min pt is: $minPt (GeV/c)\n\n";
sleep(1);

my $nc = 3.0;
my $minCf = 0.2;
my $maxCf = 0.7;
print "For regular: nc>=$nc, minCf=$minCf, maxCf=$maxCf\n\n";
sleep(1);

my $ncPythia = 3.0;
my $minCfPythia = 0.2;
my $maxCfPythia = 0.7;
print "For truth: ncPythia>=$ncPythia, minCfPythia=$minCfPythia, maxCfPythia=$maxCfPythia\n\n";
sleep(1);

my $evntdir = "/direct/phenix+scratch/arbint/JetSimEfficiency/particle";
print "Writing out events to dir: $evntdir\n\n";
sleep(3);


my $ii = 0;
for($ii=0;$ii<$segments;$ii++) {
    my $dir = "JetAnalyzerSim_$ii";
    print "\n************************************\n";
    print "Making directory $dir\n\n";
    system("cd $evntdir; rm -rf $dir; mkdir -p $dir");

#Set PYTHIA
    print "Setting PYTHIA\n\n";
    system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pythia.cfg");
    system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/generateVertexForPP.C");
    system("cd $evntdir/$dir; echo \'#!/bin/csh\' >> JetAnalyzerSim.cmd");
    system("cd $evntdir/$dir; echo \'source /opt/phenix/bin/odbcini_setup.csh\' >> JetAnalyzerSim.cmd");
    system("cd $evntdir/$dir; echo \'cd $evntdir/$dir\' >> JetAnalyzerSim.cmd");
    system("cd $evntdir/$dir; echo \'root -l -q -b \"generateVertexForPP.C($nevents)\"\' >> JetAnalyzerSim.cmd");
    system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/generatePythia.C");
    system("cd $evntdir/$dir; echo \'root -l -q -b \"generatePythia.C($R, $ncPythia, $minPt, $minCfPythia, $maxCfPythia, $nevents, true)\"\' >> JetAnalyzerSim.cmd");

    print "Creating glogon.kumac\n\n";
#Create glogon.kumac 
    system("cd $evntdir/$dir; echo \'phpythia 1 phpythia.root\nptrig $nevents\nexit\' > glogon.kumac");

#Set PISA
    print "Setting PISA\n";
    system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pisaLinker/* .");
    system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pisaToDSTLinker/* .");
    system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pcDeadArea/pad_deadroc_360934.dat pad_deadroc.dat");
    system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pcDeadArea/pad_deadch_360934.dat pad_deadch.dat");
    system("cd $evntdir/$dir; echo \'pisa <pisa.input\' >> JetAnalyzerSim.cmd");

#Set PISA to DST
    print "Setting PISA to DST\n\n";
    system("cd $evntdir/$dir; echo \'root -l -q -b \"pisaToDST.C($nevents, 360934)\"\' >> JetAnalyzerSim.cmd");

#Set the code
    print "Setting the code\n\n";
    system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/runJetAnalyzerSimPP.C");
    system("cd $evntdir/$dir; echo \'root -l -q -b \"runJetAnalyzerSimPP.C($R, $minPt, $nc, $minCf, $maxCf, $ncPythia, $minCfPythia, $maxCfPythia, 0.2,0.2,0.2, $nevents, false, true)\"\' >> JetAnalyzerSim.cmd");

#Want to delete large files
    system("cd $evntdir/$dir; echo \'rm -rf dummyFile.root phpythia.root PISAEvent.root SimDST.root\' >> JetAnalyzerSim.cmd");
    system("cd $evntdir/$dir; echo \'rm -rf phpy_xsec.root gintphnx.hbk \' >> JetAnalyzerSim.cmd");

#Submit the job
    system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/runCondor/JetSimulation/JetSimEfficiency/JetAnalyzerSim.job");
    system("cd $evntdir/$dir; chmod a+x JetAnalyzerSim.cmd");
    system("cd $evntdir/$dir; condor_submit JetAnalyzerSim.job");

#Sleep for 2 seconds
    sleep(2);
}
