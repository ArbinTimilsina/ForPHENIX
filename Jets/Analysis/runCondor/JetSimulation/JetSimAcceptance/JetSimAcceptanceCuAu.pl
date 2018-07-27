#!/usr/local/bin/perl                                                                                                                              
use strict;

my $segments = 40000;
print "\nRunning over $segments segments.\n\n";
sleep(1);

my $nevents = 100;
print "\nNumber of events in each segment is: $nevents\n\n";
sleep(1);

my $evntdir = "/direct/phenix+hhj/arbint/JetSimulation/JetSimAcceptance/CuAu";
print "Writing out events to dir: $evntdir\n\n";
sleep(5);

my $ii = 0;
for($ii=0;$ii<$segments;$ii++) {
    my $dir = "JetSimAcceptance_$ii";
    print "\n************************************\n";
    print "Making directory $dir\n\n";
    system("cd $evntdir; rm -rf $dir; mkdir -p $dir");

#Set HIJING
    print "Setting HIJING\n\n";

    system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/phhijing.cfg");
    system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/generateHijing.C");
    system("cd $evntdir/$dir; echo \'#!/bin/csh\ncd $evntdir/$dir\' >> JetSimAcceptance.cmd");
    system("cd $evntdir/$dir; echo \'root -l -q -b \"generateHijing.C($nevents)\"\' >> JetSimAcceptance.cmd");

    print "Creating glogon.kumac\n\n";
#Create glogon.kumac 
    system("cd $evntdir/$dir; echo \'phpythia 1 phhijing.root\nptrig $nevents\nexit\' > glogon.kumac");

#Set PISA
    print "Setting PISA\n";
    system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pisaLinker/* .");
    system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pisaToDSTLinker/* .");
    system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pcDeadArea/pad_deadroc_372524.dat pad_deadroc.dat");
    system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pcDeadArea/pad_deadch_372524.dat pad_deadch.dat");
    system("cd $evntdir/$dir; echo \'pisa <pisa.input\' >> JetSimAcceptance.cmd");

#Set PISA to DST
    print "Setting PISA to DST\n\n";
    system("cd $evntdir/$dir; echo \'root -l -q -b \"pisaToDST.C($nevents, 372524)\"\' >> JetSimAcceptance.cmd");

#Set the code
    print "Setting the code\n\n";
    system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/runJetSimAcceptance.C");
    system("cd $evntdir/$dir; echo \'root -l -q -b \"runJetSimAcceptance.C($nevents, true)\"\' >> JetSimAcceptance.cmd");

#Want to delete large files
    system("cd $evntdir/$dir; echo \'rm -rf dummyFile.root phhijing.root PISAEvent.root SimDST.root\' >> JetSimAcceptance.cmd");
    system("cd $evntdir/$dir; echo \'rm -rf phpy_xsec.root gintphnx.hbk \' >> JetSimAcceptance.cmd");

#Submit the job
    system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/runCondor/JetSimulation/JetSimAcceptance/JetSimAcceptance.job");
    system("cd $evntdir/$dir; chmod a+x JetSimAcceptance.cmd");
    system("cd $evntdir/$dir; condor_submit JetSimAcceptance.job");

#Sleep for 8 seconds
    sleep(8);
}
