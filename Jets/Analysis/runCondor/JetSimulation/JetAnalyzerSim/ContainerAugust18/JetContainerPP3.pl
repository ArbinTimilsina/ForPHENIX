#!/usr/local/bin/perl

use strict;

my $segments = 33000;
print "\nRunning over $segments segments.\n\n";
sleep(1);

my $nevents = 4500;
print "\nNumber of events in each segment is: $nevents\n\n";
sleep(1);

my $R = 0.2;
print "R-parameter is $R\n\n";
sleep(1);

my $minPt = 6.59;
print "Min pt is: $minPt (GeV/c)\n\n";
sleep(1);

my $ncPythia = 0.0;
my $minCfPythia = 0.0;
my $maxCfPythia = 1.0;
print "For truth: ncPythia>=$ncPythia, minCfPythia=$minCfPythia, maxCfPythia=$maxCfPythia\n\n";
sleep(1);

my $runDir = "/direct/phenix+hhj/arbint/JetSimulation/JetAnalyzerSim/Container/PP";

print "Writing out events to dir: $runDir\n\n";
sleep(2);

my $ii = 0;
for($ii=32163;$ii<$segments;$ii++) {
    my $dir = "JetContainerAnalyzer_$ii";
    print "\n************************************\n";
    print "Making directory $dir\n\n";
    system("cd $runDir; rm -rf $dir; mkdir -p $dir");

#Create glogon.kumac
    print "Creating glogon.kumac\n";
    open (GLOGON_KUMAC, ">$runDir/$dir/glogon.kumac");
    print GLOGON_KUMAC "phpythia 1 phpythia.root\n";
    print GLOGON_KUMAC "ptrig $nevents\n";
    print GLOGON_KUMAC "exit\n";
    close (GLOGON_KUMAC);

#Create JetContainerAnalyzer.cmd script
    print "Creating JetContainerAnalyzer.cmd\n";
    open (JETCONTAINERANALYZER_CMD, ">$runDir/$dir/JetContainerAnalyzer.cmd");
    print JETCONTAINERANALYZER_CMD "#!/bin/csh\n";
    print JETCONTAINERANALYZER_CMD "source /etc/csh.login\n";
    print JETCONTAINERANALYZER_CMD "foreach i (/etc/profile.d/\*.csh)\n";
    print JETCONTAINERANALYZER_CMD "     source \$i\n";
    print JETCONTAINERANALYZER_CMD "end\n";
    print JETCONTAINERANALYZER_CMD "source /opt/phenix/bin/phenix_setup.csh -n ana.475\n";
    print JETCONTAINERANALYZER_CMD "source /opt/phenix/bin/odbcini_setup.csh\n";
    print JETCONTAINERANALYZER_CMD "setenv LD_LIBRARY_PATH /direct/phenix\+u/arbint/Jets/Analysis/install/lib:\$\{LD_LIBRARY_PATH\}\n";
    print JETCONTAINERANALYZER_CMD "if (\$\?\_CONDOR_SCRATCH_DIR\) then\n";
    print JETCONTAINERANALYZER_CMD "set workDir = \$\{_CONDOR_SCRATCH_DIR\}\n";
    print JETCONTAINERANALYZER_CMD "endif\n";
    print JETCONTAINERANALYZER_CMD "cd \$workDir\n";
    print JETCONTAINERANALYZER_CMD "pwd\n";

#Set PYTHIA
    print "Setting PYTHIA\n\n";
    print JETCONTAINERANALYZER_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pythia.cfg .\n";
    print JETCONTAINERANALYZER_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/generateVertexForPP.C .\n";
    print JETCONTAINERANALYZER_CMD "root -l -q -b \"generateVertexForPP.C($nevents)\"\n";
    print JETCONTAINERANALYZER_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/generatePythia.C .\n";
    print JETCONTAINERANALYZER_CMD "root -l -q -b \"generatePythia.C($R, $ncPythia, $minPt, $minCfPythia, $maxCfPythia, $nevents)\"\n";

#Set PISA
    print "Setting PISA\n";
    print JETCONTAINERANALYZER_CMD "ln -sf $runDir/$dir/glogon.kumac .\n";
    print JETCONTAINERANALYZER_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pisaLinker/* .\n";
    print JETCONTAINERANALYZER_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pisaToDSTLinker/* .\n";
    print JETCONTAINERANALYZER_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pcDeadArea/pad_deadroc_360934.dat pad_deadroc.dat\n";
    print JETCONTAINERANALYZER_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pcDeadArea/pad_deadch_360934.dat pad_deadch.dat\n";
    print JETCONTAINERANALYZER_CMD "pisa < pisa.input\n";

#Set PISA to DST
    print "Setting PISA to DST\n\n";
    print JETCONTAINERANALYZER_CMD "root -l -q -b \"pisaToDST.C($nevents, 360934)\" \n";

#Set the code
    print "Setting the code\n\n";
    print JETCONTAINERANALYZER_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/runJetContainerAnalyzerPP.C .\n";
    print JETCONTAINERANALYZER_CMD "root -l -q -b \"runJetContainerAnalyzerPP.C($R, $minPt, $ncPythia, $minCfPythia, $maxCfPythia, $nevents)\"\n";

    close(JETCONTAINERANALYZER_CMD);

#Submit the job
    system("cd $runDir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/runCondor/JetSimulation/JetAnalyzerSim/Container/JetContainerAnalyzer.job");
    system("cd $runDir/$dir; chmod a+x JetContainerAnalyzer.cmd");
    system("cd $runDir/$dir; condor_submit JetContainerAnalyzer.job");

#Sleep for 35 seconds
    sleep(5);
}
