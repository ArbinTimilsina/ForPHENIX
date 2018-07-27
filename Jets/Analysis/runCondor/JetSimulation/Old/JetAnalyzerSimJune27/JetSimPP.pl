#!/usr/local/bin/perl

use strict;

my $segments = 25000;
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

my $nc = 3.0;
my $minCf = 0.2;
my $maxCf = 0.7;
print "For regular: nc>=$nc, minCf=$minCf, maxCf=$maxCf\n\n";
sleep(1);

my $ncPythia = 0.0;
my $minCfPythia = 0.0;
my $maxCfPythia = 1.0;
print "For truth: ncPythia>=$ncPythia, minCfPythia=$minCfPythia, maxCfPythia=$maxCfPythia\n\n";
sleep(1);

my $runDir = "/direct/phenix+hhj/arbint/JetSimulation/JetAnalyzerSim/PP";

print "Writing out events to dir: $runDir\n\n";
sleep(10);

my $ii = 0;
for($ii=0;$ii<$segments;$ii++) {
    my $dir = "JetAnalyzerSim_$ii";
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

#Create JetAnalyzerSim.cmd script
    print "Creating JetAnalyzerSim.cmd\n";
    open (JETANALYZERSIM_CMD, ">$runDir/$dir/JetAnalyzerSim.cmd");
    print JETANALYZERSIM_CMD "#!/bin/csh\n";
    print JETANALYZERSIM_CMD "source /etc/csh.login\n";
    print JETANALYZERSIM_CMD "foreach i (/etc/profile.d/\*.csh)\n";
    print JETANALYZERSIM_CMD "     source \$i\n";
    print JETANALYZERSIM_CMD "end\n";
    print JETANALYZERSIM_CMD "source /opt/phenix/bin/phenix_setup.csh -n ana.458\n";
    print JETANALYZERSIM_CMD "source /opt/phenix/bin/odbcini_setup.csh\n";
    print JETANALYZERSIM_CMD "setenv LD_LIBRARY_PATH /direct/phenix\+u/arbint/Jets/Analysis/install/lib:\$\{LD_LIBRARY_PATH\}\n";
    print JETANALYZERSIM_CMD "if (\$\?\_CONDOR_SCRATCH_DIR\) then\n";
    print JETANALYZERSIM_CMD "set workDir = \$\{_CONDOR_SCRATCH_DIR\}\n";
    print JETANALYZERSIM_CMD "endif\n";
    print JETANALYZERSIM_CMD "cd \$workDir\n";
    print JETANALYZERSIM_CMD "pwd\n";

#Set PYTHIA
    print "Setting PYTHIA\n\n";
    print JETANALYZERSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pythia.cfg .\n";
    print JETANALYZERSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/generateVertexForPP.C .\n";
    print JETANALYZERSIM_CMD "root -l -q -b \"generateVertexForPP.C($nevents)\"\n";
    print JETANALYZERSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/generatePythia.C .\n";
    print JETANALYZERSIM_CMD "root -l -q -b \"generatePythia.C($R, $ncPythia, $minPt, $minCfPythia, $maxCfPythia, $nevents)\"\n";

#Set PISA
    print "Setting PISA\n";
    print JETANALYZERSIM_CMD "ln -sf $runDir/$dir/glogon.kumac .\n";
    print JETANALYZERSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pisaLinker/* .\n";
    print JETANALYZERSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pisaToDSTLinker/* .\n";
    print JETANALYZERSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pcDeadArea/pad_deadroc_360934.dat pad_deadroc.dat\n";
    print JETANALYZERSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pcDeadArea/pad_deadch_360934.dat pad_deadch.dat\n";
    print JETANALYZERSIM_CMD "pisa < pisa.input\n";

#Set PISA to DST
    print "Setting PISA to DST\n\n";
    print JETANALYZERSIM_CMD "root -l -q -b \"pisaToDST.C($nevents, 360934)\" \n";

#Set the code
    print "Setting the code\n\n";
    print JETANALYZERSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/runJetAnalyzerSimPP.C .\n";
    print JETANALYZERSIM_CMD "root -l -q -b \"runJetAnalyzerSimPP.C($R, $minPt, $nc, $minCf, $maxCf, $ncPythia, $minCfPythia, $maxCfPythia, 0.2, $nevents)\"\n";

    close(JETANALYZERSIM_CMD);

#Submit the job
    system("cd $runDir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/runCondor/JetSimulation/JetAnalyzerSim/JetAnalyzerSim.job");
    system("cd $runDir/$dir; chmod a+x JetAnalyzerSim.cmd");
    system("cd $runDir/$dir; condor_submit JetAnalyzerSim.job");

#Sleep for 20 seconds
    sleep(20);
}
