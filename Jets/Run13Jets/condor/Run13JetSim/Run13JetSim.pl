#!/usr/local/bin/perl

use strict;

my $segments = 10000;
print "\nRunning over $segments segments.\n\n";
sleep(1);

my $nevents = 4500;
print "\nNumber of events in each segment is: $nevents\n\n";
sleep(1);

my $minPt = 5.0;
print "Min pt is: $minPt (GeV/c)\n\n";
sleep(1);

my $runDir = "/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/arbint/Run13Simulation/Run13JetSim";

print "Writing out events to dir: $runDir\n\n";
sleep(2);

my $ii = 0;
for($ii=0;$ii<$segments;$ii++) {
    my $dir = "Run13JetSim_$ii";
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

#Create Run13JetSim.cmd script
    print "Creating Run13JetSim.cmd\n";
    open (RUN13JETSIM_CMD, ">$runDir/$dir/Run13JetSim.cmd");
    print RUN13JETSIM_CMD "#!/bin/csh\n";
    print RUN13JETSIM_CMD "source /etc/csh.login\n";
    print RUN13JETSIM_CMD "foreach i (/etc/profile.d/\*.csh)\n";
    print RUN13JETSIM_CMD "     source \$i\n";
    print RUN13JETSIM_CMD "end\n";
    print RUN13JETSIM_CMD "source /opt/phenix/bin/phenix_setup.csh -n ana.516\n";
    print RUN13JETSIM_CMD "source /opt/phenix/bin/odbcini_setup.csh\n";
    print RUN13JETSIM_CMD "setenv LD_LIBRARY_PATH /direct/phenix+u/arbint/Jets/Run13Jets/install/lib:\$\{LD_LIBRARY_PATH\}\n";
    print RUN13JETSIM_CMD "setenv TSEARCHPATH /direct/phenix+u/arbint/Jets/Run13Jets/install\n";
    print RUN13JETSIM_CMD "if (\$\?\_CONDOR_SCRATCH_DIR\) then\n";
    print RUN13JETSIM_CMD "set workDir = \$\{_CONDOR_SCRATCH_DIR\}\n";
    print RUN13JETSIM_CMD "endif\n";
    print RUN13JETSIM_CMD "cd \$workDir\n";
    print RUN13JETSIM_CMD "pwd\n";

#Set PYTHIA
    print "Setting PYTHIA\n\n";
    print RUN13JETSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Run13Jets/macros/pythia.cfg .\n";
    print RUN13JETSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Run13Jets/macros/generateVertex.C .\n";
    print RUN13JETSIM_CMD "root -l -q -b \"generateVertex.C($nevents)\"\n";
    print RUN13JETSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Run13Jets/macros/generatePythia.C .\n";
    print RUN13JETSIM_CMD "root -l -q -b \"generatePythia.C($minPt, $nevents)\"\n";

#Set PISA
    print "Setting PISA\n";
    print RUN13JETSIM_CMD "ln -sf $runDir/$dir/glogon.kumac .\n";
    print RUN13JETSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Run13Jets/macros/pisaLinker/* .\n";
    print RUN13JETSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Run13Jets/macros/pisaToDSTLinker/* .\n";
    print RUN13JETSIM_CMD "pisa < pisa.input\n";

#Set PISA to DST
    print "Setting PISA to DST\n\n";
    print RUN13JETSIM_CMD "root -l -q -b \"pisaToDST.C($nevents)\" \n";

#Set the code
    print "Setting the code\n\n";
    print RUN13JETSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Run13Jets/macros/runRun13JetSim.C .\n";
    print RUN13JETSIM_CMD "root -l -q -b \"runRun13JetSim.C($minPt, $nevents)\"\n";

    close(RUN13JETSIM_CMD);

#Submit the job
    system("cd $runDir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Run13Jets/condor/Run13JetSim/Run13JetSim.job .");
    system("cd $runDir/$dir; chmod a+x Run13JetSim.cmd");
    system("cd $runDir/$dir; condor_submit Run13JetSim.job");

#Sleep for 15 seconds
    sleep(15);
}
