#!/usr/local/bin/perl                                                                                                                             
use strict;

my $segments = 25000;
print "\nRunning over $segments segments.\n\n";
sleep(1);

my $nevents = 25000;
print "\nNumber of events in each segment is: $nevents\n\n";
sleep(1);

my $eventDir = "/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/arbint/Run13Simulation/JetHadCorrection/200GeV/parton";

print "Writing out events to dir: $eventDir\n\n";
sleep(1);

my $ii = 0;
for($ii=0;$ii<$segments;$ii++) {
    my $dir = "JetHadCorrection_$ii";
    print "\n************************************\n";
    print "Making directory $dir\n\n";
    system("cd $eventDir; rm -rf $dir; mkdir -p $dir");

#Set PYTHIA
    print "Setting PYTHIA\n\n";

    open (JETHADCORRECTION_CMD, ">$eventDir/$dir/JetHadCorrection.cmd");
    print JETHADCORRECTION_CMD "#!/bin/csh\n";
    print JETHADCORRECTION_CMD "source /etc/csh.login\n";
    print JETHADCORRECTION_CMD "foreach i (/etc/profile.d/\*.csh)\n";
    print JETHADCORRECTION_CMD "     source \$i\n";
    print JETHADCORRECTION_CMD "end\n";
    print JETHADCORRECTION_CMD "source /opt/phenix/bin/phenix_setup.csh ana\n";
    print JETHADCORRECTION_CMD "source /opt/phenix/bin/odbcini_setup.csh\n";
    print JETHADCORRECTION_CMD "setenv LD_LIBRARY_PATH /direct/phenix\+u/arbint/Jets/Analysis/install/lib:\$\{LD_LIBRARY_PATH\}\n";
    print JETHADCORRECTION_CMD "if (\$\?\_CONDOR_SCRATCH_DIR\) then\n";
    print JETHADCORRECTION_CMD "set workDir = \$\{_CONDOR_SCRATCH_DIR\}\n";
    print JETHADCORRECTION_CMD "endif\n";
    print JETHADCORRECTION_CMD "cd \$workDir\n";
    print JETHADCORRECTION_CMD "pwd\n";

    print JETHADCORRECTION_CMD "ln -sf /direct/phenix+u/arbint/Jets/Run13Jets/condor/JetHadCorrection/200GeV/pythia.parton pythia.cfg\n";
    print JETHADCORRECTION_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/generatePythiaCommon.C .\n";
    print JETHADCORRECTION_CMD "root -l -q -b \"generatePythiaCommon.C($nevents)\"\n";

#Set the code
    print "Setting the code\n\n";
    print JETHADCORRECTION_CMD "ln -sf /direct/phenix+u/arbint/Jets/Run13Jets/macros/runJetHadCorrection.C .\n";
    print JETHADCORRECTION_CMD "root -l -q -b \"runJetHadCorrection.C()\"\n";

    close(JETHADCORRECTION_CMD);

#Submit the job
    system("cd $eventDir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Run13Jets/condor/JetHadCorrection//JetHadCorrection.job");
    system("cd $eventDir/$dir; chmod a+x JetHadCorrection.cmd");
    system("cd $eventDir/$dir; condor_submit JetHadCorrection.job");

#Sleep for 1 seconds
    sleep(1);
}
