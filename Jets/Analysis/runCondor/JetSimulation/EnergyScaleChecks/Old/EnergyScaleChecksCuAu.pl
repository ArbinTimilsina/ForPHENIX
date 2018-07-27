#!/usr/local/bin/perl                                                                                                                              
use strict;

my $segments = 40000;
print "\nRunning over $segments segments.\n\n";
sleep(1);

my $nevents = 2000;
print "\nNumber of events in each segment is: $nevents\n\n";
sleep(1);

my $runDir = "/direct/phenix+hhj/arbint/JetSimulation/EnergyScaleChecks/CuAu";
print "Writing out events to dir: $runDir\n\n";
sleep(5);

my $ii = 0;
for($ii=0;$ii<$segments;$ii++) {
    my $dir = "EnergyScaleChecks_$ii";
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


#Create EnergyScaleChecks.cmd script                                                                                                                                             
    print "Creating EnergyScaleChecks.cmd\n";
    open (ENERGYSCALECHECKS_CMD, ">$runDir/$dir/EnergyScaleChecks.cmd");
    print ENERGYSCALECHECKS_CMD "#!/bin/csh\n";
    print ENERGYSCALECHECKS_CMD "source /etc/csh.login\n";
    print ENERGYSCALECHECKS_CMD "foreach i (/etc/profile.d/\*.csh)\n";
    print ENERGYSCALECHECKS_CMD "     source \$i\n";
    print ENERGYSCALECHECKS_CMD "end\n";
    print ENERGYSCALECHECKS_CMD "source /opt/phenix/bin/phenix_setup.csh -n ana.475\n";
    print ENERGYSCALECHECKS_CMD "source /opt/phenix/bin/odbcini_setup.csh\n";
    print ENERGYSCALECHECKS_CMD "setenv LD_LIBRARY_PATH /direct/phenix\+u/arbint/Jets/Analysis/install/lib:\$\{LD_LIBRARY_PATH\}\n";
    print ENERGYSCALECHECKS_CMD "if (\$\?\_CONDOR_SCRATCH_DIR\) then\n";
    print ENERGYSCALECHECKS_CMD "set workDir = \$\{_CONDOR_SCRATCH_DIR\}\n";
    print ENERGYSCALECHECKS_CMD "endif\n";
    print ENERGYSCALECHECKS_CMD "cd \$workDir\n";
    print ENERGYSCALECHECKS_CMD "pwd\n";

#Set PYTHIA
    print "Setting PYTHIA\n\n";
    print ENERGYSCALECHECKS_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/runCondor/JetSimulation/EnergyScaleChecks/pythia.cfg .\n";
    print ENERGYSCALECHECKS_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/generatePythiaCommon.C .\n";
    print ENERGYSCALECHECKS_CMD "root -l -q -b \"generatePythiaCommon.C($nevents)\"\n";

#Set PISA
    print "Setting PISA\n";
    print ENERGYSCALECHECKS_CMD "ln -sf $runDir/$dir/glogon.kumac .\n";
    print ENERGYSCALECHECKS_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pisaLinker/* .\n";
    print ENERGYSCALECHECKS_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pisaToDSTLinker/* .\n";
    print ENERGYSCALECHECKS_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pcDeadArea/pad_deadroc_372524.dat pad_deadroc.dat\n";
    print ENERGYSCALECHECKS_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pcDeadArea/pad_deadch_372524.dat pad_deadch.dat\n";
    print ENERGYSCALECHECKS_CMD "pisa < pisa.input\n";

#Set PISA to DST
    print "Setting PISA to DST\n\n";
    print ENERGYSCALECHECKS_CMD "root -l -q -b \"pisaToDST.C($nevents, 372524)\" \n";

#Set the code
    print "Setting the code\n\n";
    print ENERGYSCALECHECKS_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/runEnergyScaleChecks.C .\n";
    print ENERGYSCALECHECKS_CMD "root -l -q -b \"runEnergyScaleChecks.C($nevents, true)\"\n";

    close(ENERGYSCALECHECKS_CMD);

#Submit the job
    system("cd $runDir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/runCondor/JetSimulation/EnergyScaleChecks/EnergyScaleChecks.job");
    system("cd $runDir/$dir; chmod a+x EnergyScaleChecks.cmd");
    system("cd $runDir/$dir; condor_submit EnergyScaleChecks.job");

#Sleep for 1 seconds
    sleep(1);
}
