#!/usr/local/bin/perl                                                                                                                              
use strict;

my $segments = 10000;
print "\nRunning over $segments segments.\n\n";
sleep(1);

my $nevents = 500;
print "\nNumber of events in each segment is: $nevents\n\n";
sleep(1);

my $mainDir = "/phenix/hhj/arbint/JetSimulation/CentralityStudy/VertexPlus10";
print "Writing out events to dir: $mainDir\n\n";
sleep(2);

my $ii = 0;
for($ii=0;$ii<$segments;$ii++) {
    my $dir = "CentralityStudy_$ii";
    print "\n************************************\n";
    print "Making directory $dir\n\n";
    system("cd $mainDir; rm -rf $dir; mkdir -p $dir");

#Create glogon.kumac
    print "Creating glogon.kumac\n";
    open (GLOGON_KUMAC, ">$mainDir/$dir/glogon.kumac");
    print GLOGON_KUMAC "phpythia 1 phhijing.root\n";
    print GLOGON_KUMAC "ptrig $nevents\n";
    print GLOGON_KUMAC "exit\n";
    close (GLOGON_KUMAC);

#Create CentralityStudy.cmd
    print "Creating CentralityStudy.cmd\n";
    open (CENTRALITYSTUDY_CMD, ">$mainDir/$dir/CentralityStudy.cmd");
    print CENTRALITYSTUDY_CMD "#!/bin/csh\n";
    print CENTRALITYSTUDY_CMD "source /etc/csh.login\n";
    print CENTRALITYSTUDY_CMD "foreach i (/etc/profile.d/\*.csh)\n";
    print CENTRALITYSTUDY_CMD "     source \$i\n";
    print CENTRALITYSTUDY_CMD "end\n";
    print CENTRALITYSTUDY_CMD "source /opt/phenix/bin/phenix_setup.csh -n ana.475\n";
    print CENTRALITYSTUDY_CMD "source /opt/phenix/bin/odbcini_setup.csh\n";
    print CENTRALITYSTUDY_CMD "setenv LD_LIBRARY_PATH /direct/phenix\+u/arbint/Jets/Analysis/install/lib:\$\{LD_LIBRARY_PATH\}\n";
    print CENTRALITYSTUDY_CMD "if (\$\?\_CONDOR_SCRATCH_DIR\) then\n";
    print CENTRALITYSTUDY_CMD "set workDir = \$\{_CONDOR_SCRATCH_DIR\}\n";
    print CENTRALITYSTUDY_CMD "endif\n";
    print CENTRALITYSTUDY_CMD "cd \$workDir\n";
    print CENTRALITYSTUDY_CMD "pwd\n";

#Set HIJING
    print "Setting HIJING\n\n";
    print CENTRALITYSTUDY_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/centMacros/phhijing.cfg .\n";
    print CENTRALITYSTUDY_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/centMacros/generateHijing.C .\n";
    print CENTRALITYSTUDY_CMD "root -l -q -b \"generateHijing.C($nevents)\"\n";

#Set PISA
    print "Setting PISA\n";
    print CENTRALITYSTUDY_CMD "ln -sf $mainDir/$dir/glogon.kumac .\n";
    print CENTRALITYSTUDY_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/centMacros/pisaLinker/* .\n";
    print CENTRALITYSTUDY_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/centMacros/event.par_vertexPlus10 event.par\n";
    print CENTRALITYSTUDY_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/centMacros/pisaToDSTLinker/* .\n";
    print CENTRALITYSTUDY_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pcDeadArea/pad_deadroc_372524.dat pad_deadroc.dat\n";
    print CENTRALITYSTUDY_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pcDeadArea/pad_deadch_372524.dat pad_deadch.dat\n";
    print CENTRALITYSTUDY_CMD "pisa < pisa.input\n";

#Set PISA to DST
    print "Setting PISA to DST\n\n";
    print CENTRALITYSTUDY_CMD "root -l -q -b \"pisaToDST.C($nevents, 372524)\"\n";

#Set the code
    print "Setting the code\n\n";
    print CENTRALITYSTUDY_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/centMacros/runCentralityStudy.C .\n";
    print CENTRALITYSTUDY_CMD "root -l -q -b \"runCentralityStudy.C\"\n";
    close(CENTRALITYSTUDY_CMD);

#Submit the job
    system("cd $mainDir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/runCondor/JetSimulation/CentralityStudy/CentralityStudy.job");
    system("cd $mainDir/$dir; chmod a+x CentralityStudy.cmd");
    system("cd $mainDir/$dir; condor_submit CentralityStudy.job");

#Sleep for 10 seconds
    sleep(10);
}
