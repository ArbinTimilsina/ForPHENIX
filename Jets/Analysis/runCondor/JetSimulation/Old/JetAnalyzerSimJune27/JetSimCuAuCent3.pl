#!/usr/local/bin/perl

use strict;

my $nevents = 2200;
print "\nNumber of events in each segment is: $nevents\n\n";
sleep(1);

my $R = 0.2;
print "R-parameter is $R\n\n";
sleep(1);

my $minPt = 6.59;
print "Min pt is: $minPt (GeV/c)\n\n";
sleep(1);

my $minPtData = 6.59;

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

my $counter = 0.0;

my $mainDir = "/direct/phenix+hhj/arbint/JetSimulation/JetAnalyzerSim/CuAuCent3";

#names of the filelists
my $CNTname = "/direct/phenix+u/arbint/Jets/Analysis/runFiles/filelist.list";

#open the lists
open FILELIST, "<$CNTname" or die $!;

#loop over one of the lists
while(my $line=<FILELIST>)
{
    #get rid of the newline character
    chomp($line);

    #Get run number
        if($line =~ /.*-(\d+)-(\d+).*/){
	    my $runDir = "JetAnalyzerSim-".$1."-".$2."";
	    print "\n************************************\n";
	    print "Making directory $runDir\n\n";
	    system("cd $mainDir; rm -rf $runDir; mkdir -p $runDir");
sleep(5);

my $segments = 35;
print "\nRunning over $segments segments.\n\n";
sleep(1);

#For Data events
my $skipEvents= 0;

my $ii = 0;
for($ii=0;$ii<$segments;$ii++) {
    my $dir = "JetAnalyzerSim_$ii";
    print "\n************************************\n";
    print "Making directory $dir\n\n";

    $counter++;
    print "Total count is: $counter\n\n";

    system("cd $mainDir/$runDir; rm -rf $dir; mkdir -p $dir");
    system("cd $mainDir/$runDir/$dir; echo \'$line\' > filelist.list");

#Create glogon.kumac 
    print "Creating glogon.kumac\n";
    open (GLOGON_KUMAC, ">$mainDir/$runDir/$dir/glogon.kumac");
    print GLOGON_KUMAC "phpythia 1 phpythia.root\n";
    print GLOGON_KUMAC "ptrig $nevents\n";
    print GLOGON_KUMAC "exit\n";
    close (GLOGON_KUMAC);

#Create JetAnalyzerSim.cmd script
    print "Creating JetAnalyzerSim.cmd\n";
    open (JETANALYZERSIM_CMD, ">$mainDir/$runDir/$dir/JetAnalyzerSim.cmd");
    print JETANALYZERSIM_CMD "#!/bin/csh\n";
    print JETANALYZERSIM_CMD "source /etc/csh.login\n";
    print JETANALYZERSIM_CMD "foreach i (/etc/profile.d/\*.csh)\n";
    print JETANALYZERSIM_CMD "     source \$i\n";
    print JETANALYZERSIM_CMD "end\n";
    print JETANALYZERSIM_CMD "source /opt/phenix/bin/phenix_setup.csh -n ana.458\n";
    print JETANALYZERSIM_CMD "source /opt/phenix/bin/odbcini_setup.csh\n";
    print JETANALYZERSIM_CMD "setenv DCACHE_DOOR phnxdoor1.rcf.bnl.gov:22133\n";
    print JETANALYZERSIM_CMD "setenv GSEARCHPATH \${GSEARCHPATH}:DCACHE\n";
    print JETANALYZERSIM_CMD "setenv LD_LIBRARY_PATH /direct/phenix\+u/arbint/Jets/Analysis/install/lib:\$\{LD_LIBRARY_PATH\}\n";
    print JETANALYZERSIM_CMD "set workDir = $runDir/$dir/test\n";
    print JETANALYZERSIM_CMD "if (\$\?\_CONDOR_SCRATCH_DIR\) then\n";
    print JETANALYZERSIM_CMD "set workDir = \$\{_CONDOR_SCRATCH_DIR\}\n";
    print JETANALYZERSIM_CMD "endif\n";
    print JETANALYZERSIM_CMD "cd \$workDir\n";
    print JETANALYZERSIM_CMD "pwd\n";

#Set Run 12 Cu+Au data
    print "Setting Run 12 Cu+Au data\n";
    $skipEvents = $ii*10*$nevents;
    print "Skipping to event $skipEvents \n\n";

    print JETANALYZERSIM_CMD "ln -sf $mainDir/$runDir/$dir/filelist.list .\n";
    print JETANALYZERSIM_CMD "ln -sf  /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/generateData.C .\n";
    print JETANALYZERSIM_CMD "root -l -q -b \"generateData.C($R, $nc, $minPtData, $minCf, $maxCf, $nevents, $skipEvents, 3)\"\n";

#Set PYTHIA
    print "Setting PYTHIA\n\n";
    print JETANALYZERSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pythia.cfg .\n";
    print JETANALYZERSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/generatePythia.C .\n";
    print JETANALYZERSIM_CMD "root -l -q -b \"generatePythia.C($R, $ncPythia, $minPt, $minCfPythia, $maxCfPythia, $nevents)\"\n";

#Set PISA
    print "Setting PISA\n";
    print JETANALYZERSIM_CMD "ln -sf $mainDir/$runDir/$dir/glogon.kumac .\n";
    print JETANALYZERSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pisaLinker/* .\n";
    print JETANALYZERSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pisaToDSTLinker/* .\n";
    print JETANALYZERSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pcDeadArea/pad_deadroc_372524.dat pad_deadroc.dat\n";
    print JETANALYZERSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pcDeadArea/pad_deadch_372524.dat pad_deadch.dat\n";
    print JETANALYZERSIM_CMD "pisa < pisa.input\n";

#Set PISA to DST
    print "Setting PISA to DST\n\n";
    print JETANALYZERSIM_CMD "root -l -q -b \"pisaToDST.C($nevents, 372524)\"\n";

#Set the code
    print "Setting the code\n\n";
    print JETANALYZERSIM_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/runJetAnalyzerSimCuAu.C .\n";
    print JETANALYZERSIM_CMD "root -l -q -b \"runJetAnalyzerSimCuAu.C($R, $minPt, $nc, $minCf, $maxCf, $ncPythia, $minCfPythia, $maxCfPythia, 0.2, $nevents)\"\n";
    close(JETANALYZERSIM_CMD);

#Submit the job
    system("cd $mainDir/$runDir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/runCondor/JetSimulation/JetAnalyzerSim/JetAnalyzerSim.job");
    system("cd $mainDir/$runDir/$dir; chmod a+x JetAnalyzerSim.cmd");
    system("cd $mainDir/$runDir/$dir; condor_submit JetAnalyzerSim.job");

#Sleep for 20 seconds
    sleep(20);
   }
  }
}





