#!/usr/local/bin/perl                                                                                                                              
use strict;

my $eventDir = "/direct/phenix+hhj2/arbint/JetSimWithoutDetector";

my $ii = 0;
my $segments = 55555;

my $nevents = 5555;
for($ii=0;$ii<$segments;$ii++) {
    my $dir = "JetSimWithoutDetector_$ii";
    print "\nMaking directory $dir\n";                                                                                                              
    system("cd $eventDir; rm -rf $dir; mkdir -p $dir");

#Make the sHijing.xml file
    print "Making the xml file for sHijing\n";

    my $SEED= $ii + 123456;
    print "Seed set to ${SEED}\n";

    open (SHIJING_XML, ">$eventDir/$dir/sHijing.xml");
    print SHIJING_XML "<?xml version=\"1.0\"?>\n";
    print SHIJING_XML "<HIJING>\n";
    print SHIJING_XML " <EFRM>200.0</EFRM>\n";
    print SHIJING_XML " <FRAME>CMS</FRAME>\n";
    print SHIJING_XML " <PROJ>A</PROJ>\n";
    print SHIJING_XML " <TARG>A</TARG>\n";
    print SHIJING_XML " <IAP>197</IAP>\n";
    print SHIJING_XML " <IZP>79</IZP>\n";
    print SHIJING_XML " <IAT>64</IAT>\n";
    print SHIJING_XML " <IZT>29</IZT>\n";
    print SHIJING_XML " <N>${nevents}</N>\n";
    print SHIJING_XML " <BMIN>0.0</BMIN>\n";
    print SHIJING_XML " <BMAX>20.0</BMAX>\n";
    print SHIJING_XML " <KEEP_SPECTATORS>1</KEEP_SPECTATORS>\n";
    print SHIJING_XML " <OUTPUT>sHijing.dat</OUTPUT>\n";
    print SHIJING_XML " <RANDOM>\n";
    print SHIJING_XML "  <SEED>${SEED}</SEED>\n";
    print SHIJING_XML " </RANDOM>\n";
    print SHIJING_XML " <IHPR2>\n";
    print SHIJING_XML "  <4>0</4>\n";
    print SHIJING_XML "  <21>1</21>\n";
    print SHIJING_XML " </IHPR2>\n";
    print SHIJING_XML " <FastJet>\n";
    print SHIJING_XML "  <Algorithm>\n";
    print SHIJING_XML "   <Name>AntikT</Name>\n";
    print SHIJING_XML "   <R>0.15</R>\n";
    print SHIJING_XML "   <PID>1500000</PID>\n";
    print SHIJING_XML "  </Algorithm>\n";
    print SHIJING_XML "  <Algorithm>\n";
    print SHIJING_XML "   <Name>AntikT</Name>\n";
    print SHIJING_XML "   <R>0.2</R>\n";
    print SHIJING_XML "   <PID>2000000</PID>\n";
    print SHIJING_XML "  </Algorithm>\n";
    print SHIJING_XML "  <Algorithm>\n";
    print SHIJING_XML "   <Name>AntikT</Name>\n";
    print SHIJING_XML "   <R>0.25</R>\n";
    print SHIJING_XML "   <PID>2500000</PID>\n";
    print SHIJING_XML "  </Algorithm>\n";
    print SHIJING_XML "  <Algorithm>\n";
    print SHIJING_XML "   <Name>AntikT</Name>\n";
    print SHIJING_XML "   <R>0.3</R>\n";
    print SHIJING_XML "   <PID>3000000</PID>\n";
    print SHIJING_XML "  </Algorithm>\n";
    print SHIJING_XML " </FastJet>\n";
    print SHIJING_XML "</HIJING>\n";
    close (SHIJING_XML);

#Create JetSimWithoutDetector.cmd script
    print "Creating JetSimWithoutDetector.cmd\n";
    open (JETSIMWITHOUTDETECTOR_CMD, ">$eventDir/$dir/JetSimWithoutDetector.cmd");
    print JETSIMWITHOUTDETECTOR_CMD "#!/bin/csh\n";
    print JETSIMWITHOUTDETECTOR_CMD "source /etc/csh.login\n";
    print JETSIMWITHOUTDETECTOR_CMD "foreach i (/etc/profile.d/\*.csh)\n";
    print JETSIMWITHOUTDETECTOR_CMD "     source \$i\n";
    print JETSIMWITHOUTDETECTOR_CMD "end\n";
    print JETSIMWITHOUTDETECTOR_CMD "source /opt/phenix/bin/phenix_setup.csh -n ana.458\n";
    print JETSIMWITHOUTDETECTOR_CMD "source /opt/phenix/bin/odbcini_setup.csh\n";
    print JETSIMWITHOUTDETECTOR_CMD "setenv LD_LIBRARY_PATH /direct/phenix\+u/arbint/Jets/Analysis/install/lib:\$\{LD_LIBRARY_PATH\}\n";
    print JETSIMWITHOUTDETECTOR_CMD "set workDir = $eventDir/$dir/test\n";
    print JETSIMWITHOUTDETECTOR_CMD "if (\$\?\_CONDOR_SCRATCH_DIR\) then\n";
    print JETSIMWITHOUTDETECTOR_CMD "set workDir = \$\{_CONDOR_SCRATCH_DIR\}\n";
    print JETSIMWITHOUTDETECTOR_CMD "endif\n";
    print JETSIMWITHOUTDETECTOR_CMD "cd \$workDir\n";
    print JETSIMWITHOUTDETECTOR_CMD "pwd\n";
    print JETSIMWITHOUTDETECTOR_CMD "ln -sf $eventDir/$dir/sHijing.xml .\n";
    print JETSIMWITHOUTDETECTOR_CMD "sHijing\n";

#Set Pythia
    print "Setting Pythia\n";
    print JETSIMWITHOUTDETECTOR_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/generatePythiaCommon.C .\n";
    print JETSIMWITHOUTDETECTOR_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pythiaWithoutDetector.cfg pythia.cfg\n";
    print JETSIMWITHOUTDETECTOR_CMD "root -l -q -b \"generatePythiaCommon.C($nevents)\"\n";

#Setthe code
    print "Setting the code\n";
    print JETSIMWITHOUTDETECTOR_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/runJetSimWithoutDetector.C\n";
    print JETSIMWITHOUTDETECTOR_CMD "root -l -q -b runJetSimWithoutDetector.C\n";

    close(JETSIMWITHOUTDETECTOR_CMD);

#Submit the job
    system("cd $eventDir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/runCondor/JetSimulation/JetSimWithoutDetector/JetSimWithoutDetector.job");
    system("cd $eventDir/$dir; chmod a+x JetSimWithoutDetector.cmd");
    system("cd $eventDir/$dir; condor_submit JetSimWithoutDetector.job");

#Sleep 15 seconds
    sleep(15);
}
