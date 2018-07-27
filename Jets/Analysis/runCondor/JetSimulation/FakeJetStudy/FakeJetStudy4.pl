#!/usr/local/bin/perl

use strict;

my $runDir = "/direct/phenix+hhj/arbint/JetSimulation/FakeJetStudy/FakeJetStudy4";

my $segments = 80000;

my $ii = 0;
my $nevents = 1500;
for($ii=60000;$ii<$segments;$ii++) {
    my $dir = "FakeJetStudy_$ii";
    print "\n***********************************\n";
    print "Making directory $dir\n";
    system("cd $runDir; rm -rf $dir; mkdir -p $dir");

#Make the sHijing.xml file
    print "Making the xml file for sHijing\n";

    my $SEED= $ii + 444444;
    print "Seed set to ${SEED}\n";

    open (SHIJING_XML, ">$runDir/$dir/sHijing.xml");
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
    print SHIJING_XML "   <R>0.2</R>\n";
    print SHIJING_XML "   <PID>2000000</PID>\n";
    print SHIJING_XML "  </Algorithm>\n";
    print SHIJING_XML " </FastJet>\n";
    print SHIJING_XML "</HIJING>\n";
    close (SHIJING_XML);

#Create glogon.kumac
    print "Creating glogon.kumac\n";
    open (GLOGON_KUMAC, ">$runDir/$dir/glogon.kumac");
    print GLOGON_KUMAC "phpythia 10 phhijing.root\n";
    print GLOGON_KUMAC "ptrig $nevents\n";
    print GLOGON_KUMAC "exit\n";
    close (GLOGON_KUMAC);

#Create FakeJetStudy.cmd script
    print "Creating FakeJetStudy.cmd\n";
    open (FAKEJETSTUDY_CMD, ">$runDir/$dir/FakeJetStudy.cmd");
    print FAKEJETSTUDY_CMD "#!/bin/csh\n";
    print FAKEJETSTUDY_CMD "source /etc/csh.login\n";
    print FAKEJETSTUDY_CMD "foreach i (/etc/profile.d/\*.csh)\n";
    print FAKEJETSTUDY_CMD "     source \$i\n";
    print FAKEJETSTUDY_CMD "end\n";
    print FAKEJETSTUDY_CMD "source /opt/phenix/bin/phenix_setup.csh -n ana.458\n";
    print FAKEJETSTUDY_CMD "source /opt/phenix/bin/odbcini_setup.csh\n";
    print FAKEJETSTUDY_CMD "setenv LD_LIBRARY_PATH /direct/phenix\+u/arbint/Jets/Analysis/install/lib:\$\{LD_LIBRARY_PATH\}\n";
    print FAKEJETSTUDY_CMD "set workDir = $runDir/$dir/test\n";
    print FAKEJETSTUDY_CMD "if (\$\?\_CONDOR_SCRATCH_DIR\) then\n";
    print FAKEJETSTUDY_CMD "set workDir = \$\{_CONDOR_SCRATCH_DIR\}\n";
    print FAKEJETSTUDY_CMD "endif\n";
    print FAKEJETSTUDY_CMD "cd \$workDir\n";
    print FAKEJETSTUDY_CMD "pwd\n";
    print FAKEJETSTUDY_CMD "ln -sf $runDir/$dir/sHijing.xml .\n";
    print FAKEJETSTUDY_CMD "sHijing\n";

#Convert HepMC output to PHPYTHIA output
    print "Setting up HepMCPyNodeReader\n";
    print FAKEJETSTUDY_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/runHepMCPyNodeReader.C .\n";
    print FAKEJETSTUDY_CMD "root -l -q -b runHepMCPyNodeReader.C\n";

#Set PISA
    print "Setting PISA\n";
    print FAKEJETSTUDY_CMD "ln -sf $runDir/$dir/glogon.kumac .\n";
    print FAKEJETSTUDY_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pisaLinker/* .\n";
    print FAKEJETSTUDY_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pisaToDSTLinker/* .\n";
    print FAKEJETSTUDY_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pcDeadArea/pad_deadroc_372524.dat pad_deadroc.dat\n";
    print FAKEJETSTUDY_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/pcDeadArea/pad_deadch_372524.dat pad_deadch.dat\n";
    print FAKEJETSTUDY_CMD "pisa <pisa.input\n";

#Set PISA to DST
    print "Setting PISA to DST\n";
    print FAKEJETSTUDY_CMD "root -l -q -b \"pisaToDST.C($nevents, 372524)\"\n";

#Set the code
    print "Setting the code\n";
    print FAKEJETSTUDY_CMD "ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/JetSimulation/simMacros/runFakeJetStudy.C .\n";
    print FAKEJETSTUDY_CMD "root -l -q -b runFakeJetStudy.C\n";

    close(FAKEJETSTUDY_CMD);

#Submit the job
    print "Submitting to condor\n";
    system("cd $runDir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/runCondor/JetSimulation/FakeJetStudy/FakeJetStudy.job");
    system("cd $runDir/$dir; chmod a+x FakeJetStudy.cmd");
    system("cd $runDir/$dir; condor_submit FakeJetStudy.job");

    print "Done!!\n";

#Sleep for 30 seconds
    sleep(30);
}


