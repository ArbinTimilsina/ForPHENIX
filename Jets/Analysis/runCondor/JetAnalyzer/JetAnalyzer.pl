#!/usr/local/bin/perl                                                                                                                              
use strict;

my $evntdir = "/direct/phenix+hhj2/arbint/JetAnalyzer";

#names of the filelists                                                                                                              
my $CNTname = "filelist.list";

#open the lists                                                                                                                      
open FILELIST, "<$CNTname" or die $!;

#loop over one of the lists
while(my $line=<FILELIST>)
{
       #get rid of the newline character
       chomp($line);

       #make sure the filename makes sense
        if($line =~ /.*_run12CuAu_200GeV_CA_pro99-(\d+)-(\d+).*/){

        my $dir = "JetAnalyzer_$1-$2";
        print "Making directory $dir\n";
        system("cd $evntdir; rm -rf $dir; mkdir -p $dir;");

        system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/runCondor/JetAnalyzer/run_JetAnalyzer.C"); 
        system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/runCondor/JetAnalyzer/JetAnalyzer.job"); 
        system("cd $evntdir/$dir; echo \'$line\' > $CNTname");
	system("cd $evntdir/$dir; echo \'#!/bin/csh\' >> JetAnalyzer.cmd");
	system("cd $evntdir/$dir; echo \'setenv DCACHE_DOOR phnxdoor1.rcf.bnl.gov:22133\' >> JetAnalyzer.cmd");
	system("cd $evntdir/$dir; echo \'setenv GSEARCHPATH \${GSEARCHPATH}:DCACHE\' >> JetAnalyzer.cmd");
	system("cd $evntdir/$dir; echo \'cd $evntdir/$dir\' >> JetAnalyzer.cmd");

        system("cd $evntdir/$dir; echo \'root -l -q -b run_JetAnalyzer.C\' >> JetAnalyzer.cmd; chmod a+x JetAnalyzer.cmd");
        system("cd $evntdir/$dir; condor_submit JetAnalyzer.job");

	sleep(2);
    }
}

