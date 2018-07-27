#!/usr/local/bin/perl                                                                                                                              
use strict;

my $evntdir = "/direct/phenix+hhj/arbint/CondorOutput/CuAuEMCalMatching";

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
        if($line =~ /.*_run12_online_central-(\d+)-(\d+).*/){

        my $dir = "forEMCalMatching_$1-$2";
        print "Making directory $dir\n";
        system("cd $evntdir; mkdir -p $dir;");

        #a bit of debugging output...                                                                                                            
        print "The filename for the CNT is:\n";
        print "\tCNT_run12_online_central-".$1."-".$2.".root\n";
    
        system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Calibrations/code/macro/run_forEMCalMatching.C; ln -sf /direct/phenix+u/arbint/Calibrations/runCondor/CuAuEMCalMatching/forEMCalMatching.job");
        system("cd $evntdir/$dir; echo \'$line\' > $CNTname");
        system("cd $evntdir/$dir; echo \'#!/bin/csh\ncd $evntdir/$dir\nsetenv TSEARCHPATH /direct/phenix+u/arbint/Calibrations/install\nroot -l -q -b run_forEMCalMatching.C\' > forEMCalMatching.cmd; chmod a+x forEMCalMatching.cmd");
        system("cd $evntdir/$dir; condor_submit forEMCalMatching.job");

	sleep(10);
    }
}

