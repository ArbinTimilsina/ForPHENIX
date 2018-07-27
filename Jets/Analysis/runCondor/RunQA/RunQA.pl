#!/usr/local/bin/perl                                                                                                                              
use strict;

my $evntdir = "/direct/phenix+hhj/arbint/CondorOutput/RunQA";

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

        my $dir = "RunQA_$1-$2";
        print "Making directory $dir\n";
        system("cd $evntdir; rm -rf $dir; mkdir -p $dir;");

        #a bit of debugging output...                                                                                                            
        print "The filename for the CNT is:\n";
        print "\tCNT_run12_online_central-".$1."-".$2.".root\n";
    
        system("cd $evntdir/$dir; ln -sf /direct/phenix+u/arbint/Jets/Analysis/code/macros/run_RunQA.C; ln -sf /direct/phenix+u/arbint/Jets/Analysis/runCondor/RunQA/RunQA.job");
        system("cd $evntdir/$dir; echo \'$line\' > $CNTname");
        system("cd $evntdir/$dir; echo \'#!/bin/csh\ncd $evntdir/$dir\nroot -l -q -b run_RunQA.C\' > RunQA.cmd; chmod a+x RunQA.cmd");
        system("cd $evntdir/$dir; condor_submit RunQA.job");

	sleep(10);
    }
}

