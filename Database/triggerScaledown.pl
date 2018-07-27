#!/usr/local/bin/perl
use strict;
use DBI;

#Run list
my $fileList = "/direct/phenix+u/arbint/Jets/Analysis/runFiles/goodRunListPP.txt";

#Trigger name
my $mbTriggerName = "BBCLL1(>0 tubes) narrowvtx";
my $ertTriggerName = "ERT_4x4c&BBCLL1(narrow)";

#Output name
my $outputFile = "pp.txt";

#/direct/phenix+WWW/run/daq/runcontrol/RunSummary.php
#scalerf -> scalerber -> scalerupdate
#trigger.scalerberscaled

my $pdb = DBI->connect("dbi:ODBC:daq") || die $DBI::error;

my @rowsMB = ();
my $psMB = $pdb->prepare("
select
 distinct(trigger.runnumber),
 trigger.scaledown,
 trigger.scalerfscaled,
 trigger.scalerberscaled

from run, trigger

where
 run.runnumber=?
 and (run.runtype='PHYSICS')
 and trigger.name='$mbTriggerName'
 and run.runnumber=trigger.runnumber;
");

my @rowsERT = ();
my $psERT = $pdb->prepare("
select
 distinct(trigger.runnumber),
 trigger.scaledown,
 trigger.scalerfscaled,
 trigger.scalerberscaled

from run, trigger

where
 run.runnumber=?
 and (run.runtype='PHYSICS')
 and trigger.name='$ertTriggerName'
 and run.runnumber=trigger.runnumber;
");

#open files
open FILELIST, "<$fileList" or die $!;

system("rm -rf $outputFile");

print "\nWriting the output in the format\n\n";
print "----------------------------------------------------------------------\n";
print "Run\tMB ScaleDown\tERT ScaleDown\tFinal\tBeforeEndRun\n";
print "----------------------------------------------------------------------\n\n";

#Loop over the list
my $runNumberMB = 0;
my $scaleDownMB= 0;
my $scaleDownERT= 0;
my $finalEventsMB = 0;
my $beforeEndRunEventsMB = 0;

while(my $line=<FILELIST>)
{
    chomp($line);

    $psMB->execute($line);
while(@rowsMB = $psMB->fetchrow_array())
{
    $runNumberMB = $rowsMB[0];
    $scaleDownMB = $rowsMB[1];
    if($rowsMB[2]){
    $finalEventsMB = $rowsMB[2];
    }
    if($rowsMB[3]){
    $beforeEndRunEventsMB = $rowsMB[3];
    }
}

    $psERT->execute($line);
while(@rowsERT = $psERT->fetchrow_array())
{
    $scaleDownERT = $rowsERT[1];
}
    system("echo \'$runNumberMB\t$scaleDownMB\t$scaleDownERT\t$finalEventsMB\t$beforeEndRunEventsMB\' >> $outputFile");
}

$psMB->finish();
$psERT->finish();
$pdb->disconnect;
