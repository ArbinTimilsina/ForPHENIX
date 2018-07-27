#!/usr/local/bin/perl
use strict;
use DBI;

#Trigger name
#my $triggerName = "BBCLL1(>0 tubes) narrowvtx";
#my $triggerName = "ERT_4x4c&BBCLL1";
my $triggerName = "ERT_4x4c&BBCLL1(narrow)";

my $pdb = DBI->connect("dbi:ODBC:daq") || die $DBI::error;
my $ps = $pdb->prepare("
select
 distinct(trigger.runnumber),
 trigger.scaledown,
 trigger.scalerfscaled,
 trigger.scalerberscaled

from run, trigger

where
 run.runnumber>=432637
 and run.runnumber<=436647
 and (run.runtype='PHYSICS')
 and trigger.name='$triggerName'
 and run.runnumber=trigger.runnumber;
");

#Loop over the list
my $totalRuns = 0;
my $totalEvents = 0;

$ps->execute();
my @rows = ();
while(@rows = $ps->fetchrow_array())
{
    $totalRuns++;

    if($rows[2]){
    $totalEvents += $rows[2];
    }else{
    $totalEvents += $rows[3];
    }
}

printf("\nTotal runs: %i\n", $totalRuns);
printf("Total events: %.04E\n\n", $totalEvents);


$ps->finish();
$pdb->disconnect;
