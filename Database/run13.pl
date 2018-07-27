#!/usr/local/bin/perl
use strict;
use DBI;

my $pdb = DBI->connect("dbi:ODBC:daq") || die $DBI::error;
my $ps = $pdb->prepare("
select
 run.runnumber, 
 run.runtype, 
 run.eventsinrun

from run

where
 run.runnumber<=398149
 and run.runnumber>=386773
 and run.runtype='PHYSICS'
 and run.eventsinrun > 34e6
 and (erunixtime-brunixtime)>5340;
");

#Loop over the list
print "\nWriting the output in the format\n\n";
print "------------------------\n";
print "Run number \t Events  \n";
print "------------------------\n\n";

$ps->execute();
my @rows = ();
while(@rows = $ps->fetchrow_array())
{
    printf("%i \t %i \n", $rows[0], $rows[2]);

}

$ps->finish();
$pdb->disconnect;
