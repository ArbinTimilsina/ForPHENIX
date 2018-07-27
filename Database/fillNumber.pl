#!/usr/local/bin/perl
use strict;
use DBI;

#Run list
#my $fileList = "/direct/phenix+u/arbint/Jets/Analysis/runFiles/goodRunListPP.txt";
my $fileList = "/direct/phenix+u/arbint/Jets/Analysis/runFiles/goodRunListCuAu.txt";

#Output name
#my $outputFile = "ppFill.txt";
my $outputFile = "CuAuFill.txt";

my $pdb = DBI->connect("dbi:ODBC:daq") || die $DBI::error;

my @rows = ();
my $ps = $pdb->prepare("
select
 distinct(run.runnumber),
 run.fillnumberblue

from run

where
 run.runnumber=?
 and (run.runtype='PHYSICS')
");

#open files
open FILELIST, "<$fileList" or die $!;

system("rm -rf $outputFile");

print "\nWriting the output in the format\n\n";
print "------------------------\n";
print "Run\tBlueFillNumber\n";
print "------------------------\n\n";

#Loop over the list
my $runNumber = 0;
my $blueFillNumber= 0;
while(my $line=<FILELIST>)
{
    chomp($line);

    $ps->execute($line);
while(@rows = $ps->fetchrow_array())
{
    $runNumber = $rows[0];
    $blueFillNumber = $rows[1];
}
    system("echo \'$runNumber\t$blueFillNumber\t\' >> $outputFile");
}
$ps->finish();

$pdb->disconnect;
