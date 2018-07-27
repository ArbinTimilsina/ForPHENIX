#!/usr/local/bin/perl
use strict;

my $inputFile = "pp.txt";

open INPUTFILE, "<$inputFile" or die $!;

my $runNumber = 0;
my $scaleDownMB = 0;
my $scaleDownERT = 0;
my $finalEvents = 0;
my $beforeEndRunEvents = 0;

my $totalEvents = 0;
my $totalRuns = 0;
my $totalWeightedScaleDownMB = 0;
my $totalWeightedScaleDownERT = 0;

my $totalProperEvents = 0;

while(<INPUTFILE>)
  {
    chomp;
    ($runNumber,$scaleDownMB,$scaleDownERT,$finalEvents,$beforeEndRunEvents)=split("\t");

    if($finalEvents !=0){
    $totalEvents +=$finalEvents;
    $totalWeightedScaleDownMB += ($finalEvents * $scaleDownMB);
    $totalWeightedScaleDownERT += ($finalEvents * $scaleDownERT);
    $totalProperEvents += $finalEvents * ((1 + $scaleDownMB)/(1 + $scaleDownERT));
    }else{
	$totalEvents +=$beforeEndRunEvents;
	$totalWeightedScaleDownMB += ($beforeEndRunEvents * $scaleDownMB);
	$totalWeightedScaleDownERT += ($beforeEndRunEvents * $scaleDownERT);
	$totalProperEvents += $beforeEndRunEvents * ((1 + $scaleDownMB)/(1 + $scaleDownERT));
    }
    $totalRuns++;
  }
close(INPUTFILE,);

print "\n-----------------------------------\n";
printf("Total runs analyzed: %i\n\n", $totalRuns++);

printf("Total events analyzed: %.05E\n\n", $totalEvents);

printf("Total weighted scaledown (MB): %.0f\n\n", $totalWeightedScaleDownMB);
printf("Total weighted scaledown (ERT): %.0f\n\n", $totalWeightedScaleDownERT);

printf("Average scaledown (MB): %.3f\n", $totalWeightedScaleDownMB/$totalEvents);
printf("Average scaledown (ERT): %.3f\n", $totalWeightedScaleDownERT/$totalEvents);

printf("Total events proper: %.05E\n\n", $totalProperEvents);
printf("Total events proper (with vertex cut): %.05E\n\n", $totalProperEvents*(0.57));

print "-----------------------------------\n\n";



