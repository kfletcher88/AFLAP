#!/usr/bin/perl

################################################################################
#       A script to run LepMap3 Seperate Chromosomes after generating a genotyp table with AFLAP.
################################################################################

use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_i $opt_l $opt_t $opt_d);
getopts('i:l:hd:t:');
use File::Basename;

if ( (defined($opt_h)) || !(defined($opt_i)) || !(defined($opt_l)) ) {
	print STDERR "\tThis script will run LerpMap3 Seperate Chromosomes as part of AFLAP.\n";
	print STDERR "\tRequired options are:\n";
	print STDERR "\t\t-i input Genotype table\n\t\t-l minimum LOD score\n";
	print STDERR "\tNote, if multiple LOD limits are to be tested, then please provide these to the AFLAP script. AFLAP will spawn one job for each lod limit\n\n";
	print STDERR "\tAdditional arguments:\n\t\t-t Number of Threads to use\n\t\t-d Location of LepMap3 directoy\n\n"; 
	exit;
}
my $Threads = $opt_t || "1";
my $input = $opt_i if $opt_i;
my $LOD = $opt_l if $opt_l;
my $prefix = basename $input;
$prefix =~ s{\.[^.]+$}{};
print STDERR "Input Genotype table is $input\n";
print STDERR "Minimum LOD score is $LOD\n";
my $output = 'AFLAP_Results/LOD'.$LOD.'/'.$prefix.'.LOD'.$LOD.'.txt';
my $error  = 'AFLAP_Results/LOD'.$LOD.'/'.$prefix.'.LOD'.$LOD.'.stderr';
print STDERR "Results will be written to $output\n";

my $dir = dirname $0;
print "AFLAP scripts found:\n$dir";

my $LepDir = $opt_d || $dir.'/../ThirdParty/LepMap3/bin';

print STDERR "Beginning LepMap3\n\n";
system("mkdir -p AFLAP_Results/LOD$LOD");
system("java -cp $LepDir SeparateChromosomes2 data=$input lodLimit=$LOD numThreads=$Threads > $output 2> $error");
print STDERR "LepMap3 complete with LOD limit = $LOD";
exit;

