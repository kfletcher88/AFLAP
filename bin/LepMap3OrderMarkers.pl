#!/usr/bin/perl

################################################################################
#       A script to run LepMap3 Order Markers after generating a genotyp table with AFLAP.
################################################################################

use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_i $opt_L $opt_M $opt_d);
getopts('i:L:hd:M:');
use File::Basename;

if ( (defined($opt_h)) || !(defined($opt_i)) || !(defined($opt_L)) || !(defined($opt_M)) ) {
	print STDERR "\tThis script will run LerpMap3 Seperate Chromosomes as part of AFLAP.\n";
	print STDERR "\tRequired options are:\n";
	print STDERR "\t\t-i input Genotype table\n\t\t-L Linkage Group to Order\n\t\t-M Map file generated by LepMap3 Separate Chromosomes.\n\n";
	print STDERR "\tAdditional arguments:\n\t\t-d Location of LepMap3 directoy, default is ThirdParty/LepMap3/bin as expected by AFLAP.sh\n\n"; 
	exit;
}
my $Map = $opt_M if $opt_M;
my $input = $opt_i if $opt_i;
my $LG = $opt_L if $opt_L;
my $prefix = $Map;
$prefix =~ s{\.[^.]+$}{};
print STDERR "Input Genotype table is $input\n";
print STDERR "LG being process is $LG\n";
my $output = $prefix.'.LG'.$LG.'.txt';
my $error  = $prefix.'.LG'.$LG.'.stderr';
print STDERR "Results will be written to $output\n";

my $dir = dirname $0;

my $LepDir = $opt_d || $dir.'/../ThirdParty/LepMap3/bin';

print STDERR "Beginning LepMap3\n\n";
system("java -cp $LepDir OrderMarkers2 data=$input useMorgan=1 numMergeIterations=20 chromosome=$LG map=$Map > $output 2> $error");
print STDERR "LepMap3 complete for Chromosome $LG\n";
exit;
