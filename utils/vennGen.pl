#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

veenGen.pl --file merged.vcf --patterns \"WHAM-.*,LUMPY-.*\" /
           --names WHAM,LUMPY > Rcode.txt

Description:

Generates the [R] code for a SV venn

Required:

--file,-f         - <STRING> - the vcf file to summarize.
--patterns,-p     - <STRING> - Comma separated list of pattern
--names,-n        - <STRING> - Comma seperated list of names (same order as patterns)

";


my ($help);
my $file    ;
my $patterns;
my $names   ;
my $opt_success = GetOptions('help'      => \$help,
			     "file=s"    => \$file,
			     "names=s"   => \$names,
			     "patterns=s" => \$patterns,
    );

die $usage if $help || ! $opt_success;

die $usage unless $file;
die $usage unless $patterns;
die $usage unless $names;

my @NAMES    = split /,/, $names;
my @PATTERNS = split /,/, $patterns;

open (my $IN, '<', $file) or die "Can't open $file for reading\n$!\n";

foreach my $n (@NAMES){
    print "$n\t";
}

print "type\tlength\n";

while (<$IN>) {
    chomp;
    next if $_ =~ /^\#/;
    my @line = split /\t/, $_;
    
    foreach my $p (@PATTERNS){
	if ($line[7] =~ /\Q$p/){
	    print "TRUE\t" ;
	}
	else{
	    print "FALSE\t" ;
	}
    }
    my $svlen  = $1 if $_ =~ /SVLEN=(.*?);/;
    my $svtype = $1 if $_ =~ /SVTYPE=(.*?);/;
    print "$svtype\t$svlen\n";
}

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

sub {

}

