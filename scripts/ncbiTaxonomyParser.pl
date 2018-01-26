#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

=desc
parse info from taxonomy NCBI dmp files (names.dmp and nodes.dmp)
to create a mapping file taxonNamesFull.txt that contains
ncbiID	fullName	rank	parentID
23.01.2017
=cut

sub usage {
    my $msg = shift;
    print "example: perl ncbiTaxonomyParser.pl -i dmp_folder\n";
    print "-i\tFolder contains DMP files (names.dmp and nodes.dmp)\n";
    die $msg."\n";
}

# global variables
our($opt_i);
getopts('i:');

# sanity checks;
my $input = ($opt_i) ? $opt_i : usage("ERROR: No input folder given\n");

### parse names and IDs from names.dmp
print "Parsing names.dmp...";

open(NAME,"$input/names.dmp") || die "Cannot open $input/names.dmp!\n";
my @nameIN = <NAME>;
close (NAME);

my %name;

foreach my $line(@nameIN){
	chomp($line);
	if($line =~ /scientific name/){
		my @tmp = split(/\t\|\t/,$line);	# taxID - name
		my $name = $tmp[1];
		$name =~ s/\'//g;
		$name{$tmp[0]} = $name;
	}
}
print " Finished!\n";

### parse IDs and parents IDs
print "Parsing nodes.dmp...";

open(NODE,"$input/nodes.dmp") || die "Cannot open $input/nodes.dmp!\n";
my @nodeIN = <NODE>;
close (NODE);

open(OUT,">$input/taxonNamesFull.txt");
print OUT "ncbiID	fullName	rank	parentID\n";

my $c = 0;
foreach my $line(@nodeIN){
	chomp($line);
	my @tmp = split(/\t\|\t/,$line);	# taxID - parentID - rank

	## output: ncbiID	fullName	rank	parentID
	my $rank = $tmp[2];
	$rank =~ s/\s//g; $rank =~ s/\'//g;
	print OUT "$tmp[0]\t$name{$tmp[0]}\t$rank\t$tmp[1]\n";

	$c ++;
	print $c,"/",scalar(@nodeIN),"\n";
}
print " Finished!\n Check your output at\n$input/taxonNamesFull.txt\n";
close (OUT);

exit;
