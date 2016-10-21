#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

=desc
Parsing the output file of hamstr oneseq to create input file for phyloprofile tool
(20.10.2016)
=cut

sub usage {
    my $msg = shift;
    print "example: perl phyloProfile.pl -i oneseqOutFolder -o output.matrix\n";
    print "-i\thamstr oneseq output folder (*.extended.profile)\n";
    print "-o\tOutput file\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_t,$opt_o);
getopts('i:t:o:');

# sanity checks;
my $oneseqDir = ($opt_i) ? $opt_i : usage("ERROR: No oneseq output folder given\n");
my $out = ($opt_o) ? $opt_o : usage("ERROR: No output file given\n");

### MAIN
my @allOutFiles = glob("$oneseqDir/*.extended.profile");
unless(@allOutFiles){
	usage("ERROR: No extended.profile file found in $oneseqDir!\n");
}

my %taxaList;	# list of all taxa
my %allGenes;	# list of all genes
my %fas;	# $fas{$protID#$taxonID} = FAS_score

foreach my $file(@allOutFiles){
	open(IN,$file) || die "Cannot open $file!\n";
	my @in = <IN>;
	close (IN);

	### get taxon ID, ortho ID and FAS scores from each oneseq output file
	foreach my $line(@in){
		chomp($line);	# Arath01153|aquco_5436@218851@1|28258|1	0.99512260472
#		print $line,"\n";
		$line =~ s/\|/\t/g;		# Arath01153	aquco_5436@218851@1	28258	1	0.99512260472
		my @tmp = split(/\t/,$line);

		my $geneID = $tmp[0];
		my $fas = $tmp[@tmp-1];
		my @hit = split(/\@/,$tmp[1]);
		my $taxonID = $hit[1];
		my $hitID = $hit[0];
#		print "$geneID - $taxonID - $fas";<>;

		### save to %taxaList, %allGenes and %fas
		$taxaList{"ncbi$taxonID"} = 1;
		$fas{"$geneID#ncbi$taxonID"} = $hitID."#".$fas;
		$allGenes{$geneID} = 1;
	}
}

### print output
open(OUT,">$out") || die "Cannot create $out!\n";
my @allTaxa = sort keys %taxaList;
my $allTaxa = join("\t",@allTaxa);
print OUT "geneID\t$allTaxa\n";

foreach my $gene(sort keys %allGenes){
	print OUT $gene;
	foreach my $taxon(sort @allTaxa){
		if($fas{"$gene#$taxon"}){
			print OUT "\t",$fas{"$gene#$taxon"};
		} else {
			print OUT "\t","NA";
		}
	}
	print OUT "\n";
}
close (OUT);

print "Finished! Check output file\n\t$out\n";

exit;















