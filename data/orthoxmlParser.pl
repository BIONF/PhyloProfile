#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use IO::Handle;
use Cwd;

# convert orthoXML file into wide format
# 27.06.2017

sub usage {
    my $msg = shift;
    print "example: perl orthoxmlParser.pl -i data/demo/orthoxml_test.xml\n";
    print "-i\torthoXML input file\n";
    die $msg."\n";
}

our($opt_i);
getopts('i:');

my $file = ($opt_i) ? $opt_i : usage("ERROR: No input file given\n");
open(IN,$file) || die "cannot open $file!!\n";
my @file = <IN>;
close (IN);
my $input = join("",@file);

# join lines
$input =~ s/\n/\t/g; $input =~ s/\s{2}/\t/g;
$input =~ s/\t//g;

### get taxonomy IDs and protein IDs for all species in xml file
my %taxon;	# $taxon{protID} = taxonID
my %orthoGene; # $orthoGeneID{$id} = protID

while($input =~ /<species(.)+?\/species>/){
	my $hit = $&;
	substr($input,index($input,$hit),length($hit),"");
#	print $hit;<>;

	# get taxonID
	$hit =~ /NCBITaxId=\"(\d+)\"/;
	my $taxonID = $1;

	# get protIDs
	while($hit =~ /gene id=(.)*?\/>/){
		my $hitGene = $&;
		substr($hit,index($hit,$hitGene),length($hitGene),"");

		$hitGene =~ /id=\"(.)+?\"/;
		my $orthoGeneID = $&; $orthoGeneID =~ s/id=//; $orthoGeneID =~ s/\"//g;
		$hitGene =~ /geneId=\"(.)+?\"/;
		my $geneID = $&; $geneID =~ s/geneId=//; $geneID =~ s/\"//g;
		$hitGene =~ /protId=\"(.)+?\"/;
		my $protID = $&; $protID =~ s/protId=//; $protID =~ s/\"//g;


		if(length($protID) < 1){
			$protID = $geneID;
		}
#		print $taxonID,"-",$orthoGeneID,"-",$geneID,"-",$protID;<>;

		$taxon{$protID} = $taxonID;
		$orthoGene{$orthoGeneID} = $protID;
	}
}

### get score IDs
my %scoreID; # $scoreID{id} = scoreDef
while($input =~ /<scoreDef(.)+?\/>/){
	my $hit = $&;
	substr($input,index($input,$hit),length($hit),"");

	$hit =~ /id=\"(.)+?\"/;
	my $id = $&; $id =~ s/id=//; $id =~ s/\"//g;
	$hit =~ /desc=\"((.)+)?\"/;
	my $def = $1;
#	print $id," - ",$def;<>;

	$scoreID{$id} = $def;
}

# check existing scores for ortho genes
$input =~ /<groups(.)+?\/groups>/;
my $orthoGroupBlock = $&;
my %flag; # $flag{$scoreID} = 0/1   ; 0 if not exist in geneRef blocks
foreach my $scoreID (keys %scoreID){
	$orthoGroupBlock =~ /geneRef(.)+?<\/geneRef/;
	my $hit = $&;

	if($hit =~ /score id=\"$scoreID?\"/){
		$flag{$scoreID} = 1;
		next;
	}
}

foreach my $scoreID (keys %scoreID){
	unless($flag{$scoreID}){
		delete($scoreID{$scoreID});
	}
}

### parse ortholog group and output the results
print "geneID\tncbiID\torthoID";
foreach my $scoreID (sort keys %scoreID){
	print "\t$scoreID";
}
print "\n";

while($input =~ /<orthologGroup id=\"(.)+?\"(.)+?\/orthologGroup>/){
	#	print $1,"\n",$hit;<>;
	my $hit = $&;
	substr($input,index($input,$hit),length($hit),"");

	$hit =~ /<orthologGroup id=\"(.)+?\"/;
	my $groupID = $&; $groupID =~ s/<orthologGroup id=//; $groupID =~ s/"//g;
	unless($groupID =~ /\D/){
		$groupID = "OG".$groupID;
	}

	# get all proteins for this group
	while($hit =~ /<geneRef id=\"(.)+?\"(.)+?\/geneRef>/){

		my $hitGene = $&;
		substr($hit,index($hit,$hitGene),length($hitGene),"");

		$hitGene =~ /<geneRef id=\"(.)+?\"/;
		my $orthoGeneID = $&; $orthoGeneID =~ s/<geneRef id=//; $orthoGeneID =~ s/"//g;

		print $groupID,"\t","ncbi",$taxon{$orthoGene{$orthoGeneID}},"\t",$orthoGene{$orthoGeneID};

		# get scores for this gene
		my %score;
		while($hitGene =~ /<score id=\"(.)+?\/>/){
			my $hitScore = $&;
			substr($hitGene,index($hitGene,$hitScore),length($hitScore),"");

			$hitScore =~ /id=\"(.)+?\"/;
			my $scoreID = $&; $scoreID =~ s/id=//; $scoreID =~ s/\"//g;
			$hitScore =~ /value=\"(.)+?\"/;
			my $value = $&; $value =~ s/value=//; $value =~ s/\"//g;

			$score{$scoreID} = $value;
		}

		# print scores
		foreach my $scoreID (sort keys %scoreID){
			if(defined $score{$scoreID}){
				print "\t$score{$scoreID}";
			} else {
				print "\tNA";
			}
		}
		print "\n";
	}
}

exit;
