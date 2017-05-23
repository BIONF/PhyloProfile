#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use IO::Handle;
use Cwd;

### get FULL taxonomy IDs for list of species names
### 19.05.2016
### and create a reduced list of species name
### 30.05.2016
### embedded into R script
### 17.06.2016
### including noranks
### 17.11.2016

sub usage {
    my $msg = shift;
    print "example: perl getTaxonomyInfo.pl -i taxonID_list -n /home/vinh/Desktop/data/project/taxonomy/taxonNamesFull.txt -a /home/vinh/Desktop/data/project/taxonomy/newTaxa.txt -o outputDir\n";
    print "-i\tTaxonID list\n";
    print "-n\tFile contains all NCBI IDs, their names and ranks\n";
    print "-a\tFile contains newly added taxa IDs, their names and ranks\n";
    print "-o\tOutput dir\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_n,$opt_o,$opt_a);
getopts('i:n:o:a:');

my $idList = ($opt_i) ? $opt_i : usage("ERROR: No input ID list given\n");
my $nameIN =  ($opt_n) ? $opt_n : usage("ERROR: No taxonNamesFull file given\n");
my $nameNewIN = ($opt_a) ? $opt_a : usage("ERROR: No newly added taxa file given\n");
my $outDir =  ($opt_o) ? $opt_o : usage("ERROR: No output dir given\n");

# check file exists
#open(IN,"$idList") || die "Cannot open $idList!!\n";
#my @in = <IN>;
#close (IN);

open(NAME,"$nameIN") || die "Cannot open $nameIN!!\n";
my @nameIn = <NAME>;
close (NAME);

#### get list of all species name and their NCBI taxon id
my %id;		# $id{NAME} = ID
my %name;	# $name{$id} = NAME
my %rank;	# $rank{$id} = RANK
my %parent;	# $parent{$id} = PARENT_ID
foreach my $line(@nameIn){
	chomp($line);
	my @tmp = split(/\t/,$line);	# id  name  rank  parentID
	$id{$tmp[1]} = $tmp[0];
	$name{$tmp[0]} = $tmp[1];
	$rank{$tmp[0]} = $tmp[2];
	$parent{$tmp[0]} = $tmp[3];
}

#### get the info for newly added taxa (if necessary)
open(NEW,"$nameNewIN") || die "Cannot open $nameNewIN!!\n";
my @nameNEW = <NEW>;
close (NEW);

if(scalar @nameNEW > 1){
	foreach my $line(@nameNEW){
		chomp($line);
    if(length $line > 2){
      my @tmp = split(/\t/,$line);	# id  name  rank  parentID
      $id{$tmp[1]} = $tmp[0];
      $name{$tmp[0]} = $tmp[1];
      $rank{$tmp[0]} = $tmp[2];
      $parent{$tmp[0]} = $tmp[3];
    }
	}
}

#print "Parsing taxonNamesFull.txt done!!\n";<>;

### GET LIST OF ALL TAXA from input file

#=for testing
#$idList = "geneID	ncbi272557	ncbi9360";	# example for input a genus ID instead of species or strain
#$idList = "geneID	ncbi4837";#	ncbi4932	ncbi3702	ncbi9606";
#$idList = "geneID	ncbi436017";
#$idList = "geneID	ncbi1096996	ncbi2000001	ncbi2000002	ncbi2000003	ncbi202950	ncbi62977";
#=cut

my @allTaxa = split(/\t/,$idList);
shift(@allTaxa);	# remove "geneID" tag

### get all ranking IDs for each taxon
my %info;	# $info{$ncbi}
my %reduceSpec;		# $reducedSpec{specID} = ncbiID	fullName	rank	parentID
foreach my $taxon(@allTaxa){
	my $ncbiID = $taxon;
	$ncbiID =~ s/ncbi//;

	unless($name{$ncbiID}){
		print "$ncbiID not in taxonNamesFull. CHECK AGAIN\n";
	} else {
		my $name = $name{$ncbiID};
		my $rank = $rank{$ncbiID};
		my $parentID = $parent{$ncbiID};
		my $parentRank = $rank{$parentID};

#		print "HERE: $ncbiID - $rank\n$parentID - $parentRank\n";	############ TESTING
 		$info{"$ncbiID"} = "$ncbiID#strain\t"."$ncbiID#$rank\t"."$parentID#$parentRank";
 		$reduceSpec{$ncbiID} = $ncbiID."\t".$name{$ncbiID}."\t".$rank."\t".$parentID;
 		$reduceSpec{$parentID} = $parentID."\t".$name{$parentID}."\t".$parentRank."\t".$parent{$parentID};

		unless($parentID == 1){
			do{
				$parentID = $parent{$parentID};
				$parentRank = $rank{$parentID};
#				print $parentID," - ",$parentRank," - ",$name{$parentID},"\n";	############ TESTING
				$info{"$ncbiID"} .= "\t"."$parentID#$parentRank";
				$reduceSpec{$parentID} = $parentID."\t".$name{$parentID}."\t".$parentRank."\t".$parent{$parentID};
			} until ($parentID == 1);
		}
#		print "$info{$ncbiID}\n";
#		print "NEXT...\n";		############ TESTING
#		<>;						############ TESTING
	}
}

##### create taxonNamesReduced.txt
open(NAMELIST,">$outDir/taxonNamesReduced.txt");
print NAMELIST "ncbiID	fullName	rank	parentID\n";

foreach(keys %reduceSpec){
	print NAMELIST $reduceSpec{$_},"\n";
}
close (NAMELIST);

### create output matrix
open(OUT,">$outDir/taxonID.list.fullRankID");
print OUT "No.\tabbrName\tncbiID\tfullName";

my @allRefRank =
(
"strain","norank","forma","subspecies","varietas","norank","norank","norank",
"species","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"speciessubgroup","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"speciesgroup","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"subgenus","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"genus","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"subtribe","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"tribe","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"subfamily","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"family","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"superfamily","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"parvorder","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"infraorder","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"suborder","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"order","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"superorder","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"infraclass","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"subclass","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"class","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"superclass","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"subphylum","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"phylum","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"superphylum","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"subkingdom","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"kingdom","norank","norank","norank","norank","norank","norank","norank","norank","norank","norank",
"superkingdom","norank","norank"
);

### print all ranks of refTaxon into output
#my @allRefRank = split(/\t/,$info{$refTaxon});	# all ranks of refTaxon
foreach(@allRefRank){
#	my @tmp = split(/#/,$_);
	print OUT "\t",$_;
}
print OUT "\n";

for(my $t = 0; $t < scalar(@allTaxa); $t++){
		my $ncbiID = $allTaxa[$t];
		$ncbiID =~ s/ncbi//;
		my @rankInfo = split(/\t/,$info{$ncbiID});
		print OUT $t+1,"\t",$allTaxa[$t],"\t",$ncbiID,"\t",$name{$ncbiID};

		for(my $i=0, my$j=0; $i < scalar(@allRefRank); $i++, $j++){
#			print "i=$i; j=$j (max j = ",scalar(@rankInfo),")\n";
			### rank name of refTaxon
#			my @rankRefTMP = split(/#/,$allRefRank[$i]);
#			my $rankRefName = $rankRefTMP[1];
			my $rankRefName = $allRefRank[$i];
#			print "$rankRefName\n";	############ TESTING

			### rank name of current taxon
			my @rankTMP = split(/#/,$rankInfo[$j]);
			my $rankName = $rankTMP[1];

			### if rankName = rankRefName, print ID and go to next rank ($i++, $j++)
			if($rankName eq $rankRefName){
#				print "SAME: $rankInfo[$j]\n";	############ TESTING
				my @rankInfoJ = split(/#/,$rankInfo[$j]);
				print OUT "\t$rankInfoJ[0]";
			}
			### else, increase $i, $j stays the same
			else {
#				print "noID: $rankInfo[$j] (previous: $rankInfo[$j-1])\n";	############ TESTING
				my @rankInfoJp = split(/#/,$rankInfo[$j-1]);
				print OUT "\t$rankInfoJp[0]";
				$j--;
			}
#			<>;
		}
		print OUT "\n";
}

close (OUT);
exit;
