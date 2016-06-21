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

sub usage {
    my $msg = shift;
    print "example: perl getTaxonomyInfo.pl -i taxonID_list -n /home/vinh/Desktop/data/project/taxonomy/taxonNamesFull.txt -o outputDir\n";
    print "-i\taxonID list\n";
    print "-n\tFile contains all NCBI IDs, their names and ranks\n";
    print "-o\tOutput dir\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_n,$opt_o);
getopts('i:n:o:');

my $idList = ($opt_i) ? $opt_i : usage("ERROR: No input ID list given\n");
my $nameIN =  ($opt_n) ? $opt_n : usage("ERROR: No taxonNamesFull file given\n");
my $outDir =  ($opt_o) ? $opt_o : usage("ERROR: No output dir given\n");

# check file exists
#open(IN,"$idList") || die "Cannot open $idList!!\n";
#my @in = <IN>;
#close (IN);

open(NAME,"$nameIN") || die "Cannot open $nameIN!!\n";
my @nameIn = <NAME>;
close (NAME);

#### get list of all species name and their NCBI taxon id
my %id;
my %name;
my %rank;
my %parent;
foreach my $line(@nameIn){
	chomp($line);
	my @tmp = split(/\t/,$line);	# id  name  rank  parentID
	$id{$tmp[1]} = $tmp[0];
	$name{$tmp[0]} = $tmp[1];
	$rank{$tmp[0]} = $tmp[2];
	$parent{$tmp[0]} = $tmp[3];
}

#print "Parsing taxonNamesFull.txt done!!\n";

### GET LIST OF ALL TAXA from input file

######################## NEED TO DO SOMETHING HERE!!! #########################
#my $firstLine = $in[0];
#my @allTaxa = split(/\t/,$firstLine);
my @allTaxa = split(/\t/,$idList);
shift(@allTaxa);	# remove "geneID" tag

### getting ncbi ID of a complete taxonomy hierarchy for each NCBI Id from idList
my @defaultRanks = (
	'species','speciessubgroup','speciesgroup',
	'subgenus','genus',
	'subtribe','tribe',
	'subfamily','family','superfamily',
	'infraorder','parvorder','suborder','order','superorder',
	'infraclass','subclass','class','superclass',
	'subphylum','phylum','superphylum',
	'kingdom','superkingdom'
	);
my $allRanks = join("\t",@defaultRanks);

open(OUT,">$outDir/taxonID.list.fullRankID");
print OUT "No.\tabbrName\tncbiID\tfullName\t$allRanks\n";

open(NAMELIST,">$outDir/taxonNameReduced.txt");
print NAMELIST "ncbiID	fullName	rank	parentID\n";
my %reduceSpec;		# $reducedSpec{specID} = ncbiID	fullName	rank	parentID

my $c = 1;
#foreach my $line(@in){
foreach my $taxon(@allTaxa){
#	chomp($line);
	chomp($taxon);
#	if(length($line) > 0){
#		my @tmp = split(/\t/,$line);
#		my $ncbiID = $tmp[2];
	if(length($taxon) > 0){
		my $ncbiID = $taxon;
		$ncbiID =~ s/ncbi//;

		### list contains ID of all ranks
		my %rankID = (
			"superkingdom" => "NA", "kingdom" => "NA",
			"superphylum" => "NA", "phylum" => "NA", "subphylum" => "NA",
			"superclass" => "NA", "class" => "NA", "subclass" => "NA", "infraclass" => "NA",
			"superorder" => "NA", "order" => "NA", "suborder" => "NA", "parvorder" => "NA", "infraorder" => "NA",
			"superfamily" => "NA", "family" => "NA", "subfamily" => "NA",
			"tribe" => "NA", "subtribe" => "NA",
			"genus" => "NA", "subgenus" => "NA",
			"speciesgroup" => "NA", "speciessubgroup" => "NA", "species" => "NA"
			);

		### get name, parentID and parentRank of this ncbiID
		my $name = $name{$ncbiID};
		$rankID{$rank{$ncbiID}} = $ncbiID;

		my $parentID = $parent{$ncbiID};
		my $parentRank = $rank{$parentID};

#		print $ncbiID," - ",$name," - ",$rank{$ncbiID},"\n",$parentID," - ",$parentRank,"\n";

		### get all rankIDs
		my @allParentIDs = ($parentID);
		unless($parentID == 1){
			do{
				foreach my $rankName(@defaultRanks){
					#print $rankName;<>;
					if($parentRank eq $rankName){
						$rankID{$rankName} = $parentID;
					}
				}
				$parentID = $parent{$parentID};
				$parentRank = $rank{$parentID};
				push(@allParentIDs,$parentID);
#				print $parentID," - ",$parentRank," - ",$name{$parentID},"\n";
			} until ($parentID == 1);
		}

		### if any rankID not exists, get the last "norank" ID before hitting the next upper rank
		my $allParentIDs = join(";",@allParentIDs);
		for(my $i=0; $i<scalar(@defaultRanks)-1; $i++){
			if($rankID{$defaultRanks[$i]} eq "NA"){
#				print "no ID for $defaultRanks[$i]\n";
				my $tmpRank = $defaultRanks[$i];
				$tmpRank =~ s/sub//; $tmpRank =~ s/super//; $tmpRank =~ s/infra//; $tmpRank =~ s/parv//;
				if($defaultRanks[$i] =~ /species/){
					$tmpRank =~ s/subgroup//; $tmpRank =~ s/group//;
				}
				if($rankID{$tmpRank} ne "NA"){
					$rankID{$defaultRanks[$i]} = $rankID{$tmpRank};
				} else {
#					$rankID{$rankName} = getPrevID($allParentIDs,$defaultRanks[$i],$defaultRanks[$i+1]);
					$rankID{$defaultRanks[$i]} = getPrevID($allParentIDs,$defaultRanks[$i],$defaultRanks[$i+1]);
#					if($rankID{$defaultRanks[$i]} ne "NA"){ print "NEW ID = $rankID{$defaultRanks[$i]}\n";}
				}
			}
		}
		
		### taxonID for kingdom
		if($name{$rankID{"kingdom"}} =~ /group/){
			$rankID{"kingdom"} = $rankID{"superkingdom"};
		}		
		

		### print output
	#	print "$name\t$ncbiID";
#		print OUT "$line\t$name";
		print OUT "$c\t$taxon\t$ncbiID\t$name";
		for(my $i=0; $i<scalar(@defaultRanks); $i++){
	#		print "\t",$rankID{$rankName};
			my $currentID = $rankID{$defaultRanks[$i]};
			unless($currentID eq "NA"){
				print OUT "\t",$currentID;
				
				### get spec name info
				$reduceSpec{$currentID} = $currentID."\t".$name{$currentID}."\t".$rank{$currentID}."\t".$parent{$currentID};
			} else {
				my $replaceID = "";
				my $j = $i;
				do{
					$j ++;
					$replaceID = $rankID{$defaultRanks[$j]};
					
				} until ($replaceID ne "NA");
				print OUT "\t",$replaceID;
				
				### get spec name info
				$reduceSpec{$replaceID} = $replaceID."\t".$name{$replaceID}."\t".$rank{$replaceID}."\t".$parent{$replaceID};
			}
		}
#		<>;
		print OUT "\n";

#		print $c,"/",scalar(@in),"\n"; $c++;
	}
}

close (OUT);

foreach(keys %reduceSpec){
	print NAMELIST $reduceSpec{$_},"\n";
}	
close (NAMELIST);
#print "FINSHED!! Check output at\n\t$idList.fullRankID\n\t$idList.reducedName\n";

exit;


### if ID of a taxonomy rank not exist, go for the last "norank" ID before hitting the next upper rank
sub getPrevID{
	my ($allParents,$rank,$upperrank) = @_;

	my @allParents = split(/;/,$allParents);	# all taxonomy hierarchy IDs of the taxon in this rank
	my $out = "";

	for(my $i = 1; $i<scalar(@allParents); $i++){
#		print "### $allParents[$i]\n";<>;
		if($rank{$allParents[$i]} =~ /$upperrank/){
#			print "lower for $upperrank: ",$rank{$allParents[$i-1]},"\n";
			if($rank{$allParents[$i-1]} =~ /norank/){
				$out = $allParents[$i-1];
			}
		}
	}
	
	if(length($out) < 1){$out = "NA";}
	return $out;
}

sub getID{
	my $rankInfo = $_[0];
	if($rankInfo =~ /#/){
		my @tmp = split(/#/,$rankInfo);
		return ($tmp[2]);
	} else {
		return "NA";
	}
}