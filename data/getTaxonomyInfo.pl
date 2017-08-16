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

#### get list of all species name and their NCBI taxon id
my %id;		# $id{NAME} = ID
my %name;	# $name{$id} = NAME
my %rank;	# $rank{$id} = RANK
my %parent;	# $parent{$id} = PARENT_ID

open(NAME,"$nameIN") || die "Cannot open $nameIN!!\n";
foreach my $line(<NAME>){
	chomp($line);
	my @tmp = split(/\t/,$line);	# id  name  rank  parentID
	$id{$tmp[1]} = $tmp[0];
	$name{$tmp[0]} = $tmp[1];
	$rank{$tmp[0]} = $tmp[2];
	$parent{$tmp[0]} = $tmp[3];
}
close (NAME);

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
print "Parsing taxonNamesFull.txt done!!\n";

### GET LIST OF ALL TAXA from input file

#=for testing
#$idList = "geneID	ncbi272557	ncbi9360";	# example for input a genus ID instead of species or strain
#$idList = "geneID	ncbi4837";#	ncbi4932	ncbi3702	ncbi9606";
#$idList = "geneID	ncbi436017";
#$idList = "geneID	ncbi1096996	ncbi2000001	ncbi2000002	ncbi2000003	ncbi202950	ncbi62977";
#$idList = "geneID	ncbi9606	ncbi9598	ncbi9544	ncbi10090	ncbi10116	ncbi9913	ncbi9612	ncbi9258	ncbi8364	ncbi31033	ncbi9031	ncbi7719	ncbi7739	ncbi6183	ncbi6239	ncbi7165	ncbi7227	ncbi6945	ncbi5270	ncbi5141	ncbi13616	ncbi7955	ncbi45351	ncbi5207";
#=cut

my @allTaxa = split(/\t/,$idList);
if($allTaxa[1] eq "geneID"){
	shift(@allTaxa);	# remove "geneID" tag
	if($allTaxa[1] eq "abbrName"){
		shift(@allTaxa);	# remove "abbrName" tag
	}
}

### get all ranking IDs for each taxon
my %info;	# $info{$ncbi} = rankID#rankName...
my %reduceSpec;		# $reducedSpec{specID} = ncbiID	fullName	rank	parentID
my %norankCount; # $norankCount{$norankID} = 2
my %rankCount; # count available IDs existing in each rank (except norank)

my $c = 1;

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

		$info{$ncbiID} = $name."#name";
		$rankCount{$rank} += 1;
		$rankCount{"name"} += 1;
		$rankCount{"strain"} += 1;

		unless($parentRank){
			print "parentID $parentID (of $ncbiID) not found\n";
		} else {
	#		print "HERE:\n$ncbiID\t$rank\n$parentID\t$parentRank\n";<>;	############ TESTING
			my @idArray = ($ncbiID,$parentID);
	#		print $name,"\t",$ncbiID,"\t","$parentID";
			$info{"$ncbiID"} .= "\t"."$ncbiID#strain\t"."$ncbiID#$rank\t"."$parentID#$parentRank";
	 		$reduceSpec{$ncbiID} = $ncbiID."\t".$name{$ncbiID}."\t".$rank."\t".$parentID;
	 		$reduceSpec{$parentID} = $parentID."\t".$name{$parentID}."\t".$parentRank."\t".$parent{$parentID};

			if($parentRank eq "norank"){
				unless($norankCount{$parentID}){
					$norankCount{$parentID} = 1;
				} else {
					$norankCount{$parentID} ++;
				}
			} else {
				unless($rankCount{$parentRank}){
					$rankCount{$parentRank} = 1;
				} else {
					$rankCount{$parentRank} += 1;
				}
			}


			unless($parentID == 1){
				do{
					$parentID = $parent{$parentID};
					$parentRank = $rank{$parentID};
	#				print $parentID,"\t",$parentRank,"\t",$name{$parentID},"\n";	############ TESTING
					push(@idArray,$parentID);

					$info{"$ncbiID"} .= "\t"."$parentID#$parentRank";

					if($parentRank eq "norank"){
						unless($norankCount{$parentID}){
							$norankCount{$parentID} = 1;
						} else {
							$norankCount{$parentID} ++;
						}
					} else {
						unless($rankCount{$parentRank}){
							$rankCount{$parentRank} = 1;
						} else {
							$rankCount{$parentRank} += 1;
						}
					}

					$reduceSpec{$parentID} = $parentID."\t".$name{$parentID}."\t".$parentRank."\t".$parent{$parentID};
				} until ($parentID == 1);
			}

	#		print "$info{$ncbiID}\n";
	#		print "NEXT...\n";		############ TESTING
	#		<>;						############ TESTING
		}
	}
#	print $c,"/",scalar(@allTaxa),"\n";
	$c++;
}

### print rank info
my %rankInfo;
my $mostInfo = "";	# get species with have the most taxonomy info
my %maxIndex; 	# used to save max index for each rank (and norank ID)
my %maxIndex2Rank;
my %allRankIndex = ("strain"=>1,"forma"=>2,"subspecies"=>3,"varietas"=>4,
										"species"=>5,"speciessubgroup"=>6,"speciesgroup"=>7,
										"subgenus"=>8,"genus"=>9,"subtribe"=>10,"tribe"=>11,
										"subfamily"=>12,"family"=>13,"superfamily"=>14,
										"parvorder"=>15,"infraorder"=>16,"suborder"=>17,"order"=>18,"superorder"=>19,
										"infraclass"=>20,"subclass"=>21,"class"=>22,"superclass"=>23,
										"subphylum"=>24,"phylum"=>25,"superphylum"=>26,
										"subkingdom"=>27,"kingdom"=>28,"superkingdom"=>29);
my %index2Rank = (1=>"strain",2=>"forma",3=>"subspecies",4=>"varietas",
									5=>"species",6=>"speciessubgroup",7=>"speciesgroup",
									8=>"subgenus",9=>"genus",10=>"subtribe",11=>"tribe",
									12=>"subfamily",13=>"family",14=>"superfamily",
									15=>"parvorder",16=>"infraorder",17=>"suborder",18=>"order",19=>"superorder",
									20=>"infraclass",21=>"subclass",22=>"class",23=>"superclass",
									24=>"subphylum",25=>"phylum",26=>"superphylum",
									27=>"subkingdom",28=>"kingdom",29=>"superkingdom");
#$c = 0;


### initial index for first rank "name"
$maxIndex{"name"} = 1;
$maxIndex2Rank{1} = "name";

my @allTaxaID = sort keys %info;

foreach my $taxID (sort keys %info){
	my @items = split(/\t/,$info{$taxID});
	my $c = 0;

	for(my $c = 1; $c < scalar(@items)+1; $c++){
		my @itemTMP = split(/#/,$items[$c-1]);		# $items[$c] = $taxID#$rank

		### check if this is a new rank or an existing rank in %maxIndex
		my $newRank = 1;
		if($itemTMP[1] eq "norank"){
			if($maxIndex{$itemTMP[0]}){
				$newRank = 0;
			}
		} else {
			if($maxIndex{$itemTMP[1]}){
				$newRank = 0;
			}
		}

		### for a new rank
		if($newRank == 1){
			if($itemTMP[1] eq "norank"){
				# check index of previous rank (of this species) to get index for this current norank ID
				my $p = $c-1-1;
				my $currentIndex = $c;
				while($p > -1){
					my @prevItemTMP = split(/#/,$items[$p]);
					if($prevItemTMP[1] eq "norank"){
						if($maxIndex{$prevItemTMP[0]}){
							$currentIndex = $maxIndex{$prevItemTMP[0]} + 1;
							last;
						} else {
							$p--;
						}
					} else {
						if($maxIndex{$prevItemTMP[1]}){
							$currentIndex = $maxIndex{$prevItemTMP[1]} + 1;
							last;
						} else {
							$p--;
						}
					}
				}

				# last index currently
				my $k = scalar(keys %maxIndex2Rank);
				# increase the index of higher ranks
				while($k >= $currentIndex){
					if($maxIndex2Rank{$k}){
						$maxIndex{$maxIndex2Rank{$k}} += 1;
						$maxIndex2Rank{$k+1} = $maxIndex2Rank{$k};
					}
					$k--;
				}
				# increase the current rank
				$maxIndex{$itemTMP[0]} = $currentIndex;
				$maxIndex2Rank{$currentIndex} = $itemTMP[0];
			} else {
				# if this is not the first main rank (name)
				if($itemTMP[1] ne "name"){
					# check if $c is the highest index currently
					my @sortedIndex = sort {$b<=>$a} keys %maxIndex2Rank;
					if($c > $sortedIndex[0]){
						$maxIndex{$itemTMP[1]} = $c;
						$maxIndex2Rank{$c} = $itemTMP[1];
					}
					# else, further check other ranks
					else {
						# check index of previous main rank
						my $j = $allRankIndex{$itemTMP[1]}-1;
						while($j > -1){
							# get index of previous rank if possible
							if($maxIndex{$index2Rank{$j}}){
								# last index currently
								my $k = scalar(keys %maxIndex2Rank);
								# increase the index of higher ranks
								while($k >= ($maxIndex{$index2Rank{$j}}+1)){
									if($maxIndex2Rank{$k}){
										$maxIndex{$maxIndex2Rank{$k}} += 1;
										$maxIndex2Rank{$k+1} = $maxIndex2Rank{$k};
									}
									$k--;
								}

								# then replace the next higher rank by the current rank
								$maxIndex{$itemTMP[1]} = $maxIndex{$index2Rank{$j}}+1;
								$maxIndex2Rank{$maxIndex{$index2Rank{$j}}+1} = $itemTMP[1];

								# stop the while loop
								last;
							} else {
								$j--;
							}
						}
					}
				}
			}
		}
		### for an existing rank
		else {
			if($itemTMP[1] eq "norank"){
				# if the old index smaller than $c (current position)
				if($maxIndex{$itemTMP[0]} < $c){
					# increase index for other higher ranks
					my $k = scalar(keys %maxIndex2Rank);
					while($k >= $c){
						if($maxIndex2Rank{$k}){
							$maxIndex{$maxIndex2Rank{$k}} += 1;
							$maxIndex2Rank{$k+1} = $maxIndex2Rank{$k};
						}
						$k--;
					}
					# increase the current rank
					$maxIndex{$itemTMP[0]} = $c;
					$maxIndex2Rank{$c} = $itemTMP[0];
				}
			} else {
				if($itemTMP[1] ne "name"){
					# if the old index smaller than $c (current position)
					if($maxIndex{$itemTMP[1]} < $c){
						# increase index for other higher ranks
						my $k = scalar(keys %maxIndex2Rank);
						while($k >= $c){
							if($maxIndex2Rank{$k}){
								$maxIndex{$maxIndex2Rank{$k}} += 1;
								$maxIndex2Rank{$k+1} = $maxIndex2Rank{$k};
							}
							$k--;
						}
						# increase the current rank
						$maxIndex{$itemTMP[1]} = $c;
						$maxIndex2Rank{$c} = $itemTMP[1];
					}
				}
			}
		}
	}
}


### print list of all sorted rank index
my $rankIndex;
foreach(sort {$a<=>$b} keys %maxIndex2Rank){
	$rankIndex .= "$maxIndex2Rank{$_};";
}
$rankIndex =~ s/;$//;
my @rankIndex = split(/;/,$rankIndex);



########## create taxonNamesReduced.txt
open(NAMELIST,">$outDir/taxonNamesReduced.txt");
print NAMELIST "ncbiID	fullName	rank	parentID\n";

foreach(keys %reduceSpec){
	print NAMELIST $reduceSpec{$_},"\n";
}
close (NAMELIST);


########## OUTPUT
open(OUT,">$outDir/taxonID.list.fullRankID");
my $rankIndexPrint = join("\t",@rankIndex);
$rankIndexPrint =~ s/\d+/norank/g;
$rankIndexPrint =~ s/name/fullName/;
print OUT "No.\tabbrName\tncbiID\t$rankIndexPrint\n";

$c = 1;
foreach my $taxID (sort keys %info){
	print OUT "$c\tncbi$taxID\t$taxID";

	my $items = "\t".$info{$taxID}."\t";

	# get ID for each rank in @rankIndex
	my $prevID = "";
	for(my $i = 0; $i < scalar(@rankIndex); $i++){
		my $id = "";
		if($rankIndex[$i] eq "name"){
			$items =~ /(.)+#name/;
			$id = $&; $id =~ s/#name//;
		} else{
			if($items =~ /\t(\d)+#$rankIndex[$i]?\t/){
				# get taxonomyID for this rank
				$id = $&;
				$id =~ s/#$rankIndex[$i]//;
			} elsif($items =~ /\t$rankIndex[$i]#/){
				# get taxonomyID for this "norank"
				$id = $rankIndex[$i];
			}
		}

		# print ID to output
		$id =~ s/\t//g;
		if(length($id)>0){
			#print "$id";<>;
			print OUT "\t$id";
			$prevID = $id;
		} else {
			#print $prevID;<>;
			print OUT "\t$prevID";
		}
	}
	$c++;
	print OUT "\n";
}

close (OUT);
print "FINISHED\n";
exit;
