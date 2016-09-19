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
    print "-i\tTaxonID list\n";
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

=for testing
#$idList = "geneID	ncbi33170	ncbi5476	ncbi42374	ncbi5478	ncbi5480	ncbi5482	ncbi36911	ncbi4959	ncbi28985	ncbi381046	ncbi4914	ncbi36913	ncbi4929	ncbi4922	ncbi4924	ncbi4931	ncbi27288	ncbi4932	ncbi4934	ncbi114524	ncbi114525	ncbi27291	ncbi36033	ncbi4952	ncbi4956	ncbi245562	ncbi1176130	ncbi717646	ncbi55169	ncbi794803	ncbi27339	ncbi930089	ncbi5016	ncbi665024	ncbi977863	ncbi930090	ncbi930091	ncbi1168544	ncbi112489	ncbi675120	ncbi548649	ncbi100027	ncbi5022	ncbi741139	ncbi372055	ncbi1126212	ncbi100047	ncbi1168546	ncbi1287680	ncbi147573	ncbi100044	ncbi97479	ncbi5027	ncbi45151	ncbi37885	ncbi692275	ncbi215467	ncbi5540	ncbi364715	ncbi1080233	ncbi160035	ncbi5499	ncbi665912	ncbi1150837	ncbi690899	ncbi703506	ncbi705562	ncbi671987	ncbi718229	ncbi107463	ncbi5037	ncbi5039	ncbi29001	ncbi5104	ncbi5057	ncbi36630	ncbi5059	ncbi746128	ncbi40384	ncbi162425	ncbi5061	ncbi5062	ncbi33178	ncbi40559	ncbi38033	ncbi5501	ncbi199306	ncbi5116	ncbi5518	ncbi59765	ncbi117187	ncbi148305	ncbi554155	ncbi535722	ncbi83344	ncbi1047171	ncbi70791	ncbi5141	ncbi29879	ncbi40127	ncbi121759	ncbi5076	ncbi37727	ncbi5145	ncbi5180	ncbi54788	ncbi5094	ncbi35720	ncbi63577	ncbi63418	ncbi51453	ncbi29875	ncbi39416	ncbi336963	ncbi27335	ncbi27337	ncbi321614	ncbi4897	ncbi4899	ncbi4896	ncbi5346	ncbi40410	ncbi42742	ncbi13563	ncbi29883	ncbi76773	ncbi203908	ncbi153609	ncbi5306	ncbi5322	ncbi104341	ncbi56615	ncbi5334	ncbi85982	ncbi40563	ncbi5217	ncbi5270	ncbi36080	ncbi4837	ncbi64495	ncbi109871	ncbi109760	ncbi27973	ncbi58839	ncbi6035	ncbi40302	ncbi31281	ncbi278021	ncbi70536	ncbi948595	ncbi586133	ncbi1288291	ncbi993615	ncbi10228	ncbi400682	ncbi70779	ncbi6087	ncbi45351	ncbi6183	ncbi7668	ncbi8839	ncbi7897	ncbi9669	ncbi13735	ncbi9646	ncbi28377	ncbi9913	ncbi7739	ncbi9483	ncbi9615	ncbi10141	ncbi9358	ncbi7719	ncbi51511	ncbi7955	ncbi9360	ncbi10020	ncbi9371	ncbi9796	ncbi9365	ncbi9685	ncbi31033	ncbi8049	ncbi9031	ncbi69293	ncbi9593	ncbi9606	ncbi30538	ncbi7918	ncbi9785	ncbi9315	ncbi9544	ncbi30608	ncbi13616	ncbi10090	ncbi59463	ncbi61853	ncbi9978	ncbi9258	ncbi9986	ncbi8090	ncbi30611	ncbi9598	ncbi7756	ncbi9600	ncbi9813	ncbi132908	ncbi10116	ncbi9305	ncbi42254	ncbi43179	ncbi9823	ncbi59729	ncbi9478	ncbi99883	ncbi37347	ncbi9739	ncbi8364	ncbi7868	ncbi225164	ncbi283909	ncbi6412	ncbi135651	ncbi6238	ncbi6239	ncbi281687	ncbi31234	ncbi6293	ncbi7209	ncbi54126	ncbi83485	ncbi7029	ncbi7159	ncbi7165	ncbi7460	ncbi7091	ncbi7176	ncbi6669	ncbi7217	ncbi7220	ncbi7222	ncbi7227	ncbi7230	ncbi7234	ncbi7237	ncbi7238	ncbi7240	ncbi7244	ncbi7260	ncbi7245	ncbi7425	ncbi6945	ncbi7070	ncbi121225	ncbi81824	ncbi192875	ncbi529818	ncbi227086	ncbi44689	ncbi5786	ncbi46681	ncbi5759	ncbi13642	ncbi5660	ncbi5671	ncbi347515	ncbi5691	ncbi5762	ncbi218851	ncbi59689	ncbi3702	ncbi15368	ncbi3711	ncbi81985	ncbi85681	ncbi2711	ncbi3659	ncbi71139	ncbi3847	ncbi4006	ncbi3750	ncbi3983	ncbi3880	ncbi4155	ncbi39947	ncbi3885	ncbi3218	ncbi3694	ncbi3760	ncbi3988	ncbi88036	ncbi4555	ncbi4081	ncbi4558	ncbi29760	ncbi4577	ncbi98038	ncbi554065	ncbi3055	ncbi392814	ncbi296587	ncbi242159	ncbi385169	ncbi70448	ncbi3068	ncbi574566	ncbi45157	ncbi44056	ncbi2880	ncbi186035	ncbi2850	ncbi4787	ncbi164328	ncbi67593	ncbi101203	ncbi5911	ncbi5888	ncbi31276	ncbi35128	ncbi5865	ncbi237895	ncbi5802	ncbi29176	ncbi5821	ncbi5825	ncbi5833	ncbi5849	ncbi5850	ncbi5854	ncbi5855	ncbi5861	ncbi5874	ncbi5875	ncbi5811	ncbi5807	ncbi280463	ncbi3027	ncbi55529	ncbi464988	ncbi2234	ncbi259564	ncbi2320	ncbi410358	ncbi348780	ncbi309800	ncbi269797	ncbi243232	ncbi187420	ncbi263820	ncbi53953	ncbi273075	ncbi69014	ncbi160232	ncbi374847	ncbi272557	ncbi160233	ncbi43687	ncbi70771	ncbi368408	ncbi76887	ncbi273057	ncbi311458	ncbi46770	ncbi338192	ncbi497727	ncbi693977	ncbi869210	ncbi212717	ncbi309798	ncbi485916	ncbi155978	ncbi329726	ncbi1165	ncbi1167	ncbi240292	ncbi118562	ncbi1187	ncbi1173032	ncbi1124	ncbi28069	ncbi54299	ncbi241425	ncbi379064	ncbi43988	ncbi59930	ncbi102235	ncbi713887	ncbi142864	ncbi292566	ncbi1191	ncbi1488322	ncbi33072	ncbi76023	ncbi47254	ncbi449447	ncbi44472	ncbi551115	ncbi63737	ncbi118323	ncbi482564	ncbi44475	ncbi1219	ncbi1153	ncbi34078	ncbi102116	ncbi269084	ncbi1140	ncbi321327	ncbi321332	ncbi146786	ncbi203124	ncbi443906	ncbi469383	ncbi331636	ncbi718219	ncbi511995	ncbi641892	ncbi360105	ncbi387092	ncbi387093	ncbi264462	ncbi391774	ncbi243231	ncbi448385	ncbi56780	ncbi176299	ncbi190650	ncbi269484	ncbi394221	ncbi264203	ncbi340100	ncbi159087	ncbi265072	ncbi242231	ncbi228410	ncbi292415	ncbi380703	ncbi374463	ncbi387662	ncbi360115	ncbi246195	ncbi511145	ncbi862964	ncbi717774	ncbi243233	ncbi323261	ncbi351746	ncbi413404";
#$idList = "geneID	ncbi5762";
=cut

my @allTaxa = split(/\t/,$idList);
shift(@allTaxa);	# remove "geneID" tag

### getting ncbi ID of a complete taxonomy hierarchy for each NCBI Id from idList
my @defaultRanks = (
	'strain',
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
print OUT "No.\tabbrName\tncbiID\tfullName\t$allRanks\tgroup\n";

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
			"speciesgroup" => "NA", "speciessubgroup" => "NA", "species" => "NA",
			"strain" => "NA"
			);

		### get name, parentID and parentRank of this ncbiID
		unless($name{$ncbiID}){
			print "CHECK $ncbiID\n";
		} else {
			my $name = $name{$ncbiID};
			if($rank{$ncbiID} eq "norank"){
				$rankID{"strain"} = $ncbiID;
			} else {
				$rankID{$rank{$ncbiID}} = $ncbiID;
			}

			my $parentID = $parent{$ncbiID};
			my $parentRank = $rank{$parentID};

#			print $ncbiID," - ",$name," - ",$rank{$ncbiID},"\n",$parentID," - ",$parentRank,"\n";

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
#					print $parentID," - ",$parentRank," - ",$name{$parentID},"\n";
				} until ($parentID == 1);
			}

			### if any rankID not exists, get the last identified ID before hitting this rank
			my $allParentIDs = join(";",@allParentIDs);
			for(my $i=0; $i<scalar(@defaultRanks)-1; $i++){
				if($rankID{$defaultRanks[$i]} eq "NA"){
#					print "no ID for $defaultRanks[$i]\n";<>;
					my $tmpRank = $defaultRanks[$i];
					
					### if current rank has prefix of "sub", "super", "infra" or "parv" (or "subgroup", "group" in case of "species" rank)
					### than remove that prefix (to check with "main" rank)
					$tmpRank =~ s/sub//; $tmpRank =~ s/super//; $tmpRank =~ s/infra//; $tmpRank =~ s/parv//;
					if($defaultRanks[$i] =~ /species/){
						$tmpRank =~ s/subgroup//; $tmpRank =~ s/group//;
					}
					
					### if "main" rank is defined, transfer it's ID to subrank
					if($rankID{$tmpRank} ne "NA"){
						$rankID{$defaultRanks[$i]} = $rankID{$tmpRank};
					} 
					### else, get the last identified ID before hitting this main rank
					else {
#						$rankID{$rankName} = getPrevID($allParentIDs,$defaultRanks[$i],$defaultRanks[$i+1]);
						$rankID{$defaultRanks[$i]} = getPrevID($allParentIDs,$defaultRanks[$i],$defaultRanks[$i+1]);
#						if($rankID{$defaultRanks[$i]} ne "NA"){ print "NEW ID = $rankID{$defaultRanks[$i]}\n";}
					}
				}
			}
		
			### taxonID for kingdom
			if($rankID{"kingdom"} and $name{$rankID{"kingdom"}} and $name{$rankID{"kingdom"}} =~ /group/){
				$rankID{"kingdom"} = $rankID{"superkingdom"};
			}		

			### print output
#			print "$name\t$ncbiID";
#			print OUT "$line\t$name";
			print OUT "$c\t$taxon\t$ncbiID\t$name";
			for(my $i=0; $i<scalar(@defaultRanks); $i++){
				my $currentID = $rankID{$defaultRanks[$i]};
#				print $defaultRanks[$i],"\t",$currentID;
				
				unless($currentID eq "NA"){
					print OUT "\t",$currentID;
				
					### get spec name info
					$reduceSpec{$currentID} = $currentID."\t".$name{$currentID}."\t".$rank{$currentID}."\t".$parent{$currentID};
				} else {
					my $replaceID = "";
					my $j = $i;
					do{
						$j --;
						$replaceID = $rankID{$defaultRanks[$j]};
						
					} until ($replaceID ne "NA");
#					print "\treplaced by ",$replaceID;
					print OUT "\t",$replaceID;
				
					### get spec name info
					$reduceSpec{$replaceID} = $replaceID."\t".$name{$replaceID}."\t".$rank{$replaceID}."\t".$parent{$replaceID};
				}
#				<>;
			}

			if($rankID{"superkingdom"} == 2157){
#				print "archaea 2";
				print OUT "\t2";
			} elsif($rankID{"superkingdom"} == 2759){
#				print "eukaryota 1";
				print OUT "\t1";
			} elsif($rankID{"superkingdom"} == 2){
#				print "bacteria 3";
				print OUT "\t3";
			}
#			<>;

			print OUT "\n";

#			print $c,"/",scalar(@in),"\n"; 
			$c++;
		}
	}
}

close (OUT);

foreach(keys %reduceSpec){
	print NAMELIST $reduceSpec{$_},"\n";
}	
close (NAMELIST);
#print "FINSHED!! Check output at\n\t$idList.fullRankID\n\t$idList.reducedName\n";

exit;


### if any rankID not exists, get the last identified ID before hitting this rank
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