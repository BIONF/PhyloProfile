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

#=for testing
#$idList = "geneID	ncbi329726	ncbi155978	ncbi245562	ncbi70779	ncbi7029	ncbi7159	ncbi380703	ncbi272557	ncbi176299	ncbi9646	ncbi5037	ncbi5039	ncbi29001	ncbi400682	ncbi272123	ncbi8839	ncbi46234	ncbi240292	ncbi1288291	ncbi28377	ncbi7165	ncbi278021	ncbi7460	ncbi218851	ncbi59689	ncbi3702	ncbi2234	ncbi696747	ncbi5105	ncbi33169	ncbi5057	ncbi36630	ncbi5059	ncbi746128	ncbi40384	ncbi5061	ncbi162425	ncbi5062	ncbi33178	ncbi1176130	ncbi44056	ncbi511995	ncbi5865	ncbi109871	ncbi374463	ncbi717646	ncbi264462	ncbi227086	ncbi7091	ncbi340100	ncbi9913	ncbi40559	ncbi55169	ncbi15368	ncbi7739	ncbi3711	ncbi135651	ncbi6238	ncbi6239	ncbi281687	ncbi31234	ncbi9483	ncbi76887	ncbi7868	ncbi1170562	ncbi99598	ncbi311458	ncbi360105	ncbi5476	ncbi42374	ncbi9615	ncbi5478	ncbi5480	ncbi5482	ncbi192875	ncbi81985	ncbi283909	ncbi387662	ncbi190650	ncbi10141	ncbi794803	ncbi46770	ncbi38033	ncbi1173020	ncbi211165	ncbi554065	ncbi331636	ncbi3055	ncbi184925	ncbi718219	ncbi9358	ncbi251229	ncbi7719	ncbi51511	ncbi85681	ncbi2711	ncbi5499	ncbi27339	ncbi36911	ncbi443906	ncbi212717	ncbi930089	ncbi5016	ncbi665024	ncbi5501	ncbi977863	ncbi930090	ncbi199306	ncbi665912	ncbi574566	ncbi930091	ncbi469383	ncbi5346	ncbi309798	ncbi360115	ncbi1173022	ncbi237895	ncbi40410	ncbi5807	ncbi5116	ncbi301224	ncbi1168544	ncbi3659	ncbi7176	ncbi755178	ncbi43989	ncbi292564	ncbi45157	ncbi43989	ncbi65393	ncbi395961	ncbi497965	ncbi41431	ncbi395962	ncbi292563	ncbi713887	ncbi56107	ncbi13035	ncbi7955	ncbi6669	ncbi9361	ncbi4959	ncbi159087	ncbi391774	ncbi693977	ncbi485916	ncbi44689	ncbi246195	ncbi5786	ncbi1150837	ncbi10020	ncbi112489	ncbi675120	ncbi548649	ncbi7217	ncbi7220	ncbi7222	ncbi7227	ncbi7230	ncbi7234	ncbi7237	ncbi7238	ncbi7240	ncbi7244	ncbi7260	ncbi7245	ncbi9371	ncbi2880	ncbi70536	ncbi269484	ncbi5802	ncbi280463	ncbi6035	ncbi27973	ncbi876142	ncbi31281	ncbi46681	ncbi5759	ncbi9796	ncbi9365	ncbi511145	ncbi71139	ncbi9685	ncbi372781	ncbi306281	ncbi98439	ncbi635003	ncbi31033	ncbi5518	ncbi59765	ncbi117187	ncbi8049	ncbi9031	ncbi69293	ncbi1173025	ncbi42742	ncbi243231	ncbi1173026	ncbi33072	ncbi251221	ncbi3847	ncbi9593	ncbi55529	ncbi862964	ncbi65093	ncbi309800	ncbi6412	ncbi464988	ncbi13563	ncbi9606	ncbi6087	ncbi100027	ncbi160233	ncbi6945	ncbi28985	ncbi381046	ncbi4914	ncbi374847	ncbi29883	ncbi30538	ncbi7897	ncbi5660	ncbi5671	ncbi347515	ncbi690899	ncbi5022	ncbi7918	ncbi741139	ncbi111781	ncbi83485	ncbi4006	ncbi7209	ncbi36914	ncbi372055	ncbi225164	ncbi9785	ncbi9315	ncbi9544	ncbi1126212	ncbi148305	ncbi3750	ncbi76773	ncbi3983	ncbi869210	ncbi394221	ncbi717774	ncbi3880	ncbi203908	ncbi100047	ncbi269797	ncbi259564	ncbi243233	ncbi265072	ncbi243232	ncbi2320	ncbi410358	ncbi43687	ncbi187420	ncbi449447	ncbi554155	ncbi392814	ncbi535722	ncbi30608	ncbi296587	ncbi1173027	ncbi4155	ncbi81824	ncbi13616	ncbi153609	ncbi36080	ncbi10090	ncbi9669	ncbi83344	ncbi1047171	ncbi59463	ncbi1168546	ncbi5762	ncbi160232	ncbi7425	ncbi348780	ncbi70791	ncbi242231	ncbi586133	ncbi45351	ncbi29176	ncbi1287680	ncbi5141	ncbi29879	ncbi40127	ncbi228410	ncbi497727	ncbi338192	ncbi323261	ncbi387092	ncbi61853	ncbi551115	ncbi40302	ncbi63737	ncbi317936	ncbi103690	ncbi28072	ncbi9978	ncbi9258	ncbi9986	ncbi8090	ncbi39947	ncbi56110	ncbi179408	ncbi242159	ncbi385169	ncbi70448	ncbi30611	ncbi9598	ncbi121759	ncbi5888	ncbi703506	ncbi121225	ncbi13735	ncbi5076	ncbi37727	ncbi31276	ncbi7757	ncbi5306	ncbi321614	ncbi2850	ncbi3885	ncbi4837	ncbi4787	ncbi3218	ncbi164328	ncbi67593	ncbi4929	ncbi4922	ncbi4924	ncbi263820	ncbi147573	ncbi5821	ncbi5825	ncbi5833	ncbi5849	ncbi5850	ncbi5854	ncbi5855	ncbi5861	ncbi5322	ncbi100044	ncbi118163	ncbi5145	ncbi705562	ncbi13642	ncbi9600	ncbi3694	ncbi104341	ncbi54126	ncbi9813	ncbi1219	ncbi146891	ncbi93059	ncbi93060	ncbi167546	ncbi59922	ncbi74546	ncbi74547	ncbi167542	ncbi167555	ncbi59920	ncbi167539	ncbi59919	ncbi3760	ncbi351746	ncbi82654	ncbi132908	ncbi56615	ncbi53953	ncbi70771	ncbi97479	ncbi426418	ncbi45151	ncbi10116	ncbi64495	ncbi37885	ncbi3988	ncbi373994	ncbi413404	ncbi4931	ncbi27288	ncbi4932	ncbi4934	ncbi114524	ncbi114525	ncbi27291	ncbi101203	ncbi9305	ncbi5334	ncbi4897	ncbi6183	ncbi4899	ncbi4896	ncbi653667	ncbi5180	ncbi128403	ncbi88036	ncbi692275	ncbi215467	ncbi85982	ncbi4555	ncbi671987	ncbi4081	ncbi42254	ncbi4558	ncbi448385	ncbi43179	ncbi109760	ncbi718229	ncbi40563	ncbi111780	ncbi13684	ncbi7668	ncbi641892	ncbi273057	ncbi387093	ncbi9823	ncbi56780	ncbi269084	ncbi1140	ncbi321327	ncbi321332	ncbi64471	ncbi110662	ncbi316279	ncbi195253	ncbi32049	ncbi1173263	ncbi316278	ncbi32051	ncbi84588	ncbi1148	ncbi1080228	ncbi1080229	ncbi1080230	ncbi59729	ncbi28564	ncbi9478	ncbi99883	ncbi5911	ncbi35128	ncbi273075	ncbi5874	ncbi146786	ncbi197221	ncbi98038	ncbi69014	ncbi5875	ncbi368408	ncbi5541	ncbi529818	ncbi292415	ncbi35720	ncbi5811	ncbi5217	ncbi10228	ncbi63577	ncbi7070	ncbi63418	ncbi203124	ncbi51453	ncbi29875	ncbi5691	ncbi364715	ncbi39416	ncbi37347	ncbi9739	ncbi33188	ncbi336963	ncbi5270	ncbi36033	ncbi103449	ncbi27335	ncbi27337	ncbi993615	ncbi29760	ncbi3068	ncbi6293	ncbi107463	ncbi8364	ncbi4952	ncbi1080233	ncbi4577	ncbi160035	ncbi4956	ncbi264203";
#$idList = "geneID	ncbi272557	ncbi76887	ncbi160233	ncbi43687	ncbi70771	ncbi273057	ncbi368408	ncbi2234	ncbi309800	ncbi269797	ncbi259564	ncbi243232	ncbi2320	ncbi410358	ncbi187420	ncbi348780	ncbi263820	ncbi53953	ncbi273075	ncbi69014	ncbi374847	ncbi160232	ncbi46770	ncbi497727	ncbi338192	ncbi443906	ncbi469383	ncbi176299	ncbi190650	ncbi269484	ncbi394221	ncbi264203	ncbi511995	ncbi641892	ncbi340100	ncbi718219	ncbi159087	ncbi265072	ncbi242231	ncbi228410	ncbi331636	ncbi329726	ncbi240292	ncbi43988	ncbi713887	ncbi33072	ncbi449447	ncbi551115	ncbi63737	ncbi1219	ncbi269084	ncbi321327	ncbi146786	ncbi203124	ncbi693977	ncbi869210	ncbi264462	ncbi391774	ncbi243231	ncbi448385	ncbi56780	ncbi360105	ncbi387092	ncbi387093	ncbi212717	ncbi309798	ncbi485916	ncbi380703	ncbi387662	ncbi360115	ncbi246195	ncbi511145	ncbi862964	ncbi717774	ncbi243233	ncbi323261	ncbi351746	ncbi413404	ncbi5833	ncbi237895	ncbi5807	ncbi5811	ncbi44689	ncbi5786	ncbi977863	ncbi5022	ncbi5346	ncbi40410	ncbi13563	ncbi29883	ncbi5306	ncbi104341	ncbi56615	ncbi5334	ncbi5217	ncbi5270	ncbi192875	ncbi3055	ncbi242159	ncbi109871	ncbi109760	ncbi55529	ncbi7897	ncbi9669	ncbi347515	ncbi5691	ncbi280463	ncbi5762	ncbi5741	ncbi7029	ncbi7159	ncbi9646	ncbi400682	ncbi28377	ncbi7165	ncbi7460	ncbi7091	ncbi9913	ncbi7739	ncbi6239	ncbi9483	ncbi9615	ncbi283909	ncbi10141	ncbi9031	ncbi9358	ncbi7719	ncbi51511	ncbi7176	ncbi7955	ncbi6669	ncbi9360	ncbi10020	ncbi7227	ncbi9371	ncbi9365	ncbi9685	ncbi8049	ncbi69293	ncbi9593	ncbi6412	ncbi9796	ncbi9606	ncbi6945	ncbi7918	ncbi7209	ncbi225164	ncbi9785	ncbi9315	ncbi9544	ncbi30608	ncbi13616	ncbi10090	ncbi59463	ncbi7425	ncbi45351	ncbi61853	ncbi9978	ncbi9258	ncbi8090	ncbi30611	ncbi9598	ncbi121225	ncbi13735	ncbi9823	ncbi9600	ncbi54126	ncbi9813	ncbi132908	ncbi9986	ncbi10116	ncbi9305	ncbi6183	ncbi42254	ncbi43179	ncbi7668	ncbi59729	ncbi31033	ncbi9478	ncbi99883	ncbi10228	ncbi7070	ncbi37347	ncbi9739	ncbi8364	ncbi6035	ncbi81824	ncbi36080	ncbi4837	ncbi64495	ncbi746128	ncbi5116	ncbi5518	ncbi148305	ncbi70791	ncbi5141	ncbi5076	ncbi5180	ncbi35720	ncbi39416	ncbi27337	ncbi45157	ncbi33170	ncbi5476	ncbi4959	ncbi28985	ncbi36913	ncbi4924	ncbi4952	ncbi4932	ncbi4787	ncbi3702	ncbi15368	ncbi3711	ncbi4577	ncbi3983	ncbi3880	ncbi39947	ncbi3218	ncbi3694	ncbi3760	ncbi88036	ncbi4555	ncbi4081	ncbi4558	ncbi3847	ncbi29760	ncbi4896	ncbi278021	ncbi40302	ncbi31281	ncbi136370";
#$idList = "geneID	ncbi272557	ncbi9360";	# example for input a genus ID instead of species or strain
#=cut

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

open(NAMELIST,">$outDir/taxonNamesReduced.txt");
print NAMELIST "ncbiID	fullName	rank	parentID\n";
my %reduceSpec;		# $reducedSpec{specID} = ncbiID	fullName	rank	parentID

my $c = 1;
#foreach my $line(@in){
foreach my $taxon(@allTaxa){
	chomp($taxon);
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

#			print $ncbiID," - ",$name," - ",$rank{$ncbiID},"\n",$parentID," - ",$parentRank,"\n";	############ TESTING

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
#					print $parentID," - ",$parentRank," - ",$name{$parentID},"\n";	############ TESTING
				} until ($parentID == 1);
			}

			### if input ID is not strain ID, force strain ID by this one
			if($rankID{"strain"} eq "NA"){
				$rankID{"strain"} = $ncbiID;
			}
			
			### if any rankID not exists, get the last identified ID before hitting this rank
			my $allParentIDs = join(";",@allParentIDs);
			for(my $i=1; $i<scalar(@defaultRanks)-1; $i++){
				if($rankID{$defaultRanks[$i]} eq "NA"){
#					print "no ID for $defaultRanks[$i]\n";	############ TESTING
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
						$rankID{$defaultRanks[$i]} = getPrevID($allParentIDs,$defaultRanks[$i],$defaultRanks[$i+1]);
					}
				}
			}
		
			### taxonID for kingdom
			if($rankID{"kingdom"} and $name{$rankID{"kingdom"}} and $name{$rankID{"kingdom"}} =~ /group/){
				$rankID{"kingdom"} = $rankID{"superkingdom"};
			}		

			### print output
#			print "$name\t$ncbiID";
			print OUT "$c\t$taxon\t$ncbiID\t$name";
			for(my $i=0; $i<scalar(@defaultRanks); $i++){
				my $currentID = $rankID{$defaultRanks[$i]};
#				print $defaultRanks[$i],"\t",$currentID;	############ TESTING
				
				unless($currentID eq "NA"){
					print OUT "\t",$currentID;
				
					### get spec name info
					$reduceSpec{$currentID} = $currentID."\t".$name{$currentID}."\t".$rank{$currentID}."\t".$parent{$currentID};
				} else {
					my $replaceID = "";
					my $j = $i;
					if($j == 0){$replaceID = $ncbiID;}
					else{
						do{
							$j --;
							$replaceID = $rankID{$defaultRanks[$j]};		
						} until ($replaceID ne "NA");
					}
#					print "\t($i - $j)replaced by ",$replaceID;	############ TESTING
					print OUT "\t",$replaceID;
				
					### get spec name info
					$reduceSpec{$replaceID} = $replaceID."\t".$name{$replaceID}."\t".$rank{$replaceID}."\t".$parent{$replaceID};
				}
#				<>;
			}

			if($rankID{"superkingdom"} == 2157){
				print OUT "\t2";
			} elsif($rankID{"superkingdom"} == 2759){
				print OUT "\t1";
			} elsif($rankID{"superkingdom"} == 2){
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

exit;


### if any rankID not exists, get the last identified ID before hitting this rank
sub getPrevID{
	my ($allParents,$rank,$upperrank) = @_;

	my @allParents = split(/;/,$allParents);	# all taxonomy hierarchy IDs of the taxon in this rank
	my $out = "";

	for(my $i = 1; $i<scalar(@allParents); $i++){
		if($rank{$allParents[$i]} =~ /$upperrank/){
#			print "lower for $upperrank: ",$rank{$allParents[$i-1]},"\n";	############ TESTING
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