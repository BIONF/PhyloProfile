#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;

=desc
(2016.11.01)
Parsing domain architecture of proteins from FAS XML output file

**Example input:
<?xml version="1.0"?>
<out direction="single-->set" weighting="applied">
	<single_protein id="nempa_5256_1:NEPG_00013" length="933">
		<feature type="seg_low complexity regions" evalue="NULL" weight="1.0">
			<instance inst_eval="NULL" start="323" end="337"/>
			<instance inst_eval="NULL" start="893" end="904"/>
		</feature>
	</single_protein>
	<set_protein id="mja:MJ1072" score="0.647488634" MS="0.5" PS="0.99162878" CS="0.0" length="116" mode="greedy">
		<architecture>
			<feature type="coils_coiled_coil" evalue="NULL">
				<instance inst_eval="NULL" start="52" end="72"/>
			</feature>
			<feature type="seg_low complexity regions" evalue="NULL">
				<instance inst_eval="NULL" start="36" end="48"/>
			</feature>
			<feature type="tmhmm_transmembrane" evalue="NULL">
				<instance inst_eval="NULL" start="87" end="109"/>
			</feature>
		</architecture>
		<path>
			<feature type="seg_low complexity regions" corrected_weight="0.5">
				<instance start="323" end="337"/>
				<instance start="893" end="904"/>
			</feature>
		</path>
	</set_protein>
</out>

**Output contains tabs "seedID"  "orthoID"  "feature"  "start"  "end"  "weight":

OG_1000#mja:MJ1072#nempa_5256_1:NEPG_00013	mja:MJ1072	coils_coiled_coil	52	72	NA
OG_1000#mja:MJ1072#nempa_5256_1:NEPG_00013	mja:MJ1072	seg_low complexity regions	36	48	NA
OG_1000#mja:MJ1072#nempa_5256_1:NEPG_00013	mja:MJ1072	tmhmm_transmembrane	87	109	NA
OG_1000#mja:MJ1072#nempa_5256_1:NEPG_00013	nempa_5256_1:NEPG_00013	seg_low complexity regions	323	337	0.5
OG_1000#mja:MJ1072#nempa_5256_1:NEPG_00013	nempa_5256_1:NEPG_00013	seg_low complexity regions	893	904	0.5

=cut

sub usage {
    my $msg = shift;
    print "example: perl parseArchitecture.pl -i xmlOutput.fas -o output.mDomains -g groupID\n";
    print "-i\tInput of XML file from FAS calculation\n";
    print "-o\tOutput file\n";
    print "-g\tGroup ID\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_o,$opt_g);
getopts('i:o:g:');

# sanity checks;
my $inFile = ($opt_i) ? $opt_i : usage("ERROR: No input file given\n");
my $outFile = ($opt_o) ? $opt_o : usage("ERROR: No output file given\n");
my $groupID = ($opt_g) ? $opt_g : usage("ERROR: No group ID given\n");

open(OUT,">$outFile") || die "Cannot create $outFile!\n";
open(IN,$inFile) || die "Cannot open $inFile!\n";

### MAIN
my @in = <IN>;
close (IN);
my $in = join("",@in);

### split multiple XML file into individual files if necessary
my @xml = split(/<\/out>/,$in);

foreach my $archi(@xml){
	if(length($archi) > 10){
		my @archiTMP = split(/<\/single_protein>/,$archi);

		### get search species ID (if available)
		my $searchSpec = "";
		if($archi =~ /(.)+?\.xml/){
			my $hit = $&;		# plaga_4069@5849@1_22941_fas.xml
			my @hit = split(/_/,$hit); 
			pop(@hit), pop(@hit);
			$searchSpec = join("_",@hit);
		}
#		if($archi =~ /niteu_5654/){
#			print "YES!! $searchSpec\n";<>;
#		}

		### get seed ID
		my $seedID = "";
		if($archiTMP[0] =~ /single_protein id=\"(.)+?\"/){
			$seedID = $&;
			$seedID =~ s/single_protein id=//; $seedID =~ s/\"//g;
		}

		### go through all search proteins in this xml file
		my @searchProts = split(/<\/set_protein>/,$archiTMP[1]);

		foreach my $searchProt(@searchProts){
			### get search protein ID(s)
			my $searchID = "";
			if($searchProt =~ /set_protein id=\"(.)+?\"/){
				$searchID = $&;
				$searchID =~ s/set_protein id=//; $searchID =~ s/\"//g;
				if(length($searchSpec) > 0){ $searchID = $groupID."|".$searchSpec."|".$searchID;}

				### get info
				my @info = split(/<\/architecture>/,$searchProt);

				### get protein's domain positions
				my $searchDomain = getDomainPos($groupID,$searchID,$seedID,$info[0],1);
				print OUT $searchDomain;

				### get seed's best path
				my $seedDomain = getDomainPos($groupID,$searchID,$seedID,$info[1],0);
				print OUT $seedDomain;
			}
		}
	}
}
close (OUT);

print "Finished! Check output at $outFile\n";
exit;

sub getDomainPos{
	my ($groupID,$searchID,$seedID,$block,$order) = @_;

	my @features = split(/feature/,$block);
	my $result = "";

	foreach my $feature (@features){
		if($feature =~ /start/){
			my @info = split(/\n/,$feature);

			my $firstLine = shift(@info);
			my $type = "";
			if($firstLine =~ /type=\"(.)+?\"/){
				$type = $&;
				$type =~ s/type=//; $type =~ s/\"//g;
			}

			my $weight = "NA";
			if($firstLine =~ /weight=\"(.)+?\"/){
				$weight = $&; $weight =~ s/weight=//; $weight =~ s/\"//g;
			}

			foreach my $infoLine(@info){
				chomp($infoLine);
				if($infoLine =~ /start/){
#					print $line,"\n";
					my $start = ""; my $end = "";
					if($infoLine =~ /start=\"\d+\"/){
						$start = $&; $start =~ s/start=//; $start =~ s/\"//g;
					}
					if($infoLine =~ /end=\"\d+\"/){
						$end = $&; $end =~ s/end=//; $end =~ s/\"//g;
					}

					if($order == 1){
						$result .= "$groupID#$searchID#$seedID\t$searchID\t$type\t$start\t$end\t$weight\n";
					} elsif($order == 0){
						$result .= "$groupID#$searchID#$seedID\t$seedID\t$type\t$start\t$end\t$weight\n";
					}
				}
			}
		}
	}
	return $result;
}
