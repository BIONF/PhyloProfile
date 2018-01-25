# Demo data
*PhyloProfile* comes with some test data that you can use to fully explore the functionality.

Detail about input files are expained in [Phyloprofile Wiki Page](https://github.com/BIONF/PhyloProfile/wiki/Input-Data).

## Main input files:
- `test.main.wide`: Input file wide (matrix) format. Use this as the **Main input** file on the *Input & settings* page after starting *PhyloProfile*. Each cell in the matrix contains 3 information: `Ortholog ID # Feature architecture similarity score[1] # Traceability score[2]`.
- `test.main.long`: This is the same as `test.main` but in long format.
- `test.main.xml`: This is the same input file in OrthoXML format.

## Optional input files:
- `domain_files/*.domains`: This folder contains the feature architecture data (e.g. Pfam domains) that you can optionally give under the **Additional annotation input** upload on the *Input & settings* page after startup.
- `fasta_files/*.fa`: This folder contains the fasta sequences for demo input. In order to show the sequence in **Detailed plot** you have to set path to fasta folder in **FASTA config** on the *Input & settings* page after startup.
- `other_optional_input/test.geneList`: After doing the initial plot with the files above you can use this file on the *Customized profile* tab to sub-select for only the genes present in this file.
- `test.taxaList`: This contains list of taxon names. Use this to test the function of fetching NCBI taxonomy IDs, which can be found in *Function* tab.

## Other pre-processing input files: for these files the user needs to run parsing scripts to prepare the compatible inputs for PhyloProfile.
- `oma_standalone_input/oma_example.orthoxml`: this is an output file from OMA standalone. Although it has OrthoXML format, but it lacks the taxonomy information. In order to convert it to an input file for PhyloProfile, you have to run this command `python scripts/convert_oma_standalone_orthoxml.py -x oma_example.orthoxml -m taxon_mapping_oma_orthoxml.csv > oma_example_phyloprofile_compatible.orthoxml`
- `fasta_input/test.input.fasta`: ortholog groups are saved in a multiple fasta file.
	- Main input: you can use `test.input.fasta` directly as main input file for PhyloProfile.
	- Do Pfam annotation for those sequences (you have to have PFAM-A HMM files ready):
		- Using hmmscan: `hmmscan -E 0.001 --noali --domtblout test.input.hmmscanOut path_to_Pfam-a_files/Pfam-A.hmm test.input.fasta`
		- Using pfamscan: `perl pfamscan.pl -fasta test.input.fasta -dir path_to_Pfam-A_files > test.input.pfamscanOut`
	- Convert hmmscan output to compatible domain file: `python scripts/hmmscanParser.py -i test.input.hmmscanOut > test.input.domains` (replace *hmmscanParser.py* by *pfamscanParser.py* for parsing pfamscan output file).
