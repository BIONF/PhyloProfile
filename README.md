# PhyloProfile App

PhyloProfile App is a Shiny(R)-based app for visualizing the phylogenetic profile of a list of gene sequences.
The profile can be plotted with different taxonomy ranks (species, family, phylum, etc.), which would be useful for analyzing the presence absence of sequences in a large amount of taxa.
Currently, two information represented on the profile are the Feature Architecture Similarity scores and Percentage of species that have orthologs with the reference sequences. These two information can be filtered from the plot by changing the corresponding cutoff manually.

# Usage
Download all the files to your computer and keep their original folder structure.
Open Terminal, goto the folder that contains 2 R files ui.R and server.R
Run the app by typing: R -e 'shiny::runApp('',launch.browser=TRUE)'
Note: R has to be installed on your machine. The required packages will be automatically installed if necessary.
