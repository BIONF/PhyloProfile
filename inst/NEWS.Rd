\name{PhyloProfilenews}
\title{PhyloProfile news}
\encoding{UTF-8}

\section{Version 1.6.1}{
   \itemize{
     \item modified taxonomy ranks (ropensci#875)
     \item improved rank indexing function (ropensci#874)
     \item improved x-axis label (#116)
   }
 }

\section{Version 1.4.12}{
  \itemize{
    \item domain plot works with group ID containing pipe
  }
}

\section{Version 1.4.11}{
  \itemize{
    \item added check for invalid ids in newTaxa file
  }
}

\section{Version 1.4.10}{
  \itemize{
    \item fixed bug comparing unk taxon IDs with max ncbi taxID
  }
}

\section{Version 1.4.9}{
  \itemize{
    \item fixed bug comparing taxon IDs with max ncbi taxID
  }
}

\section{Version 1.4.8}{
  \itemize{
    \item fixed bug filterProfileData (some filtered data could not be plotted)
  }
}

\section{Version 1.4.5}{
  \itemize{
    \item Add auto sizing feature
  }
}

\section{Version 1.4.4}{
  \itemize{
    \item Fixed orthoxml parsing function
  }
}

\section{Version 1.4.3}{
  \itemize{
    \item Plot background for NA orthologs
    \item Fixed highlight by categories
  }
}

\section{Version 1.4.2}{
  \itemize{
    \item Fixed error when parsing large amount of new taxa
    \item Improved the speed of parsing new taxa
    \item Apply filter to refspec orthologs when not working with strain/species
  }
}
 
\section{Version 1.4.1}{
  \itemize{
    \item Fixed bug data lost after filtering var2
    \item Add filter data to detailed plot
  }
}

\section{Version 1.3.11}{
  \itemize{
    \item Add filter number of co-orthologs to parseInfoProfile()
  }
}

\section{Version 1.3.10}{
  \itemize{
    \item Calculate percentage of present taxa after filtering of var1 and var2
    \item Fixed filter when working with high taxonomy ranks
  }
}

\section{Version 1.2.8}{
  \itemize{
    \item Added new NCBI taxonomy ranks (e.g. biotype, isolate, pathogroup, ...)
    \item Added function to reset taxonomy data
  }
}

\section{Version 1.2.7}{
  \itemize{
    \item Solved problem with new NCBI taxonomy rank "clade" by replace them with "norank"
  }
}

\section{Version 1.2.6}{
  \itemize{
    \item Fixed bug customized profile of subset of taxa not clickable
  }
}

\section{Version 1.2.4}{
  \itemize{
    \item Fixed bug checking invalid taxon IDs
  }
}

\section{Version 1.2.1}{
  \itemize{
    \item Fixed bug in rankIndexing and processNcbiTaxonomy
    \item Improved check for invalid input taxon IDs
  }
}

\section{Version 1.1.12}{
  \itemize{
    \item Added pseudo IDs for Holozoa and Holomycota clades and their non-NCBI sub-clades
  }
}

\section{Version 1.1.11}{
  \itemize{
    \item Increased the resolution of the phylostratigraphy function
  }
}

\section{Version 1.1.10}{
  \itemize{
    \item Improved the rank indexing function
  }
}

\section{Version 1.1.6}{
  \itemize{
    \item Removed grDevices, grid, stats, utils from Imports list in DESCRIPTION
  }
}

\section{Version 1.1.5}{
  \itemize{
    \item Added progress indicators for data loading and processing steps
    \item Removed shinycssloaders
  }
}

\section{Version 1.1.3}{
  \itemize{
    \item Fixed bug listing taxa by selecting a supertaxon from a higher level
    \item Added message for downloading preprocessed data
  }
}

\section{Version 1.1.2}{
  \itemize{
    \item Fixed bug sorting taxa when taxonomy rank is set as strain
  }
}

\section{Version 1.1.1}{
  \itemize{
    \item Commented examples of OMA parser functions
    \item Turned off tests for OMA parser functions
  }
}

\section{Version 0.99.31}{
  \itemize{
    \item Submitted to Bioconductor
  }
}