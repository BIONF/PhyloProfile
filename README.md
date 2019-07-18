# PhyloProfile
[![GitHub release](https://img.shields.io/badge/latest%20release-1.0.0-orange.svg)](https://github.com/BIONF/PhyloProfile/releases/tag/v1.0.0)
[![published in: Bioinformatics](https://img.shields.io/badge/published%20in-Bioinformatics-ff69b4.svg?style=flat)](https://doi.org/10.1093/bioinformatics/bty225)
[![poster at: BOSC2017](https://img.shields.io/badge/poster%20at-BOSC2017-green.svg?style=flat)](https://f1000research.com/posters/6-1782)
[![presented at: GCB2018](https://img.shields.io/badge/presented%20at-GCB2018-green.svg?style=flat)](http://gcb2018.de)
[![language: R](https://img.shields.io/badge/language-R-blue.svg?style=flat)](https://www.r-project.org/)
[![license: MIT](https://img.shields.io/badge/license-MIT-lightgrey.svg)](https://opensource.org/licenses/MIT)
[![Travis-CI Build Status](https://travis-ci.org/BIONF/PhyloProfile.svg?branch=master)](https://travis-ci.org/BIONF/PhyloProfile)

![](https://github.com/BIONF/phyloprofile-data/blob/gh-pages/www/posterSub.png)
[Click here for the full PDF version of the BOSC2017 poster](https://f1000research.com/posters/6-1782)

*PhyloProfile* is a *Shiny*-based tool for integrating, visualizing and exploring multi-layered [phylogenetic profiles](https://en.wikipedia.org/wiki/Phylogenetic_profiling).

Alongside the presence/absence patterns of [orthologs](https://en.wikipedia.org/wiki/Homology_(biology)) across large [taxon](https://en.wikipedia.org/wiki/Taxon) collections, *PhyloProfile* allows the integration of any two additional information layers. These complementary data, like [sequence similarity](https://en.wikipedia.org/wiki/Sequence_alignment) between orthologs, similarities in their [domain architecture](https://www.ncbi.nlm.nih.gov/pubmed/20221914), or differences in [functional annotations](https://en.wikipedia.org/wiki/Protein_function_prediction) enable a more informed interpretation of phylogenetic profiles.

By utilizing the [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy), *PhyloProfile* can dynamically collapse taxa into higher [systematic groups](https://en.wikipedia.org/wiki/Taxonomy_(biology)). This enables rapidly changing the resolution from the comparative analyses of proteins in individual species to that of entire kingdoms or even domains without changes to the input data.

*PhyloProfile* furthermore allows for a dynamic filtering of profiles – taking the taxonomic distribution and the additional information layers into account. This, along with functions to estimate [the age of genes](http://www.cell.com/trends/genetics/fulltext/S0168-9525(13)00111-X) and [core gene](https://en.wikipedia.org/wiki/Pan-genome) sets facilitates the exploration and analysis of large phylogenetic profiles.

Take a look at [the functionality](https://github.com/BIONF/PhyloProfile/wiki/Functionality) of *PhyloProfile* and [explore the installation-free online version](http://applbio.biologie.uni-frankfurt.de/phyloprofile/) to learn more.

# Table of Contents
* [Installation &amp; Usage](#installation--usage)
* [Input Data](#input-data)
* [Walkthrough](#walkthrough)
* [Bugs](#bugs)
* [Acknowledgements](#acknowledgements)
* [Code of Conduct &amp; License](#code-of-conduct--license)
* [Contact](#contact)

# Installation & Usage
## Using BiocManager
We are preparing to submit *PhyloProfile* to [Bioconductor](https://www.bioconductor.org/). Once the package is accepted, you can easily install it using BiocManager:

```r
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("PhyloProfile")
```

*Until then please install the dev version of *PhyloProfile* from our github repository using __[devtools](https://cran.r-project.org/web/packages/devtools/index.html)__.*
## Using devtools
The dev version of *PhyloProfile* can be installed from this github repository using `devtools`:

```r
if (!requireNamespace("devtools"))
    install.packages("devtools")
devtools::install_github("BIONF/PhyloProfile", INSTALL_opts = c('--no-lock'), build_vignettes = TRUE)
```

This step can take a while, as the tool will try do download and install all necessary dependencies automatically. *(Note: Depending on your system this sometimes fails, please check the console log for error messages concerning the dependency installation)*

## Start PhyloProfile
First install the demo data package
```r
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install(version = '3.10', ask = FALSE)
devtools::install_github(
    "bionf/PhyloProfileData",
    INSTALL_opts = c('--no-lock'),
    build_opts = c('--no-resave-data')
)
```

Then run
```r
library(PhyloProfile)
runPhyloProfile()
```
Check your web browser, *PhyloProfile* will be displayed there ;-) For the first time running, the tool will download a [pre-caculated taxonomy data](https://github.com/BIONF/phyloprofile-data). Please be patient until you see a message for uploading input files.

_**Please check our [detailed instructions](https://github.com/BIONF/PhyloProfile/wiki/Installation) if you encounter any problems while installing and starting the program.**_

# Input Data
*PhyloProfile* can read a number of different input files, including multi-FASTA files, regular tab-separated files, OMA ID list or *OrthoXML*. The additional information layers can be embedded in the OrthoXML or be provided separately.

We described all suppported input formats in section [Input Data](https://github.com/BIONF/PhyloProfile/wiki/Input-Data) in our [PhyloProfile's Wiki](https://github.com/BIONF/PhyloProfile/wiki).

# Walkthrough and Examples
Read the [walkthrough slides](https://github.com/BIONF/PhyloProfile/wiki/Walkthrough) to explore the functionality of the *PhyloProfile* GUI.

Check the vignette for learning how to use PhyloProfile's functions in some specific use-cases:
```r
browseVignettes("PhyloProfile")
```

# Bugs
Any [bug reports or comments, suggestions](https://github.com/BIONF/PhyloProfile/blob/master/CONTRIBUTING.md) are highly appreciated. Please [open an issue on GitHub](https://github.com/BIONF/PhyloProfile/issues/new) or be in touch via email.

# Acknowledgements
We would like to thank
1) [Bastian](https://github.com/gedankenstuecke) for the great initial idea and his kind support,
2) Members of [Ebersberger group](http://www.bio.uni-frankfurt.de/43045195/ak-ebersberger) for many valuable suggestions and ...bug reports :)

# Contributors
* [Vinh](https://github.com/trvinh)
* [Bastian](https://github.com/gedankenstuecke)
* [Carla](https://github.com/CarlaMoelbert)

# License
This tool is released under [MIT license](https://github.com/BIONF/PhyloProfile/blob/master/LICENSE).

# How-To Cite
Ngoc-Vinh Tran, Bastian Greshake Tzovaras, Ingo Ebersberger; PhyloProfile: Dynamic visualization and exploration of multi-layered phylogenetic profiles, Bioinformatics, , bty225, https://doi.org/10.1093/bioinformatics/bty225

or use the citation function in R CMD to have it directly in BibTex or LaTeX format
```r
citation("PhyloProfile")
```
# Contact
Vinh Tran
tran@bio.uni-frankfurt.de
