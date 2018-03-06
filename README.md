# PhyloProfile
[![language: R](https://img.shields.io/badge/language-R-blue.svg?style=flat)](https://www.r-project.org/)
[![presented at: BOSC2017](https://img.shields.io/badge/presented%20at-BOSC2017-green.svg?style=flat)](https://f1000research.com/posters/6-1782)
[![GitHub release](https://img.shields.io/badge/stable%20release-v0.1.0-orange.svg)](https://github.com/BIONF/PhyloProfile/releases/tag/v0.2.1)
[![license: MIT](https://img.shields.io/badge/license-MIT-lightgrey.svg)](https://opensource.org/licenses/MIT)

[![](www/posterSub.png)](https://f1000research.com/posters/6-1782)
[Click for the full PDF version of the poster](https://f1000research.com/posters/6-1782)

*PhyloProfile* is a *Shiny*-based tool for integrating, visualizing and exploring multi-layered phylogenetic profiles.

Alongside the presence/absence patterns of orthologs across large taxon collections, *PhyloProfile* allows the integration of any two additional information layers. These complementary data, like sequence similarity between orthologs, similarities in their domain architecture, or differences in functional annotations enable a more informed interpretation of phylogenetic profiles.

By utilizing the NCBI taxonomy, *PhyloProfile* can dynamically collapse taxa into higher systematic groups. This enables rapidly changing the resolution from the comparative analyses of proteins in individual species to that of entire kingdoms or even domains without changes to the input data.

*PhyloProfile* furthermore allows for a dynamic filtering of profiles – taking the taxonomic distribution and the additional information layers into account. This, along with functions to estimate the age of genes and core gene sets facilitates the exploration and analysis of large phylogenetic profiles.

Take a look at [the functionality](https://github.com/BIONF/PhyloProfile/wiki/Functionality) of *PhyloProfile* and [explore the installation-free online version](https://phyloprofile.shinyapps.io/phyloprofile/) to learn more.

# Table of Contents
* [Installation &amp; Usage](#installation--usage)
* [Input Data](#input-data)
* [Walkthrough](#walkthrough)
* [Bugs](#bugs)
* [Acknowledgements](#acknowledgements)
* [Code of Conduct &amp; License](#code-of-conduct--license)
* [Contact](#contact)

# Installation & Usage
*PhyloProfile* is based on the *R*-package *Shiny*, as such a recent version of *R* is needed. Once that is out of the way you can just clone this repository to get a copy of *PhyloProfile*.

`git clone https://github.com/BIONF/phyloprofile`

To start PhyloProfile simply move into the PhyloProfile directory and run the main script

```
Rscript phyloprofile.R
```

The initial start can take a while, as `phyloprofile.R` will try do download and install all necessary dependencies automatically. *(Note: Depending on your system this sometimes fails, please check the console log for error messages concerning the dependency installation)*

Once all packages are downloaded and installed your web browser will open a new tab and display the main *PhyloProfile* menu.

# Input Data
*PhyloProfile* can read a number of different input files, including multi-FASTA files, regular tab-separated files and *OrthoXML*. The additional information layers can be embedded in the OrthoXML or be provided separately.

We described all suppported input formats in section [Input Data](https://github.com/BIONF/PhyloProfile/wiki/Input-Data) in our [PhyloProfile's Wiki](https://github.com/BIONF/PhyloProfile/wiki).

In `data/demo/` you can find some test data to see how the files should look like.

# Walkthrough
Read the [walkthrough slides](https://github.com/BIONF/PhyloProfile/wiki/Walkthrough) to explore the functionality of *PhyloProfile*.

# Bugs
Any bug reports or comments, suggestions are highly appreciated. Please open an issue on GitHub or be in touch via email.

# Acknowledgements
We would like to thank
1) [Bastian](https://github.com/gedankenstuecke) for the great initial idea and his kind support,
2) Members of [Ebersberger group](http://www.bio.uni-frankfurt.de/43045195/ak-ebersberger) for many valuable suggestions and ...bug reports :)

# Code of Conduct & License
This tool is released with a [Contributor Code of Conduct](https://github.com/BIONF/PhyloProfile/blob/master/CODE_OF_CONDUCT.md) & under [MIT license](https://github.com/BIONF/PhyloProfile/blob/master/LICENSE).

# Contact
Vinh Tran
tran@bio.uni-frankfurt.de
