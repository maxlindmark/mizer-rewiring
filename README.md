# Warming alters the effect of fishing on the size spectra of an exploited temperate food web
## Description
Here we parameterize and calibrate a multispecies size-spectrum model (mizer) for the Baltic Sea. We use this as a case study to assess interactive effects of fishing and climate change on size structure and yield of a temperate food web.

**Collaborators**: Asta Audzijonyte, Julia Blanchard and Anna Gårdmark.
We use non-released version of mizer with added functionality (multiple resource background spectra and temperature-dependence of vital rates), developed also with Jon Reum and Romain Forestier. This is a fork of "astaaudzi/mizer-rewiring, ref = "rewire-temp". I have kept all files and put Baltic files in the Baltic folder.


## Mizer
mizer web page: https://sizespectrum.github.io/mizer/dev/

**Developers**

Finlay Scott, Author, copyright holder

Julia Blanchard, Author, copyright holder  

Ken Andersen, Author, copyright holder

Gustav Delius, Contributor, maintainer

**Citation for mizer** 
Scott F, Blanchard JL, Andersen KH (2014). "mizer: an R package for multispecies, trait-based and community size spectrum ecological modelling." *Methods Ecol Evol*, **5**, 1121-1125. doi: 10.1111/2041-210X.12256.


## Contents
The **baltic/data** folder contains data for estimating parameters and for compiling time series for calibration. Actually currently too big! Will need to figure this out, currently only existing locally.

The **baltic/params** folder contains .csv files with parameters that are read in modelling scripts, see **R/Parameter_estimation**

The **Baltic/figures** folder contains figures for the publication

The **baltic/R/parameter_estimation** folder contains code for estimating parameters and for compiling time series for calibration. 

The **baltic/R/MSSM** folder contains R code for model calibration code and simulations.

The **baltic/R/functions** folder contains functions for calibtration and loading packages.


## Data sources (ICES) 
### von Bertalanffy growth parameters & Weight-Length relationships (individual length-weight-age)
ICES Database of Trawl Surveys (DATRAS), Extraction 1 JAN 2018 of International Bottom Trawl Survey (BITS). ICES, Copenhagen, available at http://www.ices.dk/marine-data/data-portals/Pages/DATRAS.aspx


### Times series of Abundance and F
**Cod**:	ICES. 2013. *Report of the Baltic Fisheries Assessment Working Group (WGBFAS)*, 10 - 17 April 2013, ICES Headquarters, Copenhagen. ICES CM 2013/ACOM:10. 747 pp. [Final accepted analytical assessment by ICES; done using SAM model]

**Sprat**:	ICES. 2015. *Report of the Baltic Fisheries Assessment Working Group (WGBFAS)*, 14-21 April 2015, ICES HQ, Copenhagen, Denmark. ICES CM 2015/ACOM:10. [from XSA model] 

**Herring**: ICES. 2015. *Report of the Baltic Fisheries Assessment Working Group (WGBFAS)*, 14-21 April 2015, ICES HQ, Copenhagen, Denmark. ICES CM 2015/ACOM:10. [from XSA model]

*Eventually have a look at this guide to good "Read me's": https://gist.github.com/PurpleBooth/109311bb0361f32d87a2 and this example: https://github.com/seananderson/heavy-tails*




# mizer

mizer is a package that implements size-based ecological models.
The package has been developed to model marine ecosystems that are subject
to fishing. However, it may also be appropriate for other ecosystems.

The package contains routines and methods to allow users to set up the model
community, and then project it through time under different fishing
strategies.

Methods are included to explore the results, including plots and 
calculation of community indicators such as the slope of the size spectrum.
Size-based models can be complicated so mizer contains many default
options that can be easily changed by the user.

The package is on [CRAN](https://cran.r-project.org/package=mizer) and 
therefore available from R's build-it package manager.

See the accompanying [vignette](https://cran.r-project.org/web/packages/mizer/vignettes/mizer_vignette.pdf) 
for more details on how the package works, including detailed examples.

[![Rdoc](http://www.rdocumentation.org/badges/version/mizer)](http://www.rdocumentation.org/packages/mizer)

[![Travis-CI Build Status](https://travis-ci.org/sizespectrum/mizer.svg?branch=master)](https://travis-ci.org/sizespectrum/mizer)
