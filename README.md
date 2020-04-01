Eradicate
================

An R package for the analysis of animal abundance using data from common
survey methods employed during pest animal eradication programs.
Currently, the package contains models that can handle data collected
from occupancy surveys, encounters from remote camera surveys and
removal sampling. Parameters can be modelled with covariates, in the
style of the package ‘unmarked’. Indeed, the package uses some of the
model fitting code from ‘unmarked’. This package is primarily meant to
be used as a backend to the shiny app “eradication tools decision
support”.

## Installation

Currently the package is only available on GitHub. Install the package
using

`install_github("eradicate-dev/eradicate", build_vignettes=TRUE)`

The vignette contains a fully worked case study with examples of all the
current models. This can be accessed by

`browseVignettes(package="eradicate")`

## bug reports

The package is a work in progress so please use the [issue
tracker](https://github.com/dslramsey/eradicate/issues) to report bugs
in the code or models.
