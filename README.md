# Server ROC-GLM Functions for DataSHIELD

## Installation

Install the package on the server where DataSHIELD is running via:

```r
devtools::install.github("schalkdaniel/ds.roc.glm")
```

## Client Functions

The client functions are in a saparated file `R/client_functions.R`. Since these functions are exported in another package to make them available more easily, it is required to source this file everytime you want to use the ROC-GLM:

```r
source("R/client_functions.R")
```
