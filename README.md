# pupilMiner
Process Pupil Data in R

Install this package using devtools:

```
devtools::install_github("robbertmijn/pupilMiner")
```

Verify if installation was succesful by parsing the example data files:

```
library(pupilMiner)
dat <- parse_asc_file(system.file("extdata", "sub_1.asc", package = "pupilMiner"))
```

This package is in development and mostly for personal use
