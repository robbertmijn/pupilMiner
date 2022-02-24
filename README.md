# pupilMiner
Process Pupil Data in R

## Installation

Install this package using devtools:

```
devtools::install_github("robbertmijn/pupilMiner")
```

Verify if installation was successful by parsing the example data files:

```
library(pupilMiner)
demofile <- system.file("extdata", "sub_1.asc", package = "pupilMiner")
dat <- parse_asc_file(demofile)
```

## Usage

Pass parameters to ```parse_asc_file``` to modify options. For example:

```keep_vars``` allows you to specify which variables from the OS-experiment will be retained
```samptime``` allows downsampling to a specific sample duration in ms
```timelock_to``` sets the time markers to 0 at the onset of a specified phase (in the example the onset of the "problem" phase)
```baseline_period``` subtracts average pupil size during this period from the rest of the trace per trial

For example:

```
library(pupilMiner)
demofile <- system.file("extdata", "sub_1.asc", package = "pupilMiner")
dat <- parse_asc_file(demofile, keep_vars = c("subject_nr", "difficulty"), samptime = 40, timelock_to = "problem", baseline_period = c(-500, 0))
```

This package is in development and mostly for personal use
