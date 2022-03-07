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
# get path of demofile which is installed with pupilMiner
demofile <- system.file("extdata", "sub_1.asc", package = "pupilMiner")
# Parse the datafile
dat <- parse_asc_file(demofile,
                       keep_vars = c("subject_nr", "difficulty"),
                       samptime = 40,
                       timelock_to = "problem",
                       baseline_period = c(-100, 0))

# Pupil size during blinks up to 500ms are interpolated. 
# For visualizing, we'll interpolate the ones that are not interpolated to a temporary variable
dat[, pupil_not_interpolated := na.approx(pupil, na.rm = F), by = trial]

# Plot all trials in one plot. Red segments are interpolated blinks, blue segments are blinks
# that are too long to be interpolated and are left as NA in the data.
ggplot(dat, aes(x = time, y = pupil, gr = as.factor(trial))) +
  geom_line(alpha = .3) +
  geom_line(data = dat[!is.na(blink_id)], aes(y = pupil_not_interpolated), color = "blue") +
  geom_line(data = dat[!is.na(blink_id)], color = "red") 
```

This package is in development and mostly for personal use
