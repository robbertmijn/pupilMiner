#' smooth a trace
#'
#' @param dat pupil table object
#' @param by by
#' @param win number of samples to smooth over
#' @return pupil table object with column with smoothed signal
#' @export
smooth_trace <- function(trace, win = 21){
  if(win %% 2 != 0){
    # make winlength an uneven
    win <- win + 1
  }
  # pad around the signal half the width of the window, they will fall off when smoothing
  d <- (win-1)/2
  s <- c(trace[2:(d+1)], trace, rev(tail(trace, d+1))[-1])
  return(rollmean(s, win, align = "center", na.rm = T))
}
#' detect blinks
#'
#' @param dat pupil table object
#' @param vt A pupil velocity threshold. Lower thresholds more easily trigger blinks.
#' @param maxdur The maximum duration (in samples) for a blink. Longer blinks are not reconstructed.
#' @param margin The margin to take around missing data
#' @param smooth_winlength The width of the smoothing window. This should be an odd integer.
#' @param std_thr std_thr
#' @return A data.table with samples
#' @export
blinkreconstruct <- function(trace, vt = 15, maxdur = 500, margin = 10, smooth_winlength = 21, std_thr = 3){
  strace <- smooth_trace(trace, win = smooth_winlength)
  vtrace <- tail(strace, length(strace)-1) - head(strace, length(strace)-1)
  
  ifrom <- 0
  lblink <- NULL
  
  while(1){
    # The onset of the blink is the moment at which the pupil velocity
    # exceeds the threshold.
    l <- which(vtrace[ifrom:length(vtrace)] < -vt)[1]
    if(is.na(l)){
      break # No blink detected
    }
    istart = l[1] + ifrom
    if(ifrom == istart){
      break
    }
    # The reversal period is the moment at which the pupil starts to dilate
    # again with a velocity above threshold.
    l = which(vtrace[istart:length(vtrace)] > vt)[1]
    if(is.na(l)){
      ifrom = istart
      next
    }
    imid = l[1] + istart
    # The end blink period is the moment at which the pupil velocity drops
    # back to zero again.
    l = which(vtrace[imid:length(vtrace)] < 0)[1]
    if(is.na(l)){
      ifrom = imid
      next
    }
    iend = l[1] + imid
    ifrom = iend
    # We generally underestimate the blink period, so compensate for this
    if(istart - margin >= 0){
      istart <- istart - margin
    }
    if(iend + margin < length(trace)){
      iend <- iend + margin
    }
    # We don't accept blinks that are too long, because blinks are not
    # generally very long (although they can be).
    if(iend - istart > maxdur){
      ifrom = istart + round(maxdur/10)
      next
    }
    lblink <- rbind(lblink, c(istart, iend))
  }
  
  # Now reconstruct the trace during the blinks
  if(!is.null(lblink)){
    for(i in 1:nrow(lblink)){
      istart <- lblink[i, 1]
      iend <- lblink[i, 2]
      # First create a list of (when possible) four data points that we can
      # use for interpolation.
      dur <- iend - istart
      l = NULL
      
      l <- c(l, istart, iend)
      if(istart - dur >= 0 & 
         iend + dur < length(strace) & 
         !is.na(trace[iend + dur])){
        if(!is.na(trace[istart - dur])){
          l <- c(istart - dur, l)
          l <- c(l, iend + dur)
        }
      }
      x = l
      # If the list is long enough we use cubic interpolation, otherwise we
      # use linear interpolation
      y = trace[x]
      # cat(x, "on", y, "\n")
      xInt = istart:iend
      if(!any(is.na(y))){
        if(length(x) >= 4){
          yInt = interp1(x, y, xInt, method = "cubic")
        } else{
          yInt = interp1(x, y, xInt)
        }
        trace[xInt] <- yInt
      }
    }
  }
  # # For all remaining gaps, replace them with the previous sample if
  # # available
  # b = which(
  #   trace < (mean(trace, na.rm = T) - std_thr * sd(trace, na.rm = T))
  #   | trace > (mean(trace, na.rm = T) + std_thr * sd(trace, na.rm = T))
  #   | is.na(trace)
  # )
  # for(i in b){
  #   if (i == 0){
  #     next
  #   }
  #   trace[i] = trace[i - 1]
  # }
  return(trace)
}
#' Timelock trace
#'
#' @param dat pupil table object
#' @param phase_name of which the start will be new 0 point
#' @param by by
#' @return A data.table with samples
#' @export
timelock <- function(dat, phase_name){
  dat <- merge(dat, dat[phase == phase_name, .(onsetTime = min(time)), by = .(pp, trial)], by = c("pp", "trial"))
  # Timelock the data to the onset of the phase
  dat[, time := time - onsetTime, by = trial] 
  # Remove the temporary variables onsetTime and Baseline
  dat[, onsetTime := NULL]
}
#' baseline
#'
#' @param dat pupil table object
#' @param baselineRange timepoints to baseline on
#' @param by by
#' @return A data.table with samples
#' @export
baseline <- function(dat, baselineRange = c(-200, 0), by = NULL){
  byArgument <- unique(c(by, "pp", "trial"))
  dat <- merge(dat, dat[time %between% baselineRange, .(baseline = mean(pupil, na.rm = T)),
                        by = byArgument], by = byArgument)
  dat[, pupil := pupil - baseline]
  dat[, baseline := NULL]
  return(dat)
}
#' Downsample
#'
#' @param dat pupil table object
#' @param by by
#' @param Hz downsample to
#' @return A data.table with samples
#' @export
downsample <- function (dat, by, Hz = 100){
  sampleTime <- dat[, time[2] - time[1]]
  binSize <- 1000/Hz
  if (binSize%%sampleTime != 0) {
    warning("Sample frequency of data is not a multiple of the target frequency specified in the by argument")
  }
  dat$DS <- dat$time%/%binSize
  allF <- c(by, "DS")
  dat <- dat[, list(pupil = median(pupil), x = median(x), y = median(y),
                    phase = head(phase, 1)), by = allF]
  dat$time <- dat$DS * binSize
  dat$DS <- NULL
  return(dat)
}


