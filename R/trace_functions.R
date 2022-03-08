#' trace processor
#'
#' @param data eyetracker table object with raw samples
#' @param blinkdat blink information
#' @return Table with processed traces and list of blinks per trial
#' @export
traceprocessor <- function(data, vt = 5, maxdur = 500, margin = 100, smooth_winlength = 84){

  # copy pupil, x and y traces
  data$pupil_raw <- data$ps
  data$x_raw <- data$xp
  data$y_raw <- data$yp

  # Remove samples collected during blink period marked by velocity threshold
  data <- TP_velocity(data, vt = vt, maxdur = maxdur, margin = margin, smooth_winlength = smooth_winlength)

  sampdur <- data$time[2] - data$time[1]
  # linearly interpolate missing data up to 500 ms
  data[, ps := na.approx(ps, maxgap = round(500/sampdur), na.rm = F), by = block]

  return(data)
}
#' Extracts blinks from eyelink data
#'
#' @param dat pupil table object
#' @param margin The margin to take around missing data
#' @return A data.table with samples
#' @export
TP_eyelink <- function(data, blinkdat, margin = 100){
  # Remove data during and around blinks using eyelink data
  for(r in 1:nrow(blinkdat)){
    stime <- blinkdat[r]$stime
    etime <- blinkdat[r]$etime
    data[time %between% c(stime - margin, etime + margin), ":="(ps = NA,
                                                                xp = NA,
                                                                yp = NA)]
  }
  cat("Detected,", nrow(blinkdat), "blinks from eyelink.\n")
  return(data)
}
#' Blink detection using method described by Mathot
#' @param data pupil table object
#' @param vt A pupil velocity threshold. Lower thresholds more easily trigger blinks.
#' @param maxdur The maximum duration (in samples) for a blink. Longer blinks are not reconstructed.
#' @param margin The margin to take around missing data
#' @param smooth_winlength The width of the smoothing window. This should be an odd integer.
#' @param std_thr std_thr, TODO
#' @return A data.table with blink periods marked as NA
#' @export
TP_velocity <- function(data, vt = 5, maxdur = 500, margin = 100, smooth_winlength = 84){

  sampdur <- data$time[2] - data$time[1]
  margin <- round(margin/sampdur)
  maxdur <- round(maxdur/sampdur)
  smooth_winlength <- round(smooth_winlength/sampdur)
  data[, i := 1:.N, by = block]
  data$blink_id <- as.integer(NA)

  cat("Blinks from velocity, pars: vt =", vt,
      ", maxdur (samples) =", maxdur,
      ", sample rate =", 1000/sampdur,
      " Hz, margin (samples) =", margin,
      ", smooth (samples) =", smooth_winlength, "\n")

  for(bl in unique(data$block)){
    strace <- data[block == bl, smooth_trace(ps, win = smooth_winlength)]
    vtrace <- strace - shift(strace)

    lblink <- NULL
    ifrom <- 0

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
      if(iend + margin < length(strace)){
        iend <- iend + margin
      }
      # We don't accept blinks that are too long, because blinks are not
      # generally very long (although they can be).
      # RM: we don't interpolate, but we do want to detect these long blinks!
      # if(iend - istart > maxdur){
      #   ifrom = istart + maxdur
      #   next
      # }
      lblink <- rbind(lblink, c(istart, iend))
    }

    # Replace raw pupil/x/y traces with NA during blink periods
    if(!is.null(lblink)){
      cat("[", bl, nrow(lblink), "] ")
      for(r in 1:nrow(lblink)){
        istart <- lblink[r, 1]
        iend <- lblink[r, 2]
        data[block == bl & i %between% c(istart, iend), ":="(ps = NA,
                                                             xp = NA,
                                                             yp = NA,
                                                             blink_id = r)]
      }
    }
  }
  data$i <- NULL
  cat("\n")
  return(data)
}
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
#' Timelock trace
#'
#' @param dat pupil table object
#' @param phase_name of which the start will be new 0 point
#' @param by by
#' @return A data.table with samples
#' @export
timelock <- function(data, phase_name){
  data <- merge(data, data[phase == phase_name, .(onsetTime = min(time)), by = trial], by = "trial")
  # Timelock the data to the onset of the phase
  data[, time := time - onsetTime, by = trial]
  # Remove the temporary variables onsetTime and Baseline
  data[, onsetTime := NULL]
}
#' baseline
#'
#' @param dat pupil table object
#' @param baselineRange timepoints to baseline on
#' @param by by
#' @return A data.table with samples
#' @export
baseline <- function(data, baselineRange = c(-200, 0), by = NULL){
  byArgument <- unique(c(by, "trial"))
  data <- merge(data, data[time %between% baselineRange, .(baseline = mean(pupil, na.rm = T)),
                           by = byArgument], by = byArgument)
  data[, pupil := ifelse(is.nan(baseline), NA, pupil - baseline)]
  return(data)
}
#' Downsample
#'
#' @param dat pupil table object
#' @param samptime New sample duration
#' @param by by
#' @param ds_cols Columns with data that require downsampling
#' @return A data.table with samples
#' @export
downsample <- function(data, samptime, By = "trial", ds_cols = NULL){
  data[, DS := time %/% samptime, by = By]
  By <- c(By, "DS")
  for(col in ds_cols){
    if(col %in% names(data)){
      data[, (col) := median(get(col), na.rm = T), by = By]
    } else {
      warning(col, " not found in data")
    }
  }
  data$time <- data$DS * samptime
  data$DS <- NULL
  data <- unique(data)
  return(data)
}
