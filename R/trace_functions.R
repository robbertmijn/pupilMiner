#' trace processor
#'
#' @param data eyetracker table object with raw samples
#' @param blinkdat blink information
#' @return Table with processed traces and list of blinks per trial
#' @export
traceprocessor <- function(data, blinkparams = NULL, saccparams = NULL){
  # set general variables
  sampdur <- data$time[2] - data$time[1]

  if(is.null(blinkparams)){
    blinkparams <- list(vt = 5, maxdur = 500, margin = 100, smooth_winlength = 84, sampdur = sampdur)
  }
  if(is.null(saccparams)){
    saccparams <- list(mindur = 3, maxdur = 150, mindist = .1, maxdist = 100, msVthres = 6, smoothwin = 7, sampdur = sampdur)
  }
  data[, i := 1:.N, by = block]
  data$blink_id <- as.double(NA)
  data$sacc_freq <- as.double(0)

  # copy pupil, x and y traces
  data$pupil_raw <- data$ps
  data$x_raw <- data$xp
  data$y_raw <- data$yp

  # output parameters of the trace processor for user
  cat("Blinks from velocity, pars: vt =", blinkparams$vt,
      ", maxdur (ms) =", blinkparams$maxdur,
      ", sample rate =", 1000/sampdur,
      " Hz, margin (ms) =", blinkparams$margin,
      ", smooth (ms) =", blinkparams$smooth_winlength,
      "\nSaccades from gaze velocity, pars: vt =", saccparams$msVthres,
      ", maxdur (ms) =", saccparams$maxdur,
      ", mindur (ms) =", saccparams$mindur,
      " maxdist =", saccparams$maxdist,
      ", mindist =", saccparams$mindist,
      " msVthres =", saccparams$msVthres,
      ", smooth (ms) =", saccparams$smoothwin, "\n")

  # loop over all recorded trials (marked as "blocks" by eyelink)
  for(bl in unique(data$block)){

    # find blinks in current pupil trace trace (ps)
    blink_list <- getBlinkList(data[block == bl, ps], blinkparams)

    # Replace raw pupil/x/y traces with NA during blink periods
    if(!is.null(blink_list)){
      cat("[", bl, nrow(blink_list), "] ")
      for(r in 1:nrow(blink_list)){
        istart <- blink_list[r, 1]
        iend <- blink_list[r, 2]
        data[block == bl & i %between% c(istart, iend), ":="(ps = NA,
                                                             xp = NA,
                                                             yp = NA,
                                                             blink_id = r)]
      }
    }

    sacc_list <- getMicSaccList(data[block == bl, .(xp, yp)], saccparams)
    if(!is.null(sacc_list)){
      cat("{", bl, nrow(sacc_list), "} ")
      for(r in 1:nrow(sacc_list)){
        istart <- sacc_list[r, 1]
        iend <- sacc_list[r, 2]
        # .sacc_dist <- sqrt(sum((data[block == bl & i == istart, .(xp, yp)] - data[block == bl & i == iend, .(xp, yp)])^2, na.rm=T))
        data[block == bl & i %between% c(istart, iend), ":="(sacc_id = r,
                                                             sacc_dist = sacc_list[r, 3],
                                                             sacc_dur = (iend - istart) * sampdur)]
        data[block == bl & i == as.integer(mean(c(istart, iend))), ":="(sacc_freq = 1)]
      }
    }
  }
  data$i <- NULL
  cat("\n")

  # linearly interpolate missing data up to 500 ms
  data[, ps := na.approx(ps, maxgap = round(500/sampdur), na.rm = F), by = block]

  return(data)
}
#' return list of blinks
#'
#' @param strace pupil trace
#' @param by by
#' @param win number of samples to smooth over
#' @return pupil table object with column with smoothed signal
#' @export
getBlinkList <- function(strace, blinkparams){
  #unpack params
  sampdur <- blinkparams$sampdur
  vt <- blinkparams$vt
  maxdur <- round(blinkparams$maxdur/sampdur)
  margin <- round(blinkparams$margin/sampdur)
  smooth_winlength <- round(blinkparams$smooth_winlength/sampdur)

  strace <- smooth_trace(strace, win = smooth_winlength)
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
  return(lblink)
}
#' get Microsaccccade list
#'
#' @param xy_dat xy coordinates of gaze
#' @return A data.table with samples
#' @export
getMicSaccList <- function(xy_dat, saccparams){
  #unpack params
  sampdur <- saccparams$sampdur
  mindur <- round(saccparams$mindur/sampdur)
  maxdur <- round(saccparams$maxdur/sampdur)
  mindist <- saccparams$mindist
  maxdist <- saccparams$maxdist
  msVthres <- saccparams$msVthres
  smoothwin <- round(saccparams$smoothwin/sampdur)

  vxy_dat <- xy_dat[, .(vx = smooth_trace(xp - shift(xp), win = smoothwin),
                        vy = smooth_trace(yp - shift(yp), win = smoothwin),
                        i = 1:nrow(xy_dat))]
  vxy_dat[, v := sqrt(vx^2 + vy^2)]
  vt <- mad(vxy_dat$v, na.rm = T) * msVthres
  vtrace = vxy_dat$v
  lsacc <- NULL
  ifrom <- 0
  saccISI <- 50

  while(ifrom < length(vtrace)){
    # The onset of the saccade is the moment at which the x or y velocity
    # exceeds the threshold.
    l <- which(vtrace[ifrom:length(vtrace)] > vt)[1]

    if(is.na(l)){
      break # No sacc detected
    }

    istart = l[1] + ifrom
    if(ifrom == istart){
      break
    }

    # The end sacc is the moment at which the velocity drops below thres
    l = which(vtrace[istart:length(vtrace)] < vt)[1]
    if(!is.na(l)){
      iend <- l[1] + istart
    } else {
      iend <- length(vtrace)
    }

    # Keep a period between saccades
    ifrom <- iend + saccISI

    .sacc_dist <- sqrt(sum((xy_dat[istart, .(xp, yp)] - xy_dat[iend, .(xp, yp)])^2, na.rm=T))

    if((iend - istart) %between% c(mindur, maxdur) & .sacc_dist %between% c(mindist, maxdist)){
      lsacc <- rbind(lsacc, c(istart, iend, .sacc_dist))
    }
  }
  return(lsacc)
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
  return(data)
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
