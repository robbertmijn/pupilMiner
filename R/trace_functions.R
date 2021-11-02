#' smooth a trace
#'
#' @param dat pupil table object
#' @param by by
#' @param win number of samples to smooth over
#' @return pupil table object with column with smoothed signal
#' @export
smooth_trace <- function(dat, by = NULL, win = 10){
  gr_by <- c(by, "pp", "trial")
  # Deduce sampling frequency from time stamps
  recHz <- 1000/diff(dat$time[1:2])
  # smooth trace with rolling average of 15 ms
  dat[, pupil_s := rollmean(pupil, round(win/(recHz/1000)), align = "center", na.rm = F, fill = NA), by = gr_by]
  return(dat)
}
#' detect blinks
#'
#' @param dat pupil table object
#' @param maxDeltaDilation dilation threshold
#' @return A data.table with samples
#' @export
detect_blinks <- function(dat, maxDeltaDilation = 16){
  # Mark all samples of which the difference with the next sample exceeds a threshold as NA
  dat[, pupil_b := ifelse(abs(shift(pupil_s, -1) - pupil) > maxDeltaDilation, NA,
                          pupil_s), by = .(pp, trial)]
  return(dat)
}
#' Expand Blinks
#'
#' @param dat pupil table object
#' @param rejectionWindow ms padding
#' @return A data.table with samples
#' @export
expand_blinks<-function (dat, rejectionWindow = 50){
  # Deduce sampling frequency from time stamps
  recHz <- 1000/diff(dat$time[1:2])
  rejRegion <- rejectionWindow/1000 * recHz
  .expandblinks <- function(pupil_b) {
    rej <- which(is.na(pupil_b))
    rej <- outer(rej, -rejRegion:rejRegion, FUN = "+")
    rej <- rej[rej > 0 & rej <= length(pupil_b)]
    pupil_b[rej] <- NA
    return(pupil_b)
  }
  dat[, pupil_e := .expandblinks(pupil_b), by = .(trial, pp)]
  return(dat)
}
#' interpolate blinks
#'
#' @param dat pupil table object
#' @param type type of interpolation
#' @return A data.table with samples
#' @export
interpolateblinks <- function (dat, type = "linear", maxgap = 500) {
  recHz <- 1000/diff(dat$time[1:2])
  if (type != "linear") {
    stop("The only interpolation currently supported is linear interpolation")
  }
  .pupil.na.approx <- function(pupil_e) {
    if (all(is.na(pupil_e))) {
      return(pupil_e)
    }
    return(na.approx(pupil_e, na.rm = FALSE, maxgap = 500/1000*recHz))
  }
  dat[, pupil_i := .pupil.na.approx(pupil_e), by = .(trial, pp)]
  return(dat)
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


