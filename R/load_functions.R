#' This functions load an ASC file
#'
#' @param infile Path to the input file
#' @param keep_vars column names to keep from vars that are passed as messages to the eyelink, allows you to specify which variables from the OS-experiment will be retained.
#' @param TP Specify which type of trace processor
#' @param TP_report_dir if passed, pdfs with trace reports are saved in this folder
#' @param TP_report_period select begin and end point to report
#' @param timelock_to if passed, time variable is set to 0 at onset of a phase
#' @param samptime if passed, data is downsampled to this sample time (ms)
#' @param baseline_period pupil trace baseline
#' @return A data.table with samples
#' @export
parse_asc_file <- function(infile,
                           keep_vars = NULL,
                           TP_report_dir = NULL,
                           TP_report_period = c(-Inf, Inf),
                           timelock_to = NULL,
                           samptime = NULL,
                           baseline_period = c(NULL, NULL)){

  cat("=======\ntraceprocess", format(Sys.time(), "%Y%m%d%H%M%S"), "\n")
  # print message to console
  cat("Loading ", infile, "\n")
  # use eyelinker package to parse asc file
  dat <- read.asc(infile)

  # extract tables from eyelinker object
  tempraw <- data.table(dat$raw)
  tempphases <- data.table(dat$msg, type = "phase")[startsWith(text, "start_phase")]
  tempvars <- data.table(dat$msg)[startsWith(text, "var")]
  tempblinks <- data.table(dat$blinks)
  tempfix <- data.table(dat$fix)
  tempsacc <- data.table(dat$sacc)

  # apply trace processor
  TP_result <- traceprocessor(tempraw, tempblinks)

  # add variables and time lock to start of trial
  TP_result <- add_variables(TP_result, tempphases, tempvars, keep_vars)

  # generate pdf with trace processing report
  if(!is.null(TP_report_dir)){
    cat("creating pdf report\n")
    trace_reports(TP_result, folder = TP_report_dir, period = TP_report_period)
  }
  if(!is.null(timelock_to)){
    cat("Timelocking to", timelock_to, "\n")
    TP_result <- timelock(TP_result, timelock_to)
  }

  if(!is.null(baseline_period)){
    cat("Subtracting baseline", "\n")
    TP_result <- baseline(TP_result, baseline_period)
  }

  if(!is.null(samptime)){
    cat("Downsampling to", 1000/samptime, "Hz ... from", nrow(TP_result))
    TP_result <- downsample(TP_result, samptime, ds_cols = c("x", "y", "pupil", "pupil_raw", "x_raw", "y_raw", "t_exp"))
    cat(" to", nrow(TP_result), "rows\n")
  }
  cat("\n\n")
  return(TP_result)
}
#' Extract variables and add
#'
#' @param dat table with data, including raw pupil/x/y data
#' @param vars table with var messages
#' @param phases table with phase messages
#' @param keep_vars column names to keep from vars
#' @export
add_variables <- function(dat, phases, vars, keep_vars = NULL){
  dat <- rbind(dat, phases, fill = T)
  # sort the table based on time points (messages will be sorted with raw data)
  setkey(dat, time)
  # locf (last observation carries forward) will add info about messages to raw data
  dat[, phase := na.locf(text, na.rm = F), by = block]
  # Only keep the actual raw data, which now has info about the phase and trial
  dat <- dat[is.na(type)]
  # remove the "start_phase " bit from the phase name
  dat[, phase := gsub("start_phase ", "", phase)]
  # calculate time since trial start and time since phase start
  dat[, t_exp := time]
  dat[, time := time - min(time), by = block]
  # remove unused variables and extract pp number from filename
  # dat[, ":="(pp = gsub("\\D", "", tail(strsplit(infile, "/")[[1]], 1)), cr.info = NULL, text = NULL, type = NULL, input = NULL)]
  dat[, ":="(cr.info = NULL, text = NULL, type = NULL, input = NULL)]
  # extract variable names passed with keep_vars from messages
  if(!is.null(keep_vars)){
    keep_vars <- c("block", keep_vars)
    vars <- vars[, .(block,
                     var = strsplit(text, " ")[[1]][2],
                     value = strsplit(text, " ")[[1]][3]), by = 1:nrow(vars)]
    vars <- data.table(dcast(vars, block ~ var, value.var = "value"))
    vars <- vars[, ..keep_vars]
    dat <- merge(dat, vars, by = "block")
  }
  setnames(dat, c("block", "xp", "yp", "ps"), c("trial", "x", "y", "pupil"))

  return(dat)
}
#' load the demo data
#'
#' @return A data.table with samples
#' @export
DemoDat <- function(){
  return(parse_asc_file(system.file("extdata", "sub_1.asc", package = "pupilMiner")))
}
