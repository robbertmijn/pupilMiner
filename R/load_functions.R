#' Load a Matrix
#'
#' This functions load an ASC file
#'
#' @param infile Path to the input file
#' @return A data.table with samples
#' @export
parse_asc_file <- function(infile){
  
  # print message to console
  cat("Loading ", infile, "\n")
  
  # use eyelinker package to parse asc file
  dat <- read.asc(infile)
  
  # combine raw data with the messages we sent during the experiment in a table
  tempraw <- data.table(dat$raw)
  tempmsg <- data.table(dat$msg, type = "msg")[startsWith(text, "start_phase")]
  tempdat <- rbind(tempraw, tempmsg, fill = T)
  
  # each "start_recording" (at beginning of a trial) is counted in the variable block
  setnames(tempdat, c("block", "ps", "xp", "yp"), c("trial", "pupil", "x", "y"))
  
  # sort the table based on time points (messages will be sorted with raw data)
  setkey(tempdat, time)
  
  # locf (last observation carries forward) will add info about messages to raw data
  tempdat[, phase := na.locf(text, na.rm = F), by = trial]
  
  # Only keep the actual raw data, which now has info about the phase and trial
  tempdat <- tempdat[!is.na(phase) & is.na(type)]
  
  # remove the "start_phase " bit from the phase name
  tempdat[, phase := gsub("start_phase ", "", phase)]
  
  # calculate time since trial start and time since phase start
  tempdat[, t_exp := time]
  tempdat[, time := time - min(time), by = trial]
  tempdat[, t_phase := time - min(time), by = .(trial, phase)]
  
  # remove unused variables and extract pp number from filename
  tempdat[, ":="(pp = gsub("\\D", "", infile), cr.info = NULL, text = NULL, type = NULL, input = NULL)]
  
  return(tempdat)
}