#' load and preprocess data from an asc file
#'
#' @param dat pupil table object
#' @param accept leave just one pupil trace instead of keeping intermediate traces
#' @param timelock_to string value of phase on which to timelock to
#' @return A data.table with samples
#' @export
preprocess <- function(file, keep_raw = T, timelock_to = NULL, keep_vars = NULL, vt = 15){
  dat <- parse_asc_file(file, keep_vars)
  # pb <- progress_bar$new(total = length(unique(dat$trial)))
  # dat[, pupil_r := pupil]
  # for(tr in unique(dat$trial)){
  #   # cat(tr, "... ")
  #   pb$tick()
  #   try(
  #     dat <- dat[trial == tr, pupil := blinkreconstruct(pupil, vt = vt)]  
  #   )
  # }
  # if(!keep_raw){
  #   dat[, pupil_r := NULL]
  # }
  if(!is.null(timelock_to)){
    dat <- timelock(dat, timelock_to)
  }
  return(dat)
}
#' This functions load an ASC file
#'
#' @param infile Path to the input file
#' @return A data.table with samples
#' @export
parse_asc_file <- function(infile, keep_vars = NULL, blinkreconstruct = T){
  # print message to console
  cat("Loading ", infile, "\n")
  # use eyelinker package to parse asc file
  dat <- read.asc(infile)
  # combine raw data with the messages we sent during the experiment in a table
  tempraw <- data.table(dat$raw)
  if(blinkreconstruct){
    tempblinks <- data.table(dat$blinks)
    for(r in 1:nrow(tempblinks)){
      stime <- tempblinks[r]$stime - 50
      etime <- tempblinks[r]$etime + 50
      tempraw[time %between% c(stime, etime), ps := NA]
    }
  }
  tempphases <- data.table(dat$msg, type = "phase")[startsWith(text, "start_phase")]
  tempdat <- rbind(tempraw, tempphases, fill = T)
  # each "start_recording" (at beginning of a trial) is counted in the variable block
  setnames(tempdat, c("block", "ps", "xp", "yp"), c("trial", "pupil", "x", "y"))
  # sort the table based on time points (messages will be sorted with raw data)
  setkey(tempdat, time)
  # locf (last observation carries forward) will add info about messages to raw data
  tempdat[, phase := na.locf(text, na.rm = F), by = trial]
  # Only keep the actual raw data, which now has info about the phase and trial
  tempdat <- tempdat[is.na(type)]
  # remove the "start_phase " bit from the phase name
  tempdat[, phase := gsub("start_phase ", "", phase)]
  # calculate time since trial start and time since phase start
  tempdat[, t_exp := time]
  tempdat[, time := time - min(time), by = trial]
  # tempdat[, t_phase := time - min(time), by = .(trial, phase)]
  # remove unused variables and extract pp number from filename
  tempdat[, ":="(pp = gsub("\\D", "", tail(strsplit(infile, "/")[[1]], 1)), cr.info = NULL, text = NULL, type = NULL, input = NULL)]
  if(!is.null(keep_vars)){
    tempvars <- data.table(dat$msg)[startsWith(text, "var")]
    tempvars <- tempvars[, .(trial = block, 
                             var = strsplit(text, " ")[[1]][2], 
                             value = strsplit(text, " ")[[1]][3]), by = 1:nrow(tempvars)]
    tempvars <- data.table(dcast(tempvars, trial ~ var, value.var = "value"))
    tempvars <- tempvars[, ..keep_vars]
    tempdat <- merge(tempdat, tempvars, by = "trial")
  }
  return(tempdat)
}
#' Export a blink report
#'
#' @param dat table with pupil data, including raw pupil
#' @export
blink_report <- function(dat, folder = "blink_reports", trials_per_page = 20, period = c(-Inf, Inf)){
  theme_set(theme_classic())
  dir.create(file.path(getwd(), folder), showWarnings = FALSE)
  pdf(file = paste0(folder, "/blinks_pp", unique(dat$pp), ".pdf"))
  dat[, page := trial %/% trials_per_page]
  pb <- progress_bar$new(total = length(unique(dat$page)))
  for(p in unique(dat$page)){
    pb$tick()
    print(
      ggplot(dat[page == p & time %between% period], aes(x = time, group = trial)) +
        geom_line(aes(y = pupil_r), alpha = .4) +
        geom_line(aes(y = pupil), alpha = .4, color= "red") +
        facet_wrap(~trial) +
        labs(title = paste0("pp ", unique(dat$pp), " page ", p, " of ", max(dat$page)))
    )
  }
  dev.off()
}
#' load the demo data
#'
#' @return A data.table with samples
#' @export
DemoDat <- function(){
  return(parse_asc_file(system.file("extdata", "sub_1.asc", package = "pupilMiner")))
}
