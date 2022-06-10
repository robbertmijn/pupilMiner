#' Export a trace processing report for a single pp
#'
#' @param dat table with data, including raw pupil/x/y data
#' @param folder where to save reports
#' @param tag add this tag to filenames
#' @param trials_per_page trials per page
#' @param period start and end time to plot
#' @export
trace_reports <- function(dat,
                          folder = "trace_reports",
                          tag = "traces",
                          trials_per_page = 20,
                          period = c(-Inf, Inf),
                          incl_sacc = T){
  theme_set(theme_classic())
  dir.create(file.path(getwd(), folder), showWarnings = FALSE)
  pdf(file = paste0(folder, "/", tag, "_subject_nr", unique(dat$subject_nr), ".pdf"))
  dat[, page := trial %/% trials_per_page]
  pb <- progress_bar$new(total = length(unique(dat$page)))
  for(p in unique(dat$page)){
    pb$tick()
    print(
      ggplot(dat[page == p & time %between% period], aes(x = time)) +
        geom_line(aes(y = pupil_raw), color = "red") +
        geom_path(data = dat[page == p & time %between% period & !is.na(blink_id)], aes(y = pupil_raw, grp = as.factor(blink_id)), color = "blue") +
        facet_wrap(~trial, scales = "free") +
        labs(title = paste0("RAW: subject_nr ", unique(dat$subject_nr), ", page ", p, " of ", max(dat$page)))
    )
    print(
      ggplot(dat[page == p & time %between% period], aes(x = time)) +
        geom_line(aes(y = pupil)) +
        facet_wrap(~trial) +
        labs(title = paste0("PROCESSED: subject_nr ", unique(dat$subject_nr), ", page ", p, " of ", max(dat$page)))
    )
    print(
      ggplot(dat[page == p & time %between% period]) +
        geom_path(aes(x = x, y = y, color = time), alpha = .3) +
        geom_path(data = dat[page == p & time %between% period & !is.na(sacc_id)], aes(y = y, x = x, grp = as.factor(sacc_id)), color = "red") +
        facet_wrap(~trial, scales = "free") +
        labs(title = paste0("MICROSACCADES: subject_nr ", unique(dat$subject_nr), ", page ", p, " of ", max(dat$page)))
    )

  }
  dev.off()
  dat[, page := NULL]
}

