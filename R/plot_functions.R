#' Export a trace processing report
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
                          period = c(-Inf, Inf)){
  theme_set(theme_classic())
  dir.create(file.path(getwd(), folder), showWarnings = FALSE)
  pdf(file = paste0(folder, "/", tag, "_pp", unique(dat$pp), ".pdf"))
  dat[, page := trial %/% trials_per_page]
  pb <- progress_bar$new(total = length(unique(dat$page)))
  for(p in unique(dat$page)){
    pb$tick()
    print(
      ggplot(dat[page == p & time %between% period], aes(x = time, group = trial)) +
        geom_line(aes(y = pupil), size = .5, color= "red") +
        geom_line(aes(y = pupil_raw), size = .2, ) +
        facet_wrap(~trial) +
        labs(title = paste0("pp ", unique(dat$pp), " page ", p, " of ", max(dat$page)))
    )
  }
  dev.off()
  dat[, page := NULL]
}
