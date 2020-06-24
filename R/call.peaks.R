# call peaks from density and return data frame of with x,y location of peaks
call.peaks <- function(den, neighlim=0.05) {
  # call peaks from density
  # return data frame with x,y location of peaks
  peaks.input <- cbind(x=den$x, y=den$y)
  peaks <- peakpick(mat=peaks.input, neighlim=neighlim)
  peaks.input[peaks[, 2], ] %>%
    as_tibble() %>%
    mutate(x=round(x, 2)) %>%
    group_by(x) %>%
    summarise(y=round(median(y), 2)) %>%
    ungroup() %>%
    group_by(y) %>%
    summarise(x=round(median(x), 2)) %>%
    ungroup() %>%
    arrange(desc(y)) %>%
    dplyr::select(x, y) %>%
    as.data.frame()
}
