# tpmFun ------------------------------------------------------------------

## tpmFun: functions to convert counts to tpm

# inset function
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

# main function
tpmFun <-function(df){
  df %>%
    pivot_longer(7:ncol(df), names_to="sample", values_to="cnt") %>%
    group_by(sample) %>%
    mutate(tpm=tpm(cnt, Length)) %>%
    dplyr::select(-cnt) %>%
    pivot_wider(names_from=sample, values_from=tpm)
}

# rpkmFun -----------------------------------------------------------------

rpkm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(counts) * 1e9
}

rpkmFun <- function(df) {
  df %>% 
    pivot_longer(7:ncol(df), names_to="sample", values_to="cnt") %>%
    group_by(sample) %>%
    mutate(rpkm=rpkm(cnt, Length)) %>%
    select(-cnt) %>%
    pivot_wider(names_from=sample, values_from=rpkm)
}
