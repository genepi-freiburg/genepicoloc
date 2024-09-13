#' @title Read sumstats 1
#' @param sumstats_file path to sumstats.
#' @param CHR_var chromosome (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats.
#' @export
query_sumstats_1 <- function(sumstats_file, CHR_var, BP_START_var, BP_STOP_var) {
  tabix_cmd <- paste0("tabix -h ", sumstats_file, " ", CHR_var, ":", BP_START_var, "-", BP_STOP_var)
  sumstats <- read.table(text=system(tabix_cmd, intern = T), sep = "\t", header = T)
  if (nrow(sumstats) == 0) { return(sumstats) }
  sumstats <- subset(sumstats, CHR == CHR_var & POS >= BP_START_var & POS <= BP_STOP_var)
  return(sumstats)
}

