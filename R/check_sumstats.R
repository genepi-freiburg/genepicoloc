#' check_sumstats
#' check input sumstats for possible issues
check_sumstats <- function(sumstats, Name = "Name") {
  attributes(sumstats)$QC <- ""
  # Duplicates
  Name_dup <- duplicated(sumstats[[Name]])
  if (any(Name_dup)) {
    sumstats <- subset(sumstats, !Name_dup)
    attributes(sumstats)$QC <- paste0(attributes(sumstats)$QC, "sumstats_2_duplicated_names")
  }
  return(sumstats)
}

#' check_sumstats_full
#' check input sumstats for possible issues
check_sumstats_full <- function(sumstats, fix=T,
                           AF = "AF", BETA = "BETA", SE = "SE", P = "P",
                           Name = "Name") {
  AF_1 <- sumstats[[AF]] == 1
  if (any(AF_1)) {warning("AF = 1 detected")}
  AF_0 <- sumstats[[AF]] == 0
  if (any(AF_0)) {warning("AF = 0 detected")}
  # INF
  BETA_INF <- is.infinite(sumstats[[BETA]])
  if (any(BETA_INF)) {warning("Infinite BETA detected")}
  SE_INF <- is.infinite(sumstats[[SE]])
  if (any(SE_INF)) {warning("Infinite SE detected")}
  P_INF <- is.infinite(sumstats[[P]])
  if (any(P_INF)) {warning("Infinite P detected")}
  # 0
  BETA_0 <- sumstats[[BETA]] == 0
  if (any(BETA_0)) {warning("BETA = 0 detected")}
  SE_0 <- sumstats[[SE]] == 0
  if (any(SE_0)) {warning("SE = 0 detected")}
  P_0 <- sumstats[[P]] == 0
  if (any(P_0)) {warning("P = 0 detected")}
  # Duplicates
  Name_dup <- duplicated(sumstats[[Name]])
  if (any(Name_dup)) {warning("Duplicated Names detected")}
  # To fix
  if (fix) {
    cols_to_check <- c()
    if (length(AF_1 != 0)) { sumstats <- subset(sumstats, !AF_1) }
    if (length(AF_0 != 0)) { sumstats <- subset(sumstats, !AF_0) }
    if (length(BETA_INF != 0)) { sumstats <- subset(sumstats, !BETA_INF) }
    if (length(SE_INF != 0)) { sumstats <- subset(sumstats, !SE_INF) }
    if (length(P_INF != 0)) { sumstats <- subset(sumstats, !P_INF) }
    if (length(BETA_0 != 0)) { sumstats <- subset(sumstats, !BETA_0) }
    if (length(SE_0 != 0)) { sumstats <- subset(sumstats, !SE_0) }
    if (length(P_0 != 0)) { sumstats <- subset(sumstats, !P_0) }
    if (length(Name_dup != 0)) { sumstats <- subset(sumstats, !Name_dup) }
  }
  return(sumstats)
}

#' sumstats_QC
sumstats_QC <- function(sumstats, BETA="BETA", SE="SE", P="P") {
  l1 <- !is.na(sumstats[[BETA]])
  l2 <- !is.na(sumstats[[SE]])
  l3 <- !is.na(sumstats[[P]])
  l4 <- !is.infinite(sumstats[[BETA]])
  l5 <- !is.infinite(sumstats[[SE]])
  l_all <- l1 & l2 & l3 & l4 & l5
  sumstats <- subset(sumstats, l_all)
  return(sumstats)
}

