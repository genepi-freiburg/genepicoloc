#' Read and format input summary statistics (under development)
#' @param sumstats_file path to sumstats file.
#' @param sumstats sumstats data.frame loaded into R.
#' @param Name update development
#' @param CHR update development
#' @param POS update development
#' @param A1 update development
#' @param A2 update development
#' @param BETA update development
#' @param SE update development
#' @param nlog10p_value update development
#' @param AF update development
#' @param N update development
#' @param other_columns vector with names of other columns
#' @return data frame with formatted columns
#' @examples
#' under development
#' @export
read_sumstats <- function(sumstats, sumstats_file=NULL,
                          Name, Name_new = "Name",
                          A1, A1_new = "A1",
                          BETA, BETA_new = "BETA",
                          SE, SE_new = "SE",
                          nlog10p_value, nlog10p_value_new = "nlog10P",
                          rsID = NULL, rsID_new = "rsID",
                          CHR = NULL, CHR_new = "CHR",
                          POS = NULL, POS_new = "POS",
                          A2 = NULL, A2_new = "A2",
                          AF = NULL, AF_new = "AF",
                          N = NULL, N_new = "N",
                          p_value=NULL,
                          other_columns = NULL) {
  if (!is.null(p_value)) {
    stop(paste("This function now accepts only nlog10p_value, please add it to sumstats first.\n",
               " Please ensure that underflow is handled properly (use handle_overflow() if needed)."))
  }
  if (!is.null(sumstats_file)) {
    stop(paste0("The function supports only 'sumstats' input.",
                "Please read the file into memory first and then pass it here."))
  }
  ### Checks
  # TODO identify necessary and optional columns
  if (missing(Name)) {stop("Name is required as CHR:POS:REF:ALT, please create it first.")}
  if (missing(A1)) {stop("A1 (effect allele) is required.")}
  if (missing(BETA)) {stop("BETA is required.")}
  if (missing(SE)) {stop("SE is required.")}
  if (missing(nlog10p_value)) {stop("nlog10p_value is required.")}
  ### Read
  sumstats <- as.data.frame(sumstats)
  ### Select
  all_cols <- c(Name, rsID, CHR, POS,
                A1, A2, BETA, SE,
                nlog10p_value, AF, N, other_columns)
  if (!all(all_cols %in% colnames(sumstats))) {
    stop("Not all provided colnames match colnames in the sumstats")
  }
  sumstats <- sumstats[,all_cols]
  ### Mandatory columns
  colnames(sumstats)[colnames(sumstats) == Name] <- Name_new
  colnames(sumstats)[colnames(sumstats) == A1] <- A1_new
  if (class(sumstats[[BETA]]) != "numeric") {
    stop("Beta column is not of class 'numeric', please check the data")
  }
  colnames(sumstats)[colnames(sumstats) == BETA] <- BETA_new
  if (class(sumstats[[SE]]) != "numeric") {
    warning("SE column is not of class 'numeric', converting to numeric")
    sumstats[[SE]] <- as.numeric(sumstats[[SE]])
  }
  colnames(sumstats)[colnames(sumstats) == SE] <- SE_new
  if (class(sumstats[[nlog10p_value]]) != "numeric") {
    stop(paste0("The nlog10p column is not 'numeric' class, please check sumstats."))
  }
  colnames(sumstats)[colnames(sumstats) == nlog10p_value] <- nlog10p_value_new
  ### Optional columns
  if (!is.null(rsID)) {
    colnames(sumstats)[colnames(sumstats) == rsID] <- rsID_new
  }
  if (!is.null(CHR)) {
    colnames(sumstats)[colnames(sumstats) == CHR] <- CHR_new
  }
  if (!is.null(POS)) {
    if (class(sumstats[[POS]]) != "integer") {
      warning("POS column is not of the 'integer' class, converting to integer")
      sumstats[[POS]] <- as.integer(sumstats[[POS]])
    }
    colnames(sumstats)[colnames(sumstats) == POS] <- POS_new
  }
  if (!is.null(A2)) {
    colnames(sumstats)[colnames(sumstats) == A2] <- A2_new
  }
  if (!is.null(AF)) {
    colnames(sumstats)[colnames(sumstats) == AF] <- AF_new
  }
  if (!is.null(N)) {
    colnames(sumstats)[colnames(sumstats) == N] <- N_new
  }
  return(sumstats)
}
