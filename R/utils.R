# Tools and utils to facilitate other genepicoloc processes

#' handle_underflow
#' Process p-values to handle underflow in case of very small p-values (p<1e-320)
#' @importFrom Rmpfr mpfr
#' @export
handle_underflow <- function(pvalue_vec,
                             return_nlog10P=F) {
  message("Converting to numeric using Rmpfr::mpfr() ... ")
  if (is.numeric(pvalue_vec)) {
    pvalue_vec <- Rmpfr::mpfr(pvalue_vec, precBits=100)
  } else {
    pvalue_vec <- Rmpfr::mpfr(pvalue_vec)
  }
  if (return_nlog10P) {
    message("Converting p-values to negative log10 scale")
    pvalue_vec <- as.numeric(-log10(pvalue_vec))
  }
  return(pvalue_vec)
}

#' flip_alleles
#' flip alleles according to complementary base pairing 
flip_alleles <- function(vec) {
  vec_out <- vec
  vec_out <- toupper(vec_out)
  vec_out[vec == "A"] <- "T"
  vec_out[vec == "T"] <- "A"
  vec_out[vec == "C"] <- "G"
  vec_out[vec == "G"] <- "C"
  return(vec_out)
}



simulate_sumstats <- function(seed=12345, POS=NULL, A1=NULL, A2=NULL,
                              CHR_var="16", BP_START_var=19881010,
                              BP_STOP_var=20881010, N_var=5000,
                              rsID=NA,
                              nucleotides=c("A", "T", "G", "C"), N=10000,
                              BETA_min = -1, BETA_max=1, SE_min=0.1, SE_max=1,
                              AF_min=0.01, AF_max=0.99) {
  set.seed(seed)
  if (is.null(POS)) {
    POS <- sort(sample(BP_START_var:BP_STOP_var, N_var))
  } else { N_var <- length(POS) }
  if (is.null(A1)) A1 <- sample(nucleotides, N_var, replace = T)
  if (is.null(A2)) A2 <- sample(nucleotides, N_var, replace = T)
  A2[A2 == A1] <- sapply(A2[A2 == A1], function(x) nucleotides[x != nucleotides][1])
  stopifnot(all(A2 != A1))
  Name <- paste0("chr", CHR_var, ":", POS, ":", A2, ":", A1)
  BETA <- runif(N_var, min = BETA_min, max=BETA_max)
  SE <- runif(N_var, min = SE_min, max=SE_max)
  P <- pnorm(BETA/SE, lower.tail=F)*2
  nlog10P <- -log10(P)
  AF <- runif(N_var, min = AF_min, max=AF_max)
  data.frame(Name=Name, rsID=rsID, CHR=CHR_var, POS=POS, #A1=A1, A2=A2,
             BETA=BETA, SE=SE, P=P, nlog10P=nlog10P, AF=AF, N=N)
}

write_simulated_sumstats <- function(N=1, seed=12345, POS=NULL, A1=NULL, A2=NULL) {
  tmp_file <- tempfile()
  # simulate sumstats, save, index
  sumstats1 <- simulate_sumstats(seed=seed, POS=POS, A1=A1, A2=A2)
  write.table(sumstats1, paste0(tmp_file, "_ss", seed, ".tsv"), sep="\t", row.names = F, quote=F)
  system(paste0("bgzip ", tmp_file, "_ss", seed, ".tsv"))
  system(paste0("tabix -s3 -b4 -e4 ", tmp_file, "_ss", seed, ".tsv.gz -c Name"))
  if (N>1) {
    sumstats_N <- sapply(2:N, function(i) {
      seed <- seed+(i-1)
      sumstats <- simulate_sumstats(seed=seed, POS=POS, A1=A1, A2=A2)
      sumstats$Name <- sumstats1$Name
      sumstats$A1 <- sumstats1$A1
      sumstats$A2 <- sumstats1$A2
      sumstats$AF <- sumstats1$AF
      write.table(sumstats, paste0(tmp_file, "_ss", seed, ".tsv"), sep="\t", row.names = F, quote=F)
      system(paste0("bgzip -f ", tmp_file, "_ss", seed, ".tsv"))
      system(paste0("tabix -f -s3 -b4 -e4 ", tmp_file, "_ss", seed, ".tsv.gz -c Name"))
    }, simplify = F)
  }
  return(tmp_file)
}

unlink_simulated_sumstats <- function(tmp_file) {
  lapply(c(paste0(tmp_file, "_ss1.tsv.gz"), paste0(tmp_file, "_ss1.tsv.gz.tbi"), 
           paste0(tmp_file, "_ss2.tsv.gz"), paste0(tmp_file, "_ss2.tsv.gz.tbi")),
         unlink)
  return(NULL)
}


#'@export
retrieve_sumstats <- function(query_fun, sumstats_file,
                              CHR_var, BP_START_var, BP_STOP_var,
                              gene=NULL, nlog10P=NULL) {
  # retrieve sumstats
  sumstats <- do.call(query_fun, list(sumstats_file=sumstats_file,
                                      CHR_var=CHR_var,
                                      BP_START_var=BP_START_var,
                                      BP_STOP_var=BP_STOP_var))
  # if reading eQTL sumstats, select the required gene from the list
  if (!is.null(gene)) {
    sumstats <- list_to_df_eQTL(sumstats, gene)
  }
  # create "P" column
  ## When reading sumstats with "nlog10P" column, convert to "P" using mpfr
  if (!missing(nlog10P)) {
    if (nlog10P) {
      sumstats[["P"]] <- 10^-sumstats[["nlog10P"]]
    } else {
      ## When reading sumstats with character "P" column, convert to numeric "P" using mpfr
      if (is.character(sumstats[["P"]])) {
        sumstats[["nlog10P"]] <- -log10(Rmpfr::mpfr(sumstats[["P"]]))
      } else {
        sumstats[["nlog10P"]] <- -log10(sumstats[["P"]])
      }
    }
  }
  return(sumstats)
}

#' list_to_df_eQTL
#' processing eQTL datasets (more than 1 gene in the same sumstats file)
list_to_df_eQTL <- function(sumstats, gene,
                            Phenotype_col="Phenotype") {
  if (is.data.frame(sumstats)) {stop("input sumstats is data.frame but it should be list")}
  ind <- which(sapply(sumstats, function(x) {
    all(grepl(gene, x[[Phenotype_col]]))
  }))
  if (length(ind) == 0) {
    stop("No genes were matched, please check manually")
  } else if (length(ind) > 1) {
    stop("More than one gene was matched, please check manually")
  } else {
    sumstats <- sumstats[[ind]]
  }
}


