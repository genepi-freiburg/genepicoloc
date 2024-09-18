#' query_sumstats_1
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

#' query_eQTL_Catalogue
#' @importFrom httr GET
#' @importFrom httr accept_json
#' @importFrom httr content
#' @importFrom httr status_code
#' @importFrom jsonlite fromJSON
#' @export
query_eQTL_Catalogue <- function(sumstats_file,
                                 CHR_var, BP_START_var, BP_STOP_var,
                                 start=0, size=1000) {
  pos <- paste0(CHR_var, ":", BP_START_var, "-", BP_STOP_var)
  message("Querying eQTL Catalog.", " File: ", sumstats_file)
  while (T) {
    sumstats_URL <- paste0("https://www.ebi.ac.uk/eqtl/api/v2/datasets/",
                           sumstats_file, "/associations?size=",
                           size, "&start=", start, "&",
                           "pos=", pos)
    message("Iterations: ", start, "-", start+size, "... ", appendLF = F)
    r <- GET(sumstats_URL, accept_json())
    cont <- content(r, "text", encoding = "UTF-8")
    Sys.sleep(2) # to avoid too many API requests
    if (status_code(r) == 200) {
      message("Success, continuing queries.")
    } else {
      if (start == 0) {
        message("No lines were available, returning empty data.frame.")
        return(data.frame())
        stop(paste0("Error ", status_code(r)))
      }
      message("Reached the end of the table.")
      message("All lines were successfully queried.")
      break
    }
    sumstats_chunk <- fromJSON(cont)
    if (start == 0) {
      sumstats <- sumstats_chunk
    } else {
      sumstats <- rbind(sumstats, sumstats_chunk)
    }
    # add counter
    start <- start + size
  }
  # format
  sumstats <- format_eQTL_Catalogue(sumstats)
  return(sumstats)
}

#' get_datasets_eQTL_Catalogue
#' @importFrom httr GET
#' @importFrom httr accept_json
#' @importFrom httr content
#' @importFrom httr status_code
#' @importFrom jsonlite fromJSON
#' @export
get_datasets_eQTL_Catalogue <- function() {
  message("Querying datasets_eQTL_Catalogue table.")
  datasets_url <- "https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size=1000&start=0"
  r <- GET(datasets_url, accept_json())
  cont <- content(r, "text", encoding = "UTF-8")
  datasets_eQTL_Catalogue <- fromJSON(cont)
  return(datasets_eQTL_Catalogue)
}

# helpers
#' format_eQTL_Cat
format_eQTL_Catalogue <- function(sumstats) {
  sumstats$Name <- paste0("chr", sumstats$chromosome, ":", sumstats$position, ":",
                          sumstats$ref, ":", sumstats$alt)
  cols <- c("Name", "rsid", "chromosome", "position", "alt", "ref", "beta", "se", "nlog10p", "maf", "an", "molecular_trait_id")
  sumstats <- sumstats[,cols]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "nlog10P", "AF", "N", "Phenotype")
  # format by phenotype ID
  sumstats <- sapply(unique(sumstats$Phenotype), function(x) {
    sumstats <- subset(sumstats, Phenotype == x)
  }, simplify = F)
  return(sumstats)
}

