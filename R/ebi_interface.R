# Ontology Lookup Service ----
# 1. Function to search EFO terms from OLS
search_efo_terms <- function(query) {
  base_url <- "https://www.ebi.ac.uk/ols/api/search"
  
  response <- GET(base_url, 
                  query = list(
                    q = query,
                    ontology = "efo",
                    rows = 50
                  ))
  
  if (status_code(response) == 200) {
    content <- fromJSON(rawToChar(response$content))
    
    if (length(content$response$docs) > 0) {
      results <- data.frame(
        label = content$response$docs$label,
        iri = content$response$docs$iri,
        short_form = content$response$docs$short_form,
        stringsAsFactors = FALSE
      )
      return(results[grepl("^EFO_", results$short_form), ])
    }
  }
  return(NULL)
}

# GWAS Catalog ----
# 2. Function to get associations for an EFO term
get_gwas_associations <- function(efo_id, efo_id_label) {
  
  message("Querying ", efo_id, ",", efo_id_label)
  base_url <- paste0("https://www.ebi.ac.uk/gwas/rest/api/efoTraits/", efo_id, "/associations")
  
  response <- GET(base_url)
  
  if (status_code(response) == 200) {
    content <- fromJSON(rawToChar(response$content))
    
    if (length(content$`_embedded`$associations) > 0) {
      associations <- content$`_embedded`$associations
      
      results <- list()
      
      for(i in seq_along(associations$loci)) {
        # Extract variant information
        risk_allele <- associations$loci[[i]]$strongestRiskAlleles[[1]]$riskAlleleName[1]
        
        # Get location based on whether it's rsID or chr:pos
        location <- if(grepl("rs", risk_allele)) {
          get_snp_details(risk_allele)
        } else {
          gsub("-\\?", "", risk_allele)  # Clean up chr:pos format
        }
        location <- gsub("chr", "", location)
        
        # Format p-value in scientific notation
        pvalue <- as.numeric(sprintf("%.0fe%d", 
                                     associations$pvalueMantissa[i],
                                     associations$pvalueExponent[i]))
        
        # Format beta with +/- sign
        beta_formatted <- as.numeric(if(!is.na(associations$betaNum[i])) {
          direction_sign <- if(associations$betaDirection[i] == "decrease") "-" else "+"
          paste0(direction_sign, sprintf("%.2f", associations$betaNum[i]))
        } else {
          "-"
        })
        
        # Create row with simplified columns
        row <- data.frame(
          variant_risk_allele = risk_allele,
          location = location,
          pvalue = pvalue,
          raf = ifelse(is.na(associations$riskFrequency[i]), "NR", 
                       associations$riskFrequency[i]),
          beta = beta_formatted,
          ci = ifelse(is.na(associations$range[i]), "-", associations$range[i]),
          mapped_gene = ifelse(length(associations$loci[[i]]$authorReportedGenes[[1]]) == 0,
                               "-", 
                               paste(unlist(associations$loci[[i]]$authorReportedGenes), 
                                     collapse = ", ")),
          efo_id = efo_id,
          efo_name = efo_id_label,
          stringsAsFactors = FALSE
        )
        
        results[[i]] <- row
      }
      
      return(do.call(rbind, results))
    }
  } else {
    data.frame(efo_id = efo_id,
               efo_name = efo_id_label)
  }
}

# 3. Helper function to get SNP details
get_snp_details <- function(rsid) {
  rsid <- strsplit(rsid, "-")[[1]][1]
  snp_url <- paste0("https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/", rsid)
  
  response <- GET(snp_url)
  if (status_code(response) == 200) {
    content <- fromJSON(rawToChar(response$content))
    return(paste(content$locations$chromosomeName, content$locations$chromosomePosition, sep = ":"))
  }
  return(NA)
}

#' download_GCST
#' @importFrom RCurl getURL
download_GCST <- function(url_in, dest_folder) {
  if (missing(dest_folder)) stop("dest_folder not provided")
  if (!dir.exists(dest_folder)) dir.create(dest_folder)
  url_in <- gsub("^http://", "ftp://", url_in)
  url_in <- paste0(url_in, "/")
  tryCatch({
    files <- RCurl::getURL(url_in, ftp.use.epsv = FALSE, dirlistonly = TRUE)
    file_list <- strsplit(files, "\n|\r")[[1]]
    if ("harmonised" %in% file_list) {
      url_in <- paste0(url_in, "harmonised/")
      files <- RCurl::getURL(url_in, ftp.use.epsv = FALSE, dirlistonly = TRUE)
      file_list <- strsplit(files, "\n|\r")[[1]]
    }
    lapply(file_list, function(dest_file) {
      if (!file.exists(file.path(dest_folder, dest_file))) {
        system(paste0("wget ", paste0(url_in, "/", dest_file), " -O ", file.path(dest_folder, dest_file)))
      } else {
        message("dest_file ", dest_file, " file found, not rewriting")
      }
    })
  }, error = function(e) {
    message("Error accessing ", url_in, ": ", e$message)
  })
}


#' process_GCST
#' @importFrom data.table fread
process_GCST <- function(file_in, test_mode=F) {
  
  if (test_mode) {
    sumstats <- fread(cmd=paste0("zcat ", file_in, " | head -2"))
  } else {
    sumstats <- fread(file_in)
  }
  
  # find pval column
  p_value_col <- "pval"
  if (! "pval" %in% colnames(sumstats)) {
    if ("p_value" %in% colnames(sumstats)) {
      p_value_col <- "p_value"
    } else {
      stop("p-value column not found")
    }
  }
  
  # find effect size column
  es_col <- "beta"
  if (! "beta" %in% colnames(sumstats)) {
    if (all(c("odds_ratio", "ci_lower", "ci_upper") %in% colnames(sumstats))) {
      # calculate standard error from OR and CI
      sumstats[[es_col]] <- log(sumstats$odds_ratio)
      # deal with a special case when all(ci_lower == T), e.g., GCST90476253.h
      if (all(sumstats$ci_lower == T & sumstats$ci_upper == T) & p_value_col %in% colnames(sumstats)) {
        z_score <- qnorm(sumstats$p_value/2, lower.tail = FALSE)
        sumstats$standard_error <- abs(sumstats[[es_col]]) / z_score
      } else {
        sumstats$standard_error <- (log(sumstats$ci_upper) - log(sumstats$ci_lower)) / (2 * 1.96)
      }
    } else {
      stop("effect size column not found")
    }
  }
  
  # find sample size
  N_col <- "n"
  if ( ! N_col %in% colnames(sumstats)) {
    if ("num_cases" %in% colnames(sumstats) & "num_controls" %in% colnames(sumstats)) {
      sumstats[[N_col]] <- sumstats$num_cases + sumstats$num_controls
    } else if (file.exists(paste0(file_in, "-meta.yaml"))) {
      ss_txt <- system(paste0("grep sample_size ", paste0(file_in, "-meta.yaml")), intern = T)
      n_tmp <- sum(unique(as.numeric(gsub("[a-z_:\ ]", "", ss_txt))))
      sumstats[[N_col]] <- n_tmp
      message("Getting sample size from the yaml file, please check manually: ",
              n_tmp, "; ", file_in)
    } else {
      stop("sample size column not found")
    }
  }
  
  # reorder and output
  cols <- c("chromosome", "base_pair_location", "effect_allele", "other_allele",
            es_col, "standard_error", p_value_col, "effect_allele_frequency", N_col)
  stopifnot(cols %in% colnames(sumstats))
  sumstats <- sumstats[,..cols]
  return(sumstats)
}

# eQTL Catalogue ----

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

#' annotate_eQTL_Catalog
#' @export
annotate_eQTL_Catalog <- function(coloc_out, annotation_df,
                                  datasets_eQTL_Catalogue) {
  coloc_out$dataset_id <- gsub("(.*)_.*", "\\1", coloc_out$sumstats_2_file)
  coloc_out$gene_id <- gsub(".*_(.*)", "\\1", coloc_out$sumstats_2_file)
  coloc_out <- merge(coloc_out, annotation_df, all.x=T)[, union(names(coloc_out), names(annotation_df))]
  coloc_out <- merge(coloc_out, datasets_eQTL_Catalogue, all.x=T)[, union(names(coloc_out), names(datasets_eQTL_Catalogue))]
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

#' make_eQTL_Catalogue_args
#' @export
make_eQTL_Catalogue_args <- function(QTS="QTS000015") {
  datasets_eQTL_Catalogue <- get_datasets_eQTL_Catalogue()
  datasets_eQTL_Catalogue <- subset(datasets_eQTL_Catalogue,
                                    quant_method == "ge" & study_id %in% QTS)
  eQTL_Catalogue <- list()
  eQTL_Catalogue$sumstats_2_files <- datasets_eQTL_Catalogue$dataset_id
  eQTL_Catalogue$sumstats_2_function <- "query_eQTL_Catalogue"
  eQTL_Catalogue$sumstats_2_type <- "quant"
  eQTL_Catalogue$sumstats_2_sdY <- NA
  eQTL_Catalogue
}



