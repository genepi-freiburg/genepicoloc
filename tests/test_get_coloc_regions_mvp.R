library(data.table)

# --- Load test data ---
# Usage: Rscript tests/test_get_coloc_regions_mvp.R [chr21|chr22|...|full]
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) > 0) args[1] else "chr21"

mvp_file <- "~/Work/bioinfo/nextcloud/Downloads/MVP_R4.1000G_AGR.Creat_BSP_Mean_INT.EUR.GIA.dbGaP.txt.gz"

if (grepl("^chr", mode)) {
  chr_num <- as.integer(sub("^chr", "", mode))
  cache_file <- paste0("tests/mvp_chr", chr_num, ".rds")
  if (!file.exists(cache_file)) {
    message("Extracting chr", chr_num, " from full MVP file...")
    dt <- fread(mvp_file, select = c("chrom", "pos", "pval"))
    dt <- dt[chrom == chr_num]
    dt[, CHR := chrom][, POS := pos]
    dt[, nlog10P := -log10(pval)]
    dt[!is.finite(nlog10P), nlog10P := 350]
    ss <- as.data.frame(dt[, .(CHR, POS, nlog10P)])
    saveRDS(ss, cache_file)
    message("Saved ", nrow(ss), " chr", chr_num, " rows to ", cache_file)
    rm(dt)
  } else {
    message("Loading cached chr", chr_num, " data...")
    ss <- readRDS(cache_file)
  }
} else {
  message("Loading FULL MVP file...")
  dt <- fread(mvp_file, select = c("chrom", "pos", "pval"))
  dt[, CHR := chrom][, POS := pos]
  dt[, nlog10P := -log10(pval)]
  dt[!is.finite(nlog10P), nlog10P := 350]
  ss <- as.data.frame(dt[, .(CHR, POS, nlog10P)])
  rm(dt)
}

message("Rows: ", nrow(ss))
message("Max -log10P: ", round(max(ss$nlog10P, na.rm = TRUE), 2))
message("N significant (p < 5e-8): ", sum(ss$nlog10P > 7.30103, na.rm = TRUE))
gc(verbose = FALSE)

# Source the NEW version
source("R/get_coloc_regions.R")
get_coloc_regions_new <- get_coloc_regions

# Source the OLD version from git
old_code <- system("git show HEAD:R/get_coloc_regions.R", intern = TRUE)
old_code <- sub("^get_coloc_regions <-", "get_coloc_regions_old <-", old_code)
tmp <- tempfile(fileext = ".R")
writeLines(old_code, tmp)
source(tmp)

message("\n=== NEW (sorted-peak) ===")
t1 <- proc.time()
res_new <- get_coloc_regions_new(ss, CHR_name = "CHR", POS_name = "POS",
                                 nlog10p_value_name = "nlog10P")
t_new <- (proc.time() - t1)[3]
n_new <- if (is.data.frame(res_new$coloc_regions_PASS)) nrow(res_new$coloc_regions_PASS) else 0
message("NEW: ", n_new, " PASS regions in ", round(t_new, 2), "s")

message("\n=== OLD (original) ===")
t2 <- proc.time()
res_old <- get_coloc_regions_old(ss, CHR_name = "CHR", POS_name = "POS",
                                 nlog10p_value_name = "nlog10P")
t_old <- (proc.time() - t2)[3]
n_old <- if (is.data.frame(res_old$coloc_regions_PASS)) nrow(res_old$coloc_regions_PASS) else 0
message("OLD: ", n_old, " PASS regions in ", round(t_old, 2), "s")

message("\n=== COMPARISON ===")
message("Speedup: ", round(t_old / t_new, 1), "x")
message("N PASS regions - new: ", n_new, ", old: ", n_old, ", match: ", n_new == n_old)

# Compare PASS region coordinates (sorted by chr, start)
pass_new <- res_new$coloc_regions_PASS
pass_old <- res_old$coloc_regions_PASS
pass_new <- pass_new[order(pass_new$CHR_var, pass_new$BP_START_var), ]
pass_old <- pass_old[order(pass_old$CHR_var, pass_old$BP_START_var), ]
rownames(pass_new) <- rownames(pass_old) <- NULL
message("PASS coordinates identical: ", identical(pass_new, pass_old))

if (!identical(pass_new, pass_old)) {
  # Check if the sets of intervals match even if numbering differs
  key_new <- paste(pass_new$CHR_var, pass_new$BP_START_var, pass_new$BP_STOP_var, sep = ":")
  key_old <- paste(pass_old$CHR_var, pass_old$BP_START_var, pass_old$BP_STOP_var, sep = ":")
  message("Same set of intervals: ", setequal(key_new, key_old))
  only_new <- setdiff(key_new, key_old)
  only_old <- setdiff(key_old, key_new)
  if (length(only_new) > 0) message("Only in NEW (", length(only_new), "): ", paste(head(only_new, 5), collapse = ", "))
  if (length(only_old) > 0) message("Only in OLD (", length(only_old), "): ", paste(head(only_old, 5), collapse = ", "))
}

# Compare filtered sumstats
nfilt_new <- nrow(res_new$sumstats_filt)
nfilt_old <- nrow(res_old$sumstats_filt)
message("\nFiltered rows - new: ", nfilt_new, ", old: ", nfilt_old, ", match: ", nfilt_new == nfilt_old)
