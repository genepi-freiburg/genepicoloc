library(data.table)

mvp_file <- "~/Work/bioinfo/nextcloud/Downloads/MVP_R4.1000G_AGR.BUN_BSP_Mean_INT.EUR.GIA.dbGaP.txt.gz"
old_pass  <- "~/Work/bioinfo/nextcloud/Downloads/BUN_BSP_Mean_INT_EUR_coloc_regions_PASS.tsv"
old_all   <- "~/Work/bioinfo/nextcloud/Downloads/BUN_BSP_Mean_INT_EUR_coloc_regions.tsv"

# Load old results
old_regions <- fread(old_pass)
message("Old PASS regions: ", nrow(old_regions))

# Load and format sumstats
message("Loading BUN MVP file...")
t0 <- proc.time()
dt <- fread(mvp_file, select = c("chrom", "pos", "pval"))
message("Loaded ", nrow(dt), " rows in ", round((proc.time() - t0)[3], 1), "s")
dt[, CHR := chrom][, POS := pos]
dt[, nlog10P := -log10(pval)]
dt[!is.finite(nlog10P), nlog10P := 350]
ss <- as.data.frame(dt[, .(CHR, POS, nlog10P)])
rm(dt); gc(verbose = FALSE)

message("N significant: ", sum(ss$nlog10P > 7.30103, na.rm = TRUE))

# Run new version
source("R/get_coloc_regions.R")
message("\nRunning new get_coloc_regions...")
t1 <- proc.time()
res <- get_coloc_regions(ss, CHR_name = "CHR", POS_name = "POS",
                         nlog10p_value_name = "nlog10P")
t_new <- (proc.time() - t1)[3]
new_regions <- as.data.table(res$coloc_regions_PASS)
message("New PASS regions: ", nrow(new_regions), " in ", round(t_new, 2), "s")

# Compare
message("\n=== COMPARISON ===")
key_old <- paste(old_regions$CHR_var, old_regions$BP_START_var, old_regions$BP_STOP_var, sep = ":")
key_new <- paste(new_regions$CHR_var, new_regions$BP_START_var, new_regions$BP_STOP_var, sep = ":")

n_exact <- sum(key_new %in% key_old)
message("Exact matches: ", n_exact, " / ", nrow(new_regions))
message("Only in NEW: ", length(setdiff(key_new, key_old)))
message("Only in OLD: ", length(setdiff(key_old, key_new)))

# Check chromosome-level match (same chr, overlapping intervals)
setnames(old_regions, c("CHR_var","BP_START_var","BP_STOP_var"), c("chr","start_old","stop_old"))
setnames(new_regions, c("CHR_var","BP_START_var","BP_STOP_var"), c("chr","start_new","stop_new"))
merged <- merge(old_regions[, .(chr, start_old, stop_old, id_old = .I)],
                new_regions[, .(chr, start_new, stop_new, id_new = .I)],
                by = "chr", allow.cartesian = TRUE)
merged[, overlap := pmin(stop_old, stop_new) - pmax(start_old, start_new)]
matched <- merged[overlap > 0, .SD[which.max(overlap)], by = id_new]

message("\nNew regions with overlapping old region: ", nrow(matched), " / ", nrow(new_regions))
boundary_diffs <- matched[start_old != start_new | stop_old != stop_new]
message("Boundary differences (merge-related): ", nrow(boundary_diffs))

if (nrow(boundary_diffs) > 0) {
  message("\nFirst 10 boundary differences:")
  print(head(boundary_diffs[, .(chr, start_old, start_new, stop_old, stop_new,
                                 d_start = start_new - start_old,
                                 d_stop = stop_new - stop_old)], 10))
}

# Unmatched
unmatched_new <- setdiff(seq_len(nrow(new_regions)), matched$id_new)
unmatched_old <- setdiff(seq_len(nrow(old_regions)), matched$id_old)
if (length(unmatched_new) > 0) {
  message("\nNew regions without old match:")
  print(new_regions[unmatched_new])
}
if (length(unmatched_old) > 0) {
  message("\nOld regions without new match:")
  print(old_regions[unmatched_old])
}
