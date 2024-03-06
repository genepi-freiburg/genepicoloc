#' add_cumsum
#' Add columns with the cumulative sum of BP position for plotting
#' @param coloc_summary Output from genepicoloc pipeline (summary.RDS)
#' @export
add_cumsum <- function(coloc_summary,
                       position="position",
                       chr="CHR_var") {
  coloc_summary <- readRDS("summary.RDS")
  position <- "position"
  chr <- "CHR_var"
  if (any(coloc_summary[[chr]] == "X")) {
    warning("X chromosome coded as 'X', changing values to '23' and class to numeric")
    coloc_summary[[chr]][coloc_summary[[chr]] == "X"] <- "23"
    coloc_summary[[chr]] <- as.numeric(coloc_summary[[chr]])
    coloc_summary <- coloc_summary[order(coloc_summary[[chr]]),]
  }
  if (class(coloc_summary[[chr]]) == "character") {
    warning("Chromosome column class is character, changing to numeric")
    coloc_summary[[chr]] <- as.numeric(coloc_summary[[chr]])
    coloc_summary <- coloc_summary[order(coloc_summary[[chr]]),]
  }
  
  coloc_summary[[position]] <- (coloc_summary[["BP_START_var"]]+coloc_summary[["BP_STOP_var"]])/2
  coloc_summary_cumsum <- tapply(coloc_summary[[position]], coloc_summary[[chr]], max)
  coloc_summary_cumsum <- data.frame(CHR = as.numeric(names(coloc_summary_cumsum)),
                                     BPmax = coloc_summary_cumsum)
  colnames(coloc_summary_cumsum)[1] <- chr
  cumsum_vector <- c(0,cumsum(as.numeric(coloc_summary_cumsum$BPmax)))
  coloc_summary_cumsum$BP_cumsum <- cumsum_vector[-length(cumsum_vector)]
  coloc_summary_cumsum$BPmax <- NULL
  coloc_summary_merged <- merge(coloc_summary, coloc_summary_cumsum, by=chr)
  coloc_summary_merged$BP_cumsum_updated <- coloc_summary_merged[[position]] + coloc_summary_merged$BP_cumsum
  coloc_summary_merged$BP_cumsum <- NULL
  
  axis_x <- tapply(coloc_summary_merged$BP_cumsum_updated, coloc_summary_merged[[chr]], median)
  axis_x <- data.frame(CHR = names(axis_x),
                       BP_cumsum_updated_median = axis_x)
  coloc_summary_merged[[chr]] <- as.factor(coloc_summary_merged[[chr]])
  return(coloc_summary_merged)
}

#' plot_FinnGen
#' Visualize colocalization results with FinnGen study
#' @param coloc_summary_merged Output from add_cumsum function
#' @export
plot_FinnGen <- function(coloc_summary_merged,
                         FinnGen_r9_var="FinnGen_r9",
                         chr="CHR_var",
                         max_char=30) {
  if (!"ggplot2" %in% rownames(installed.packages())) {
    stop("'ggplot2' is required to run this function'")
  }
  if (!"RColorBrewer" %in% rownames(installed.packages())) {
    stop("'RColorBrewer' is required to run this function'")
  }
  if (!"ggrepel" %in% rownames(installed.packages())) {
    stop("'ggrepel' is required to run this function'")
  }
  library(ggplot2)
  library(RColorBrewer)
  library(ggrepel)
  # subset
  FinnGen_r9 <- coloc_summary_merged[coloc_summary_merged[["Dataset"]] == FinnGen_r9_var,]
  FinnGen_r9 <- subset(FinnGen_r9, grepl("^I.*|^X.*|^V.*", FinnGen_r9_category))
  FinnGen_r9 <- subset(FinnGen_r9, grepl("^I.*|^X.*|^V.*", FinnGen_r9_category))
  FinnGen_r9 <- subset(FinnGen_r9, FinnGen_r9_category %in% names(table(FinnGen_r9$FinnGen_r9_category)[table(FinnGen_r9$FinnGen_r9_category) > 5]))
  # TODO fix
  colnames(FinnGen_r9)[colnames(FinnGen_r9) == "FinnGen_r9_X...phenotype"] <- "FinnGen_r9_phenotype"
  # cut too long phenotype names
  to_update <- sapply(FinnGen_r9[["FinnGen_r9_phenotype"]], nchar) > max_char
  FinnGen_r9[["FinnGen_r9_phenotype"]][to_update] <- substr(FinnGen_r9[["FinnGen_r9_phenotype"]][to_update], 1, max_char-3)
  FinnGen_r9[["FinnGen_r9_phenotype"]][to_update] <- paste0(FinnGen_r9[["FinnGen_r9_phenotype"]][to_update], "...")
  # Update manually
  to_update <- FinnGen_r9[["FinnGen_r9_category"]] == "III Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism (D3_)"
  FinnGen_r9[["FinnGen_r9_category"]][to_update] <- FinnGen_r9[["FinnGen_r9_category"]][to_update] <- "III Diseases of the blood and blood-forming organs"
  # gg object
  FinnGen_r9_gg <- ggplot(FinnGen_r9, aes_string(x="BP_cumsum_updated",
                                                 y="FinnGen_r9_category",
                                                 size="sumstats_2_max_nlog10P")) +
    geom_point(aes_string(color=chr)) +
    scale_color_manual(values = rep(brewer.pal(8, "Dark2"), 3)) +
    scale_x_continuous(expand = c(0.02, 0.02),
                       label = axis_x[["CHR"]],
                       breaks= axis_x[["BP_cumsum_updated_median"]]) +
    labs(x = "Chromosome", y = "Disease categories") +
    ggtitle("FinnGen study (Phenome-wide association study)") +
    theme_bw() +
    theme(
      legend.position="none"
    ) +
    geom_text_repel(data = subset(FinnGen_r9, sumstats_2_max_nlog10P > 30), aes(label=FinnGen_r9_phenotype), size = 4)
  return(FinnGen_r9_gg)
}
