#' manhattan.plot
manhattan.plot <- function(sumstats, manhattan_title = NA, manhattan_subtitle = NA,
                           position = "position", chr = "chr", pval = "pval",
                           SNP = "SNP", pval_filtering, pval_threshold,
                           nlog10 = F, annotate_all = F) {
  # inspired by RColorBrewer::brewer.pal(8, "Dark2"), not importing to limit dependencies:
  favorite_colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                       "#66A61E", "#E6AB02", "#A6761D", "#666666")
  if (any(sumstats[[chr]] == "X")) {
    warning("X chromosome coded as 'X', changing values to '23' and class to numeric")
    sumstats[[chr]][sumstats[[chr]] == "X"] <- "23"
    sumstats[[chr]] <- as.numeric(sumstats[[chr]])
    sumstats <- sumstats[order(sumstats[[chr]]),]
  }
  if (class(sumstats[[chr]]) == "character") {
    warning("Chromosome column class is character, changing to numeric")
    sumstats[[chr]] <- as.numeric(sumstats[[chr]])
    sumstats <- sumstats[order(sumstats[[chr]]),]
  }
  if (nlog10 == F) {
    sumstats <- sumstats[sumstats[[pval]] < pval_filtering,]
  } else {
    sumstats <- sumstats[sumstats[[pval]] >= pval_filtering,]
  }
  selected_cols <- c(position, chr, pval, SNP)
  sumstats <- sumstats[,selected_cols]
  sumstats_cumsum <- tapply(sumstats[[position]], sumstats[[chr]], max)
  sumstats_cumsum <- data.frame(CHR = as.numeric(names(sumstats_cumsum)),
                                BPmax = sumstats_cumsum)
  colnames(sumstats_cumsum)[1] <- chr
  cumsum_vector <- c(0,cumsum(as.numeric(sumstats_cumsum$BPmax)))
  sumstats_cumsum$BP_cumsum <- cumsum_vector[-length(cumsum_vector)]
  sumstats_cumsum$BPmax <- NULL
  sumstats_merged <- merge(sumstats, sumstats_cumsum, by=chr)
  sumstats_merged$BP_cumsum_updated <- sumstats_merged[[position]] + sumstats_merged$BP_cumsum
  sumstats_merged$BP_cumsum <- NULL
  # Find top SNP for each chromosome
  if (nlog10 == F) {
    sumstats_merged_sign <- sumstats_merged[sumstats_merged[[pval]] < pval_threshold,]
  } else {
    sumstats_merged_sign <- sumstats_merged[sumstats_merged[[pval]] >= pval_threshold,]
  }
  if (nrow(sumstats_merged_sign) > 0) {
    sumstats_merged_sign[[chr]] <- as.factor(sumstats_merged_sign[[chr]])
    top_snps <- by(sumstats_merged_sign, sumstats_merged_sign[[chr]], function(sub_df) {
      if (nlog10 == F) {
        sub_df[which.min(sub_df[[pval]]),]
      } else {
        sub_df[which.max(sub_df[[pval]]),]
      }
    })
    top_snps <- do.call(rbind, top_snps)
    if (nlog10 == F) {
      top_snps[[pval]] <- -log10(top_snps[[pval]])
    }
    if (annotate_all) {
      top_snps <- sumstats_merged_sign
    }
  }
  # create x axis
  axis_x <- tapply(sumstats_merged$BP_cumsum_updated, sumstats_merged[[chr]], median)
  axis_x <- data.frame(CHR = names(axis_x),
                       BP_cumsum_updated_median = axis_x)
  # to discrete
  sumstats_merged[[chr]] <- as.factor(sumstats_merged[[chr]])
  if (nlog10 == F) {
    sumstats_merged[[pval]] <- -log10(sumstats_merged[[pval]])
  }
  # plot
  manh_gg <- ggplot2::ggplot(sumstats_merged, aes_string(x="BP_cumsum_updated", y=pval)) +
    geom_point(aes_string(color=chr), size=1) + # alpha=0.8, 
    scale_color_manual(values = rep(favorite_colors, 3)) +
    scale_x_continuous(expand = c(0.02, 0.02),
                       label = axis_x$CHR,
                       breaks= axis_x$BP_cumsum_updated_median) +
    scale_y_continuous(expand = c(0.05, 0.05)) +
    labs(x = "Chromosome", y = "-log10(P)") +
    geom_hline(yintercept = -log10(5e-8)) +
    geom_hline(yintercept = -log10(1e-5), linetype="dashed") +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.subtitle = element_text(hjust = 0.5)
    )
  if (nrow(sumstats_merged_sign) > 0) {
    manh_gg <- manh_gg +
      geom_text(data = top_snps, aes_string(x = "BP_cumsum_updated", y = pval, label = SNP), size = 2, vjust = -1)
    # disable to limit dependencies
    # manh_gg <- manh_gg + geom_label_repel(data = top_snps, aes_string(x = "BP_cumsum_updated", y = pval, label = SNP),
    #                                       min.segment.length = unit(0, 'lines'))
  }
  if (!is.na(manhattan_title)) {
    manh_gg <- manh_gg + ggtitle(label=manhattan_title)
  }
  if (!is.na(manhattan_title) & !is.na(manhattan_subtitle)) {
    manh_gg <- manh_gg + ggtitle(label=manhattan_title, subtitle=manhattan_subtitle)
  }
  return(manh_gg)
}
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

locus_input <- function(sumstats, index_var=NULL,
                        plink2="plink2", bfile=NULL,
                        do_locuszoom_input=T, add_r_mat=NULL,
                        rm_tmp=T) {
  if (!is.null(add_r_mat)) {stop("Please use LD_matrix() function to retrieve LD matrix")}
  tmp <- tempfile()
  write.table(sumstats[,"Name"], tmp, row.names = F, col.names = F, quote=F)
  # bfile for plink2 LD estimation
  if (is.null(bfile)) bfile <- get_bfile(CHR=CHR_var)
  plink2_cmd <- paste0(plink2, " --bfile ", bfile, " --extract ", tmp)
  # output template
  locus_input_list <- list()
  # option 1
  if (do_locuszoom_input) {
    # find index variant if it was not provided
    if (is.null(index_var)) {
      index_var <- sumstats[which.max(sumstats[["nlog10P"]]),][["Name"]]
      index_var_nlog10 <- sumstats[which.max(sumstats[["nlog10P"]]),][["nlog10P"]]
      message("Checking if index variant is present in reference panel. If not - using second most significant, third, etc ... ")
      bim <- read.table(paste0(bfile, ".bim"))
      while (! index_var %in% bim[["V2"]]) {
        sumstats <- subset(sumstats, ! Name %in% index_var)
        index_var <- sumstats[which.max(sumstats[["nlog10P"]]),][["Name"]]
        index_var_nlog10 <- sumstats[which.max(sumstats[["nlog10P"]]),][["nlog10P"]]
      }
      message(paste0("Selected index variant is ", index_var, ": -log10P=", index_var_nlog10))
    }
    #
    system(paste0(plink2_cmd, " --out ", tmp, "_r2",
                  " --r2-unphased --ld-snp ", index_var,
                  " --ld-window-kb 9999999 --ld-window-r2 0 "))
    message("Adding r2 to the input sumstats")
    r2.vcor <- read.table(paste0(tmp, "_r2.vcor"), comment.char = "", header=T)
    r2.vcor <- r2.vcor[,c("ID_B", "UNPHASED_R2")]
    colnames(r2.vcor) <- c("Name", "ld")
    locus_input <- merge(sumstats, r2.vcor, all.x=T)
    #
    locus_input_list[["locus_input"]] <- locus_input
    locus_input_list[["index_var"]] <- index_var
  }
  # output
  if (rm_tmp) { message("Cleaning tmp files"); unlink(paste0(tmp, "*"))}
  return(locus_input_list)
}

genepicoloc_rap <- function(locus_input_list, gene, ens_db="ensDb_v106",
                            p="P", chrom="CHR", pos="POS", labs="Name",
                            filter_gene_name = NULL) {
  if (!requireNamespace("AnnotationHub", quietly = TRUE)) {
    stop("The AnnotationHub package must be installed to use this functionality")
  }
  if (!requireNamespace("locuszoomr", quietly = TRUE)) {
    stop("The locuszoomr package must be installed to use this functionality")
  }
  library(AnnotationHub)
  library(locuszoomr)
  ah <- AnnotationHub()
  query(ah, c("EnsDb", "Homo sapiens"))
  ensDb_v106 <- ah[["AH89426"]] # AH100643
  # input variables
  locus_input <- locus_input_list[["locus_input"]]
  index_snp <- locus_input_list[["index_var"]]
  message("locuszoomr processing")
  # reformat
  locus_input <- locus(locus_input,
                       ens_db = ensDb_v106,
                       gene = gene,
                       p=p, chrom=chrom, pos=pos, labs=labs,
                       index_snp=index_snp)
  # correct underflow - tmp
  # locus_input[["data"]][["logP"]] <- asNumeric(locus_input[["data"]][["logP"]])
  p <- locus_ggplot(locus_input, filter_gene_name=filter_gene_name)
  return(p)
}
