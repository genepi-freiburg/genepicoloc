test_that("coloc_abf produces valid output structure", {
  set.seed(42)
  n_snps <- 100

  d1 <- list(
    beta = rnorm(n_snps, 0, 0.1),
    varbeta = rep(0.01, n_snps),
    snp = paste0("snp", 1:n_snps),
    type = "quant",
    sdY = 1
  )

  d2 <- list(
    beta = rnorm(n_snps, 0, 0.1),
    varbeta = rep(0.01, n_snps),
    snp = paste0("snp", 1:n_snps),
    type = "quant",
    sdY = 1
  )

  result <- coloc_abf(d1, d2)

  # Check structure
  expect_true(is.list(result))
  expect_true(all(c("summary", "results", "priors") %in% names(result)))

  # Check summary
  expect_equal(result$summary["nsnps"], c(nsnps = 100))
  pp_names <- paste0("PP.H", 0:4, ".abf")
  expect_true(all(pp_names %in% names(result$summary)))

  # PPs should sum to 1
  pp <- result$summary[pp_names]
  expect_equal(sum(pp), 1, tolerance = 1e-10)

  # All PPs should be between 0 and 1
  expect_true(all(pp >= 0 & pp <= 1))

  # Check results data.frame
  expect_equal(nrow(result$results), 100)
  expect_true("snp" %in% names(result$results))
  expect_true("SNP.PP.H4" %in% names(result$results))

  # Per-SNP H4 should sum to PP.H4
  expect_equal(sum(result$results$SNP.PP.H4), as.numeric(pp["PP.H4.abf"]),
               tolerance = 1e-10)
})


test_that("coloc_abf detects colocalization with shared signal", {
  set.seed(42)
  n_snps <- 200

  # Null SNPs
  beta1 <- rnorm(n_snps, 0, 0.01)
  beta2 <- rnorm(n_snps, 0, 0.01)

  # One shared causal SNP (strong signal in both)
  beta1[50] <- 0.5
  beta2[50] <- 0.3

  d1 <- list(
    beta = beta1,
    varbeta = rep(0.005, n_snps),
    snp = paste0("snp", 1:n_snps),
    type = "quant",
    sdY = 1
  )

  d2 <- list(
    beta = beta2,
    varbeta = rep(0.005, n_snps),
    snp = paste0("snp", 1:n_snps),
    type = "quant",
    sdY = 1
  )

  result <- coloc_abf(d1, d2)

  # Should strongly favor H4 (shared causal variant)
  expect_gt(result$summary["PP.H4.abf"], 0.8)

  # Top SNP should be snp50
  top_snp <- result$results$snp[which.max(result$results$SNP.PP.H4)]
  expect_equal(top_snp, "snp50")
})


test_that("coloc_abf detects independent signals (H3)", {
  set.seed(42)
  n_snps <- 200

  beta1 <- rnorm(n_snps, 0, 0.01)
  beta2 <- rnorm(n_snps, 0, 0.01)

  # Different causal SNPs
  beta1[50] <- 0.5
  beta2[150] <- 0.5

  d1 <- list(
    beta = beta1,
    varbeta = rep(0.005, n_snps),
    snp = paste0("snp", 1:n_snps),
    type = "quant",
    sdY = 1
  )

  d2 <- list(
    beta = beta2,
    varbeta = rep(0.005, n_snps),
    snp = paste0("snp", 1:n_snps),
    type = "quant",
    sdY = 1
  )

  result <- coloc_abf(d1, d2)

  # Should favor H3 (independent signals) over H4
  expect_gt(result$summary["PP.H3.abf"], result$summary["PP.H4.abf"])
})


test_that("coloc_abf handles partial SNP overlap", {
  set.seed(42)

  d1 <- list(
    beta = rnorm(100, 0, 0.1),
    varbeta = rep(0.01, 100),
    snp = paste0("snp", 1:100),
    type = "quant",
    sdY = 1
  )

  # Only 50 shared SNPs
  d2 <- list(
    beta = rnorm(80, 0, 0.1),
    varbeta = rep(0.01, 80),
    snp = paste0("snp", 51:130),
    type = "quant",
    sdY = 1
  )

  result <- coloc_abf(d1, d2)
  expect_equal(result$summary["nsnps"], c(nsnps = 50))
})


test_that("coloc_abf handles case-control datasets", {
  set.seed(42)
  n_snps <- 100

  d1 <- list(
    beta = rnorm(n_snps, 0, 0.1),
    varbeta = rep(0.01, n_snps),
    snp = paste0("snp", 1:n_snps),
    type = "cc",
    s = 0.3
  )

  d2 <- list(
    beta = rnorm(n_snps, 0, 0.1),
    varbeta = rep(0.01, n_snps),
    snp = paste0("snp", 1:n_snps),
    type = "quant",
    sdY = 1
  )

  result <- coloc_abf(d1, d2)
  pp <- result$summary[paste0("PP.H", 0:4, ".abf")]
  expect_equal(sum(pp), 1, tolerance = 1e-10)
})


test_that("coloc_abf matches coloc::coloc.abf when available", {
  skip_if_not_installed("coloc")

  set.seed(42)
  n_snps <- 200

  beta1 <- rnorm(n_snps, 0, 0.05)
  beta2 <- rnorm(n_snps, 0, 0.05)
  beta1[50] <- 0.4
  beta2[50] <- 0.3

  d1 <- list(
    beta = beta1,
    varbeta = rep(0.005, n_snps),
    snp = paste0("snp", 1:n_snps),
    type = "quant",
    sdY = 1
  )

  d2 <- list(
    beta = beta2,
    varbeta = rep(0.005, n_snps),
    snp = paste0("snp", 1:n_snps),
    type = "quant",
    sdY = 1
  )

  # Our implementation
  ours <- coloc_abf(d1, d2)

  # Reference implementation
  ref <- suppressWarnings(suppressMessages(
    coloc::coloc.abf(dataset1 = d1, dataset2 = d2)
  ))

  # Compare posteriors
  pp_names <- paste0("PP.H", 0:4, ".abf")
  for (h in pp_names) {
    expect_equal(ours$summary[h], ref$summary[h], tolerance = 1e-6,
                 label = paste("Our", h), expected.label = paste("coloc", h))
  }
})
