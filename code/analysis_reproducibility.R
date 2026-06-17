# Reproduce the main-text tables of the MS200V vs IDEXX Catalyst One method
# comparison study. Pearson r, Deming regression, and Bland-Altman.
#
# The CSV data/paired_results.csv is already filtered to the per-analyte
# cohorts used in the published analysis; no further exclusion is required.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble)
})

root <- dirname(dirname(normalizePath(sys.frame(1)$ofile %||% "code/analysis_reproducibility.R")))
data_path <- file.path(root, "data", "paired_results.csv")
TEA <- c(ALB=15, ALP=25, AMY=25, Ca=10, CHOL=20, CREA=20, GLU=20, PHOS=15, TBIL=25, TP=10, UREA=15)
HYBRID_ABS <- c(TBIL=3.0)

deming <- function(x, y, lam=1.0) {
  n <- length(x); mx <- mean(x); my <- mean(y)
  Sxx <- sum((x-mx)^2); Syy <- sum((y-my)^2); Sxy <- sum((x-mx)*(y-my))
  slope <- (Syy - lam*Sxx + sqrt((Syy - lam*Sxx)^2 + 4*lam*Sxy^2)) / (2*Sxy)
  intercept <- my - slope*mx
  r <- Sxy / sqrt(Sxx*Syy)
  list(slope=slope, intercept=intercept, r=r, n=n)
}

jk_ci <- function(x, y, alpha=0.05) {
  n <- length(x)
  slopes <- sapply(seq_len(n), function(i) deming(x[-i], y[-i])$slope)
  ints   <- sapply(seq_len(n), function(i) deming(x[-i], y[-i])$intercept)
  full <- deming(x, y)
  se_s <- sqrt((n-1)/n * sum((slopes - mean(slopes))^2))
  se_b <- sqrt((n-1)/n * sum((ints - mean(ints))^2))
  tval <- qt(1 - alpha/2, n - 2)
  list(slope=full$slope, slope_lo=full$slope - tval*se_s, slope_hi=full$slope + tval*se_s,
       intercept=full$intercept, intercept_lo=full$intercept - tval*se_b, intercept_hi=full$intercept + tval*se_b,
       r=full$r, n=n)
}

df <- read_csv(data_path, show_col_types=FALSE) %>% filter(!is.na(IDEXX), !is.na(MS200V))
out <- list()
for (a in sort(unique(df$Analyte))) {
  d <- df %>% filter(Analyte == a)
  res <- jk_ci(d$IDEXX, d$MS200V)
  pct <- (d$MS200V - d$IDEXX) / d$IDEXX * 100
  mb <- mean(pct); sdb <- sd(pct)
  tea <- TEA[[a]]
  within <- abs(pct) <= tea
  if (!is.null(HYBRID_ABS[[a]])) {
    within <- within | (abs(d$MS200V - d$IDEXX) <= HYBRID_ABS[[a]])
  }
  out[[a]] <- tibble(
    Analyte=a, N=res$n,
    Deming_slope=round(res$slope, 3),
    Deming_intercept=round(res$intercept, 2),
    Pearson_r=round(res$r, 3),
    R2=round(res$r^2, 3),
    Mean_bias_pct=round(mb, 2),
    SD_bias_pct=round(sdb, 2),
    LoA_lo_pct=round(mb - 1.96*sdb, 1),
    LoA_hi_pct=round(mb + 1.96*sdb, 1),
    Pct_within_TEa=round(mean(within)*100, 1),
    ASVCP_TEa_pct=tea
  )
}
out_df <- bind_rows(out)
write_csv(out_df, file.path(root, "data", "expected_outputs.csv"))
print(out_df, n=Inf)
