# ==============================================================================
# Reproducible R 4.6.2 Verification Script for Three-Track Systematic Review
# Title: Subunit and multi-epitope vaccines against carbapenem-resistant
#        Acinetobacter baumannii: a systematic review and meta-analysis
# Registration: PROSPERO CRD420261572661 | Zenodo: 10.5281/zenodo.24332515
# ==============================================================================

args_all <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", args_all, value = TRUE)
if (length(script_flag) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_flag[1])))
} else {
  script_dir <- getwd()
}
data_dir <- normalizePath(file.path(script_dir, "..", "data"), mustWork = FALSE)

# 1. Track 1: Preclinical Murine Challenge Meta-Analysis (k = 22 experiments, N = 574 mice)
t1 <- read.csv(file.path(data_dir, "Track1_22_experiments_verified.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(t1) == 22, sum(t1$nv) == 390, sum(t1$nc) == 184, sum(t1$ev) == 280, sum(t1$ec) == 13)

# Haldane-Anscombe +0.5 continuity correction on zero-cell 2x2 tables
zero_cell <- (t1$ev == 0) | (t1$ev == t1$nv) | (t1$ec == 0) | (t1$ec == t1$nc)
a <- ifelse(zero_cell, t1$ev + 0.5, t1$ev)
b <- ifelse(zero_cell, (t1$nv - t1$ev) + 0.5, t1$nv - t1$ev)
c <- ifelse(zero_cell, t1$ec + 0.5, t1$ec)
d <- ifelse(zero_cell, (t1$nc - t1$ec) + 0.5, t1$nc - t1$ec)

lnRR <- log((a / (a + b)) / (c / (c + d)))
var_lnRR <- (1 / a) - (1 / (a + b)) + (1 / c) - (1 / (c + d))

# DerSimonian-Laird + HKSJ on lnRR
dl_pool <- function(yi, vi) {
  k <- length(yi)
  wi <- 1 / vi
  mu_fe <- sum(wi * yi) / sum(wi)
  Q <- sum(wi * (yi - mu_fe)^2)
  C_val <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- max(0, (Q - (k - 1)) / C_val)
  wi_re <- 1 / (vi + tau2)
  mu <- sum(wi_re * yi) / sum(wi_re)
  se_dl <- sqrt(1 / sum(wi_re))
  q_hksj <- sum(wi_re * (yi - mu)^2) / ((k - 1) * sum(wi_re))
  se_hksj <- sqrt(max(1 / sum(wi_re), q_hksj))
  q_mult <- sum(wi_re * (yi - mu)^2) / (k - 1)
  se_hksj_std <- se_dl * sqrt(max(1, q_mult))
  tcrit <- qt(0.975, df = k - 1)
  I2 <- max(0, (Q - (k - 1)) / Q) * 100
  list(mu = mu, se_dl = se_dl, se_hksj = se_hksj_std, q = q_mult,
       wald_ci = mu + c(-1, 1) * qnorm(0.975) * se_dl,
       hksj_ci = mu + c(-1, 1) * tcrit * se_hksj_std,
       Q = Q, p_Q = pchisq(Q, df = k - 1, lower.tail = FALSE), I2 = I2, tau2 = tau2)
}

res_rr <- dl_pool(lnRR, var_lnRR)
cat(sprintf("Track 1 Pooled RR_DL = %.4f [Wald 95%% CI: %.4f - %.4f; HKSJ: %.4f - %.4f], Q = %.2f (p = %.4f), I2 = %.1f%%\n",
            exp(res_rr$mu), exp(res_rr$wald_ci[1]), exp(res_rr$wald_ci[2]),
            exp(res_rr$hksj_ci[1]), exp(res_rr$hksj_ci[2]), res_rr$Q, res_rr$p_Q, res_rr$I2))

# Additive Risk Difference (RevMan zero-cell +0.5 table adjustment)
p1_h <- a / (a + b)
p0_h <- c / (c + d)
rd_h <- (p1_h - p0_h) * 100
var_rd_h <- ((p1_h * (1 - p1_h) / (a + b)) + (p0_h * (1 - p0_h) / (c + d))) * 10000
res_rd <- dl_pool(rd_h, var_rd_h)
cat(sprintf("Track 1 Pooled RD_DL = +%.2f pp [Wald: %.2f - %.2f; HKSJ: %.2f - %.2f], Q_RD = %.2f, I2_RD = %.1f%%\n",
            res_rd$mu, res_rd$wald_ci[1], res_rd$wald_ci[2], res_rd$hksj_ci[1], res_rd$hksj_ci[2], res_rd$Q, res_rd$I2))

# 2. Track 2: In Silico MEVC Quality Rubric (k = 22 studies)
t2 <- read.csv(file.path(data_dir, "Track2_22_insilico_rubric_verified.csv"), stringsAsFactors = FALSE)
cat(sprintf("Track 2 Rubric Score (k = %d): Mean = %.2f +/- %.2f / 14 (Median = %.1f, Range = %d-%d)\n",
            nrow(t2), mean(t2$total_verified_14), sd(t2$total_verified_14),
            median(t2$total_verified_14), min(t2$total_verified_14), max(t2$total_verified_14)))

# 3. Track 3: Regional ICU Epidemiology (k = 17 cohorts)
t3 <- read.csv(file.path(data_dir, "Track3_17_cohorts_verified.csv"), stringsAsFactors = FALSE)
eg_unsel <- subset(t3, region_group == "Egypt" & selection_category == "Unselected_Surveillance" & carbapenem_reported == "True")
wt_eg <- sum(as.numeric(eg_unsel$carbapenem_r_n_derived)) / sum(as.numeric(eg_unsel$n_isolates_analyzed)) * 100
cat(sprintf("Track 3 Strictly Unselected Egypt CRAB (k = %d, N = %d): %.2f%%\n",
            nrow(eg_unsel), sum(as.numeric(eg_unsel$n_isolates_analyzed)), wt_eg))
