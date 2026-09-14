###############################################################################
# Interrogation Duration and the Estimation of False Confession Wrongful
# Conviction Risk: A Reply to Smith and Colleagues
# Mourtgos & Adams -- replication code
###############################################################################

# Reproduces every number and figure in the reply, in the order they appear.
# The duration model is Smith et al.'s (2026, J Crim Justice 107:102742). The
# study inputs and every model function are copied from their deposited script,
#   smith_et_al_2026_deposited_code/2_fcwc_risk_conditional_on_differential_interrogation_length.r
# (github.com/thomasbryansmith/FCWC_Replication, files dated 2026-09-08).
# Nothing in the model is re-specified. We change which studies enter the two
# likelihoods and hold the base rate fixed; that is all.
#
# Run from the package folder (RStudio: Session > Set Working Directory >
# To Source File Location). Writes figures/ and output/results.txt. ~2 min.

# ── Libraries ----------------------------------------------------------------
library(dplyr)   # bind_rows() inside Smith et al.'s bootstrap function
library(ragg)    # 300-dpi PNGs

options(pillar.sigfig = 4)   # so 14.75 prints as 14.75, not 14.8

dir.create("figures", showWarnings = FALSE)
dir.create("output",  showWarnings = FALSE)
sink("output/results.txt", split = TRUE)   # console output is also saved


###############################################################################
# 1. SMITH ET AL. INPUTS ------------------------------------------------------
###############################################################################

# Study-level duration summaries in hours. Their script 2, ll. 284-303, verbatim.
fcwc_studies <- data.frame(study = c("Leo_2004", "Redlich_2011"),
                           mean = c(16.23, 3.07),
                           sd = c(15.05, 4.21),
                           n = c(44, 35))

nonfcwc_studies <- data.frame(study = c("Kassin_2007", "Brimbal_2024",
                                        "Redlich_2011", "Cleary&Bull_2021",
                                        "Kelly_et_al_2016", "Cleary_2014",
                                        "Russano_et_al_2024_ncf",
                                        "Russano_et_al_2024_pcf",
                                        "Russano_et_al_2024_fcf",
                                        "Russano_et_al_2026_ncf",
                                        "Russano_et_al_2026_pcf",
                                        "Russano_et_al_2026_fcf"),
                              mean = c(1.6, 1.5, 1.94, 1.49, 1.53028, 1.093333, 1.1345, 1.2991667,
                                       1.3268333, 1.658333, 1.73233333, 1.81616667),
                              sd = c(0.89, 1.17, 3.14, 1.84, 1.046667, 0.983333, 0.77066667, 0.9425,
                                     1.0845, 1.72066667, 1.73416667, 1.4491667),
                              n = c(601, 484, 30, 249, 29, 57, 40, 63, 38, 34, 78, 36))

# Their l. 288 drops Drizin & Leo before anything is fitted. The published
# FCWC likelihood is Redlich et al. (2011) alone.
fcwc_published <- fcwc_studies[which(fcwc_studies$study == "Redlich_2011"), ]

# Base rate: their posterior median at full attribution (p. 4). Held fixed
# throughout so that only the duration likelihoods move.
pi_fcwc <- 0.019

# What each non-FCWC row counts. Checked against each article.
nonfcwc_studies$source <- c("officer survey", "officer survey",
                            "suspect self-report", "suspect self-report",
                            "recorded", "recorded",
                            rep("recorded", 6))

# Drizin & Leo (2004) Table 7, p. 949: reported interrogation length, N = 44
dl_table7 <- data.frame(hours = c("<6", "6-12", "12-24", "24-48", "48-72", "72-96"),
                        n = c(7, 15, 17, 3, 1, 1))


###############################################################################
# 2. MODEL FUNCTIONS ----------------------------------------------------------
###############################################################################

# ── Smith et al., script 2, ll. 306-519 -- verbatim, only re-indented --------

### Convert to lognormal
lognormal <- function(mean, sd){
  sdlog <- sqrt(log(1 + (sd / mean)^2))
  meanlog <- log(mean) - 0.5 * sdlog^2
  list(meanlog = meanlog,
       sdlog = sdlog)
}

### Mixture distributions
mixture_dist <- function(studies, weighting = c("sample_size",
                                                "equal")){
  weighting <- match.arg(weighting)
  pars <- lognormal(mean = studies$mean,
                    sd = studies$sd)
  studies$meanlog <- pars$meanlog
  studies$sdlog   <- pars$sdlog
  if(weighting == "sample_size"){
    studies$weight <- studies$n /
      sum(studies$n)
  }else{
    studies$weight <- rep(
      1 / nrow(studies),
      nrow(studies)
    )
  }
  return(studies)
}

#### Truncated lognormal distribution
rlnorm_trunc <- function(n, meanlog, sdlog, lower = 0, upper = 48) {
  if(lower < 0) stop("Lower bound must be >= 0 for log-normal")

  # CDF bounds
  p_lower <- if(lower == 0) 0 else plnorm(lower,
                                          meanlog = meanlog,
                                          sdlog = sdlog)
  p_upper <- plnorm(upper,
                    meanlog = meanlog,
                    sdlog = sdlog)
  if (p_upper <= p_lower)
    stop("Invalid bounds: upper CDF <= lower CDF")

  # Draw from truncated uniform
  u <- runif(n, min = p_lower, max = p_upper)

  # Inverse CDF
  return(qlnorm(u, meanlog, sdlog))
}

#### Truncated lognormal kernel density
dlnorm_trunc <- function(x, meanlog, sdlog, lower = 0, upper = 48){
  Z <- plnorm(upper, meanlog = meanlog,
              sdlog = sdlog) -
    plnorm(lower, meanlog = meanlog,
           sdlog = sdlog)

  dens <- dlnorm(x, meanlog = meanlog,
                 sdlog = sdlog) / Z

  return(ifelse(x <= lower | x >= upper, 0, dens))
}

#### Component likelihoods
study_likelihood_H <- function(h, study_mix,
                               lower = 0,upper = 48){
  likelihood <- mapply(FUN = dlnorm_trunc,
                       meanlog = study_mix$meanlog,
                       sdlog = study_mix$sdlog,
                       MoreArgs = list(
                         x = h,
                         lower = lower,
                         upper = upper))

  return(data.frame(study = study_mix$study,
                    likelihood = likelihood,
                    weight = study_mix$weight,
                    weighted_likelihood =
                      likelihood * study_mix$weight))
}

#### Marginal study-mixture likelihood
mixture_likelihood_H <- function(h, study_mix,
                                 lower = 0, upper = 48){

  component_likelihoods <- study_likelihood_H(h = h,
                                              study_mix = study_mix,
                                              lower = lower,
                                              upper = upper)

  return(sum(component_likelihoods$weighted_likelihood))
}

# Bayes theorem with study-mixture likelihoods
posterior_fcwc_given_H_mixture <- function(h, pi,
                                           fcwc_mix, nonfcwc_mix,
                                           lower = 0, upper = 48){
  L1 <- mixture_likelihood_H(h = h, study_mix = fcwc_mix,
                             lower = lower, upper = upper)

  L0 <- mixture_likelihood_H(h = h, study_mix = nonfcwc_mix,
                             lower = lower, upper = upper)

  posterior <- (L1 * pi) / (L1 * pi + L0 * (1 - pi))

  data.frame(L1 = L1, L0 = L0,
             LR = L1 / L0, pi = pi,
             posterior = posterior)
}

# maximum likelihood for banding
fit_tlnorm_mle <- function(x, lower = 0, upper = 48){
  x <- x[is.finite(x) & x > lower & x < upper]

  nll <- function(par){
    meanlog <- par[1]
    sdlog   <- exp(par[2])

    p_lower <- if (lower == 0) 0 else{
      plnorm(lower,
             meanlog = meanlog,
             sdlog = sdlog)
    }

    p_upper <- plnorm(upper,
                      meanlog = meanlog,
                      sdlog = sdlog)

    Z <- p_upper - p_lower

    if (!is.finite(Z) || Z <= 0) return(Inf)

    -sum(dlnorm(x, meanlog = meanlog,
                sdlog = sdlog,
                log = TRUE) -
           log(Z))
  }

  log_x <- log(x)

  initial <- c(mean(log_x), log(sd(log_x)))

  fit <- optim(par = initial,
               fn = nll,
               method = "BFGS")
  c(meanlog = fit$par[1],
    sdlog = exp(fit$par[2]))
}

# Bootstrapped Mixture Dist
bootstrap_mixture_likelihood <- function(hours, study_mix,
                                         n_rep = 1000,
                                         lower = 0, upper = 48){
  results <- vector("list", n_rep)

  for (r in seq_len(n_rep)){

    boot_pars <- lapply(seq_len(nrow(study_mix)),
                        function(j){
                          x_boot <- rlnorm_trunc(n = study_mix$n[j],
                                                 meanlog = study_mix$meanlog[j],
                                                 sdlog = study_mix$sdlog[j],
                                                 lower = lower,
                                                 upper = upper)

                          fit_tlnorm_mle(x = x_boot, lower = lower, upper = upper)})

    boot_pars <- do.call(rbind, boot_pars)

    component_density <- vapply(seq_len(nrow(study_mix)),
                                function(j) {
                                  dlnorm_trunc(x = hours,
                                               meanlog = boot_pars[j, "meanlog"],
                                               sdlog = boot_pars[j, "sdlog"],
                                               lower = lower,
                                               upper = upper)},
                                FUN.VALUE = numeric(length(hours)))

    mixture_density <- as.vector(component_density %*% study_mix$weight)

    results[[r]] <- data.frame(iteration = r, hour = hours, likelihood = mixture_density)}

  dplyr::bind_rows(results)
}

# ── Ours ---------------------------------------------------------------------

# posterior at h, base rate fixed
post_h <- function(h, fcwc_mix, nonfcwc_mix) {
  posterior_fcwc_given_H_mixture(h = h, pi = pi_fcwc,
                                 fcwc_mix = fcwc_mix,
                                 nonfcwc_mix = nonfcwc_mix)$posterior
}

# share of a fitted (0-48 h truncated) mixture at or above h
mass_above <- function(h, mix) {
  sum(mix$weight * (1 - plnorm(h, mix$meanlog, mix$sdlog) /
                          plnorm(48, mix$meanlog, mix$sdlog)))
}

# posterior over the whole 0-48 h range on a 0.01 h grid
h_grid <- seq(0.05, 47.95, by = 0.01)

posterior_curve <- function(fcwc_mix, nonfcwc_mix, h = h_grid) {
  sapply(h, post_h, fcwc_mix = fcwc_mix, nonfcwc_mix = nonfcwc_mix)
}

# first duration at which the posterior rises through a threshold. Scanned,
# not uniroot: the curve is not monotone for some specifications (it can sit
# above 10% at very short durations as well)
crossing <- function(p, threshold = .10, h = h_grid) {
  up <- which(diff(sign(p - threshold)) > 0) + 1
  if (length(up) == 0) NA_real_ else h[up[1]]
}

fcwc_mix    <- mixture_dist(fcwc_published, weighting = "sample_size")
nonfcwc_mix <- mixture_dist(nonfcwc_studies, weighting = "sample_size")


###############################################################################
# 3. FIDELITY CHECK -----------------------------------------------------------
###############################################################################

# Their published posteriors (pp. 9-10): 0.012 at 1.53 h, 0.027 at 3.07 h,
# 0.10 at 6.28 h, 0.25 at 12.06 h, 0.50 at 33.96 h. They average over the
# base-rate posterior; pi is fixed at its median here, hence the small offsets.
fidelity <- tibble(hours = c(1.53, 3.07, 6.28, 12.06, 33.96),
                   published = c("0.012", "0.027", "0.10", "0.25", "0.50"),
                   fixed_pi = round(sapply(hours, post_h, fcwc_mix, nonfcwc_mix), 4))
print(fidelity)


###############################################################################
# 4. THE FALSE CONFESSION SAMPLE ----------------------------------------------
###############################################################################

# Redlich et al. (2011), total questioning: false confessors are the published
# FCWC arm; the true confessors are one row of the non-FCWC arm
print(rbind(fcwc_published[, 1:4],
            nonfcwc_studies[which(nonfcwc_studies$study == "Redlich_2011"), 1:4]))

# Fitted Redlich distribution: 12.1% of its mass at 6 h or more, 3.2% at 12 h
# or more. Drizin & Leo Table 7: 84% in the 6 h+ categories, 50% in the 12 h+
# categories. Their 6-12 h band includes 12, so 12 h+ is "at least 50%".
fc_tails <- tibble(hours = c(6, 12),
                   fitted_redlich_pct = round(100 * sapply(hours, mass_above, mix = fcwc_mix), 1),
                   drizin_leo_pct = round(100 * c(sum(dl_table7$n[-1]),
                                                  sum(dl_table7$n[-(1:2)])) / sum(dl_table7$n), 1))
print(fc_tails)


###############################################################################
# 5. THE COMPARISON SAMPLE ----------------------------------------------------
###############################################################################

# Sample-size weight by source type: 62.4% officer surveys, 16.0% suspect
# self-reports, 21.6% recorded interrogations
weights_by_source <- nonfcwc_mix %>%
  group_by(source) %>%
  summarise(rows = n(),
            n = sum(n),
            weight_pct = round(100 * sum(weight), 1),
            .groups = "drop") %>%
  arrange(desc(weight_pct))
print(weights_by_source)

# Kassin (601 officers) against Kelly (29 recorded interrogations): ~21x
print(nonfcwc_mix$weight[nonfcwc_mix$study == "Kassin_2007"] /
        nonfcwc_mix$weight[nonfcwc_mix$study == "Kelly_et_al_2016"])


###############################################################################
# 6. THE REPORTED THRESHOLDS DEPEND ON THE INPUTS ------------------------------
###############################################################################

# Published specification with pi fixed: 10.6% at 6.28 h, crosses 10% at 6.07 h
post_published <- posterior_curve(fcwc_mix, nonfcwc_mix)
print(c(at_6.28_h = post_h(6.28, fcwc_mix, nonfcwc_mix),
        crossing_10pct_h = crossing(post_published)))

# Swap only the FCWC source. "Both" uses their own n-weighting:
# Drizin & Leo 44/79 = 0.557, Redlich 35/79 = 0.443
fcwc_arms <- list(
  "Redlich only (published)" = fcwc_mix,
  "Both, n-weighted"         = mixture_dist(fcwc_studies, weighting = "sample_size"),
  "Drizin & Leo only"        = mixture_dist(fcwc_studies[which(fcwc_studies$study == "Leo_2004"), ],
                                            weighting = "sample_size"))
print(fcwc_arms[["Both, n-weighted"]][, c("study", "n", "weight")])

# At 12.06 h: 26.3% / 62.1% / 72.6% (Figure 2)
source_choice <- bind_rows(lapply(names(fcwc_arms), function(a) {
  p <- posterior_curve(fcwc_arms[[a]], nonfcwc_mix)
  tibble(fcwc_source = a,
         at_6.28_h = round(100 * post_h(6.28, fcwc_arms[[a]], nonfcwc_mix), 1),
         at_12.06_h = round(100 * post_h(12.06, fcwc_arms[[a]], nonfcwc_mix), 1),
         crossing_10pct_h = crossing(p))
}))
print(source_choice)


###############################################################################
# 7. COMPARING DURATIONS WITHIN THE SAME STUDY --------------------------------
###############################################################################

# Redlich's 30 true confessors (M = 1.94, SD = 3.14) as the comparison arm.
# 3.6% at 6.28 h, 4.0% at 12.06 h; never above ~4.2%, never crosses 10%
redlich_tc_mix <- mixture_dist(nonfcwc_studies[which(nonfcwc_studies$study == "Redlich_2011"), ],
                               weighting = "sample_size")
post_within <- posterior_curve(fcwc_mix, redlich_tc_mix)

within_study <- tibble(comparison_arm = c("Smith et al. pooled", "Redlich true confessors"),
                       at_6.28_h = round(100 * c(post_h(6.28, fcwc_mix, nonfcwc_mix),
                                                 post_h(6.28, fcwc_mix, redlich_tc_mix)), 1),
                       at_12.06_h = round(100 * c(post_h(12.06, fcwc_mix, nonfcwc_mix),
                                                  post_h(12.06, fcwc_mix, redlich_tc_mix)), 1))
print(within_study)
print(c(within_max_pct_0_48_h = round(100 * max(post_within), 1),
        within_crossing_10pct_h = crossing(post_within)))


###############################################################################
# 8. WHERE THE DATA ARE NOT OBSERVED ------------------------------------------
###############################################################################

# Fitted non-FCWC mass above 12.06 h (0.116%) and above 34 h (0.0024%)
nonfcwc_tails <- tibble(hours = c(12.06, 34),
                        nonfcwc_above_pct = signif(100 * sapply(hours, mass_above, mix = nonfcwc_mix), 3))
print(nonfcwc_tails)

# Dispersion: every non-FCWC sdlog scaled by k, means unchanged.
# 10% crossing 3.93 h at 0.75, 6.07 h at 1, 14.75 h at 1.25, none at 1.5
scale_sdlog <- function(mix, k) {
  mix$sdlog <- mix$sdlog * k
  mix
}

dispersion <- bind_rows(lapply(c(0.75, 1, 1.25, 1.5), function(k) {
  p <- posterior_curve(fcwc_mix, scale_sdlog(nonfcwc_mix, k))
  tibble(sdlog_x = k,
         crossing_10pct_h = crossing(p),
         max_pct_0_48_h = round(100 * max(p), 1))
}))
print(dispersion)

# Composition of the non-FCWC arm. 10% crossing: 7.20 h with equal weights,
# 8.69 h without the officer surveys, 6.62 h with recorded interrogations only
nonfcwc_arms <- list(
  "equal weights"      = mixture_dist(nonfcwc_studies, weighting = "equal"),
  "no officer surveys" = mixture_dist(nonfcwc_studies[nonfcwc_studies$source != "officer survey", ],
                                      weighting = "sample_size"),
  "recorded only"      = mixture_dist(nonfcwc_studies[nonfcwc_studies$source == "recorded", ],
                                      weighting = "sample_size"))

composition <- bind_rows(lapply(names(nonfcwc_arms), function(a) {
  p <- posterior_curve(fcwc_mix, nonfcwc_arms[[a]])
  tibble(nonfcwc_arm = a,
         rows = nrow(nonfcwc_arms[[a]]),
         crossing_10pct_h = crossing(p))
}))
print(composition)

# Likelihood uncertainty. Smith et al.'s parametric bootstrap (the shading in
# their Fig. 2a) run on both arms and pushed through Bayes' rule with pi fixed.
# 400 replications. Seed and call order are fixed: the intervals in the text
# came from this run, and 3.07 h (Redlich's mean) goes first -- dropping it
# changes the draws for the other two.
# 6.28 h: 5.7-14.6%; 12.06 h: 7.1-42.3%
set.seed(20260912)
boot <- bind_rows(lapply(c(3.07, 6.28, 12.06), function(h) {
  L1 <- bootstrap_mixture_likelihood(hours = h, study_mix = fcwc_mix,    n_rep = 400)$likelihood
  L0 <- bootstrap_mixture_likelihood(hours = h, study_mix = nonfcwc_mix, n_rep = 400)$likelihood
  p  <- (L1 * pi_fcwc) / (L1 * pi_fcwc + L0 * (1 - pi_fcwc))
  tibble(hours = h,
         point = post_h(h, fcwc_mix, nonfcwc_mix),
         boot_median = median(p),
         lower_95 = unname(quantile(p, .025)),
         upper_95 = unname(quantile(p, .975)))
}))
print(boot %>% mutate(across(-hours, ~ round(100 * .x, 1))))   # in percent


###############################################################################
# 9. FIGURES ------------------------------------------------------------------
###############################################################################

# Okabe-Ito, colour-blind safe; every series also gets its own line type
col_pool <- "#000000"; col_red <- "#009E73"; col_hi <- "#0072B2"; col_lo <- "#D55E00"; col_15 <- "#999999"

# --- Figure 1: the inverse-probability problem ---
agg_png("figures/Figure1_inverse_probability_problem.png", width = 7.2, height = 4.9, units = "in", res = 300)
par(mar = c(0.2, 0.2, 0.2, 0.2), family = "sans")
plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

x0 <- 0.33; x1 <- 0.56; x2 <- 0.79          # column edges
yT <- 0.90; yM <- 0.68; yB <- 0.46          # row edges
blue <- "#0072B2"; orange <- "#D55E00"; grey_fill <- "#E6E6E6"

# cells; B is the unobserved one
rect(x0, yM, x1, yT, col = "white", border = "grey40", lwd = 1.2)           # A
rect(x1, yM, x2, yT, col = grey_fill, border = "grey40", lwd = 1.2)         # B
rect(x1, yM, x2, yT, col = "grey55", border = NA, density = 14, angle = 45)  # hatch B
rect(x0, yB, x1, yM, col = "white", border = "grey40", lwd = 1.2)           # C
rect(x1, yB, x2, yM, col = "white", border = "grey40", lwd = 1.2)           # D

text((x0+x1)/2, (yM+yT)/2, "A", cex = 2.0, font = 2)
text((x1+x2)/2, (yM+yT)/2 + 0.035, "B", cex = 2.0, font = 2, col = "grey20")
text((x1+x2)/2, (yM+yT)/2 - 0.055, "not observed", cex = 0.85, col = "grey20", font = 3)
text((x0+x1)/2, (yB+yM)/2, "C", cex = 2.0, font = 2)
text((x1+x2)/2, (yB+yM)/2, "D", cex = 2.0, font = 2)

# headers
text((x0+x1)/2, yT + 0.045, "FCWC", cex = 1.05, font = 2)
text((x1+x2)/2, yT + 0.045, "No FCWC", cex = 1.05, font = 2)
text(x0 - 0.045, (yM+yT)/2, "Prolonged\ninterrogation", cex = 1.0, font = 2, adj = c(1, 0.5))
text(x0 - 0.045, (yB+yM)/2, "Other\ninterrogation", cex = 1.0, font = 2, adj = c(1, 0.5))

# column A+C = what outcome-selected studies estimate; row A+B = what the claim needs
pad <- 0.012
rect(x0 - pad, yB - pad, x1 + pad, yT + pad, border = blue,   lwd = 3.2, lty = "22")
rect(x0 - 2*pad, yM - pad, x2 + 2*pad, yT + 2*pad, border = orange, lwd = 3.2, lty = "solid")

# annotations
segments(0.05, 0.335, 0.11, 0.335, col = blue, lwd = 3.2, lty = "22")
text(0.13, 0.335, expression(paste("What outcome-selected studies estimate:   ",
     italic(P), "(prolonged | FCWC) = ", italic(A), " / (", italic(A), " + ", italic(C), ")")),
     adj = c(0, 0.5), cex = 0.95)
segments(0.05, 0.245, 0.11, 0.245, col = orange, lwd = 3.2)
text(0.13, 0.245, expression(paste("What the risk claim requires:   ",
     italic(P), "(FCWC | prolonged) = ", italic(A), " / (", italic(A), " + ", italic(B), ")")),
     adj = c(0, 0.5), cex = 0.95)
rect(0.05, 0.13, 0.11, 0.17, col = grey_fill, border = "grey40"); rect(0.05, 0.13, 0.11, 0.17, col = "grey55", border = NA, density = 14, angle = 45)
text(0.13, 0.175, expression(paste(italic(B), ":  prolonged interrogations that did not produce an FCWC,")),
     adj = c(0, 0.5), cex = 0.95)
text(0.13, 0.125, "the quantity outcome-selected data cannot supply", adj = c(0, 0.5), cex = 0.95)
invisible(dev.off())


# --- Figure 2: which false-confession source builds the FCWC likelihood ---
# pi, the pooled non-FCWC arm, the lognormal family and the 0-48 h truncation
# are all held; only the FCWC arm changes
src <- list(
  list(lab = "Published specification: Redlich et al. (2011) only", A = fcwc_arms[["Redlich only (published)"]], col = col_pool, lty = "solid",   lwd = 2.6),
  list(lab = "Both sources retained, sample-size weighted",         A = fcwc_arms[["Both, n-weighted"]],         col = col_hi,   lty = "dotdash", lwd = 2.4),
  list(lab = "Drizin & Leo (2004) only",                            A = fcwc_arms[["Drizin & Leo only"]],        col = col_lo,   lty = "dashed",  lwd = 2.4))
h4 <- seq(0.25, 20, by = 0.02)   # curves start at 15 minutes; crossings still come from the full grid
for (i in seq_along(src)) {
  src[[i]]$y    <- posterior_curve(src[[i]]$A, nonfcwc_mix, h = h4)
  src[[i]]$p628 <- post_h(6.28,  src[[i]]$A, nonfcwc_mix)
  src[[i]]$p12  <- post_h(12.06, src[[i]]$A, nonfcwc_mix)
  src[[i]]$xc   <- crossing(posterior_curve(src[[i]]$A, nonfcwc_mix))
}

agg_png("figures/Figure2_fcwc_source_choice.png", width = 10, height = 5.8, units = "in", res = 300)
par(mar = c(4.6, 5.0, 1.2, 1.0), mgp = c(2.9, 0.75, 0), las = 1, family = "sans")
plot(NA, xlim = c(0, 20), ylim = c(0, 1), xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
     xlab = "Interrogation duration (hours)", ylab = "Estimated probability of FCWC", cex.lab = 1.15)
axis(1, at = seq(0, 20, 2), cex.axis = 1.0); axis(2, at = seq(0, 1, 0.2), labels = paste0(seq(0, 100, 20), "%"), cex.axis = 1.0)
abline(h = seq(0.2, 0.8, 0.2), col = "grey92"); abline(v = seq(2, 18, 2), col = "grey92")
abline(h = 0.10, col = "grey40", lwd = 1.4, lty = "dashed")
segments(c(6.28, 12.06), 0, c(6.28, 12.06), 0.78, col = "grey55", lwd = 1.2, lty = "dotted")   # stop below the legend table
text(6.28,  0.022, "6.28 h",  adj = c(0.5, 0), cex = 0.95, col = "grey30")
text(12.06, 0.022, "12.06 h", adj = c(0.5, 0), cex = 0.95, col = "grey30")
text(19.8, 0.10, "10% threshold", adj = c(1, -0.5), cex = 0.95, col = "grey30")
for (s_ in src) lines(h4, s_$y, col = s_$col, lty = s_$lty, lwd = s_$lwd)

# value labels; the offsets keep the three 6.28 h labels apart and off the dashed line
off628 <- c(0.022, -0.004, 0.022)
for (k in seq_along(src)) { s_ <- src[[k]]
  points(c(6.28, 12.06), c(s_$p628, s_$p12), pch = 21, bg = "white", col = s_$col, cex = 1.4, lwd = 2.2)
  text(6.28  + 0.25, s_$p628 + off628[k], sprintf("%.1f%%", 100*s_$p628), adj = c(0, 0.5), cex = 1.0, font = 2, col = s_$col)
  text(12.06 + 0.25, s_$p12,              sprintf("%.1f%%", 100*s_$p12),  adj = c(0, 0.5), cex = 1.0, font = 2, col = s_$col)
}

# legend as a small table: 6.28 h, 12.06 h and the crossing side by side
short <- c("Redlich et al. (2011) only — published", "Both sources, n-weighted", "Drizin & Leo (2004) only")
hdr   <- c("at 6.28 h", "at 12.06 h", "reaches 10%")
cx <- 0.95
# column right edges from rendered text widths so nothing overlaps
w1 <- max(strwidth(hdr[1], cex = cx, font = 2), strwidth("19.5%", cex = cx, font = 2))
w2 <- max(strwidth(hdr[2], cex = cx, font = 2), strwidth("72.6%", cex = cx, font = 2))
w3 <- max(strwidth(hdr[3], cex = cx, font = 2), strwidth("6.07 h", cex = cx))
gap <- 0.55
xlab <- 1.45; xl0 <- 0.45; xl1 <- 1.25
xc1 <- xlab + max(strwidth(short, cex = cx)) + gap + w1
xc2 <- xc1 + gap + w2
xc3 <- xc2 + gap + w3
yh <- 0.975; dy <- 0.048
text(xlab, yh, "FCWC likelihood source", adj = c(0, 0.5), cex = cx, font = 2)
text(xc1, yh, hdr[1], adj = c(1, 0.5), cex = cx, font = 2)
text(xc2, yh, hdr[2], adj = c(1, 0.5), cex = cx, font = 2)
text(xc3, yh, hdr[3], adj = c(1, 0.5), cex = cx, font = 2)
segments(xl0, yh - dy/2 - 0.004, xc3 + 0.15, yh - dy/2 - 0.004, col = "grey60")
for (k in seq_along(src)) { s_ <- src[[k]]; yk <- yh - k*dy
  segments(xl0, yk, xl1, yk, col = s_$col, lty = s_$lty, lwd = s_$lwd)
  text(xlab, yk, short[k], adj = c(0, 0.5), cex = cx)
  text(xc1, yk, sprintf("%.1f%%", 100*s_$p628), adj = c(1, 0.5), cex = cx, font = 2, col = s_$col)
  text(xc2, yk, sprintf("%.1f%%", 100*s_$p12),  adj = c(1, 0.5), cex = cx, font = 2, col = s_$col)
  text(xc3, yk, sprintf("%.2f h", s_$xc),       adj = c(1, 0.5), cex = cx)
}
invisible(dev.off())


# --- Figure 3: sensitivity (A) and likelihood uncertainty (B) ---
# Panel A curves from the model, drawn from 15 minutes; crossings and maxima
# from the full 0-48 h grid
hgrid <- seq(0.25, 16, by = 0.02)
specs <- list(
  list(lab = "Smith et al. pooled specification",        m0 = nonfcwc_mix,                    col = col_pool, lty = "solid",    lwd = 2.6),
  list(lab = "Within-study comparison (Redlich 2011)",   m0 = redlich_tc_mix,                 col = col_red,  lty = "dashed",   lwd = 2.2),
  list(lab = "Non-FCWC dispersion + 25%",                m0 = scale_sdlog(nonfcwc_mix, 1.25), col = col_hi,   lty = "dotdash",  lwd = 2.2),
  list(lab = "Non-FCWC dispersion + 50%",                m0 = scale_sdlog(nonfcwc_mix, 1.50), col = col_15,   lty = "longdash", lwd = 2.0),
  list(lab = "Non-FCWC dispersion - 25%",                m0 = scale_sdlog(nonfcwc_mix, 0.75), col = col_lo,   lty = "dotted",   lwd = 2.6))
for (i in seq_along(specs)) {
  p_full        <- posterior_curve(fcwc_mix, specs[[i]]$m0)
  specs[[i]]$y  <- posterior_curve(fcwc_mix, specs[[i]]$m0, h = hgrid)
  specs[[i]]$xc <- crossing(p_full)
  specs[[i]]$mx <- max(p_full)
}

# Panel B: the bootstrap table from section 8, unrounded
pts <- boot$point[boot$hours %in% c(6.28, 12.06)]
los <- boot$lower_95[boot$hours %in% c(6.28, 12.06)]
his <- boot$upper_95[boot$hours %in% c(6.28, 12.06)]

agg_png("figures/Figure3_sensitivity_and_uncertainty.png", width = 11, height = 4.9, units = "in", res = 300)
layout(matrix(1:2, 1, 2), widths = c(1.55, 1))
par(family = "sans")

# Panel A
par(mar = c(4.4, 4.6, 2.2, 1.0), mgp = c(2.6, 0.7, 0), las = 1)
plot(NA, xlim = c(0, 16), ylim = c(0, 1), xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
     xlab = "Interrogation duration (hours)", ylab = "Estimated probability of FCWC", cex.lab = 1.05)
axis(1, at = seq(0, 16, 2), cex.axis = 0.95); axis(2, at = seq(0, 1, 0.2), labels = paste0(seq(0, 100, 20), "%"), cex.axis = 0.95)
abline(h = seq(0.2, 0.8, 0.2), col = "grey92"); abline(v = seq(2, 14, 2), col = "grey92")
abline(h = 0.10, col = "grey40", lwd = 1.3, lty = "solid")
text(8.2, 0.10, "10% threshold", adj = c(0, -0.45), cex = 0.85, col = "grey30")
for (s in specs) lines(hgrid, s$y, col = s$col, lty = s$lty, lwd = s$lwd)
# crossing markers
for (s in specs) if (!is.na(s$xc) && s$xc <= 16) {
  points(s$xc, 0.10, pch = 21, bg = "white", col = s$col, cex = 1.35, lwd = 2)
  text(s$xc, 0.10, sprintf("%.2f h", s$xc), pos = 3, offset = 0.55, cex = 0.85, col = s$col, font = 2)
}
# the two that never cross are labelled in the legend, not on the crowded right edge
leg_lab <- sapply(specs, `[[`, "lab")
leg_lab[2] <- sprintf("%s: no crossing, max %.1f%%", leg_lab[2], 100*specs[[2]]$mx)
leg_lab[4] <- sprintf("%s: no crossing, max %.1f%%", leg_lab[4], 100*specs[[4]]$mx)
legend("topleft", inset = c(0.012, 0.015), bty = "n", cex = 0.86, seg.len = 3.2,
       legend = leg_lab, col = sapply(specs, `[[`, "col"),
       lty = sapply(specs, `[[`, "lty"), lwd = sapply(specs, `[[`, "lwd"))
mtext("A", side = 3, line = 0.5, adj = 0, font = 2, cex = 1.2)

# Panel B
par(mar = c(4.4, 4.6, 2.2, 1.0), mgp = c(2.6, 0.7, 0), las = 1)
xs <- c(1, 2)
plot(NA, xlim = c(0.45, 2.8), ylim = c(0, 0.5), xaxt = "n", yaxt = "n", yaxs = "i",
     xlab = "Interrogation duration", ylab = "Estimated probability of FCWC", cex.lab = 1.05)
axis(1, at = xs, labels = c("6.28 hours", "12.06 hours"), cex.axis = 0.95)
axis(2, at = seq(0, 0.5, 0.1), labels = paste0(seq(0, 50, 10), "%"), cex.axis = 0.95)
abline(h = seq(0.1, 0.4, 0.1), col = "grey92")
abline(h = 0.10, col = "grey40", lwd = 1.3)
text(2.77, 0.10, "10% threshold", adj = c(1, -0.45), cex = 0.85, col = "grey30")
arrows(xs, los, xs, his, angle = 90, code = 3, length = 0.07, lwd = 2.2, col = col_pool)
points(xs, pts, pch = 21, bg = "white", col = col_pool, cex = 1.7, lwd = 2.2)
text(xs + 0.13, pts, sprintf("%.1f%%", 100*pts), adj = c(0, 0.5), cex = 0.9, font = 2)
text(xs + 0.13, his, sprintf("%.1f%%", 100*his), adj = c(0, 0.5), cex = 0.8, col = "grey25")
text(xs + 0.13, los, sprintf("%.1f%%", 100*los), adj = c(0, 0.5), cex = 0.8, col = "grey25")
text(1.5, 0.468, "Point: fixed base rate (0.019)\nBar: 95% interval from the authors' bootstrap\nover study-level duration parameters", cex = 0.78, col = "grey25")
mtext("B", side = 3, line = 0.5, adj = 0, font = 2, cex = 1.2)
invisible(dev.off())

# what went into the figures
print(tibble(figure_3_panel_A = sapply(specs, `[[`, "lab"),
             crossing_10pct_h = sapply(specs, `[[`, "xc"),
             max_pct_0_48_h = round(100 * sapply(specs, `[[`, "mx"), 1)))
list.files("figures")

sink()
