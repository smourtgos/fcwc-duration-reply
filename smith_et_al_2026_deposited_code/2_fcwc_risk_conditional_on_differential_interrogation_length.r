#==============================================================================#

# load packages (simulations) ====
library(dplyr)
library(brms)
library(tidyr)
library(ggplot2)
library(ggridges)
library(posterior)
library(bayestestR)
library(bayesplot)
library(purrr)
library(truncdist)
library(ggh4x)

#==============================================================================#
# Read Data                                                                ====#
#==============================================================================#

# Catlin (2026)
## Structure
### MAref: ID in Mourtgos & Adams' data
### Ref: ID
### StudyYear: publication year
### DataYear: range of years that data were collected
### Included: included in Mourtgos & Adams' data
### Estimate: FC ("FCshare") v. WC ("error")
### Location: data geolocation
### Participants: respondent type
### N: sample size
### bin: bin i in study j
### nbin: respondents in bin i
### %error: % WC
### nerror: N WC
### weight: evidence-type for weighting
### shareddata: data intersections between studies
### notes: notes about data
### comments: additional observations

mar_dat <- readxl::read_xlsx("./data/Data_2026_02_25_cleaned.xlsx")
mar_dat <- mar_dat[which(mar_dat$include),]

wc_dat <- mar_dat[which(mar_dat$Estimate == "Error"),]

fc_dat <- mar_dat[which(mar_dat$Estimate == "FCshare"),]
fc_dat <- fc_dat %>%
  group_by(Ref, Location) %>%
  mutate(p_est_ungr = p_est,
         n_ungr = n,
         N_ungr = N,
         n = sum(n),
         N = sum(N),
         p_est = weighted.mean(p_est, N)) %>%
  ungroup() %>%
  filter(US == 0,
         Ref != "Henkel et al. (2008a)")

#==============================================================================#
# 3.1.1 Estimating False Confession - Wrongful Conviction (FCWC)           ====#
#==============================================================================#

# 3.1.1.1. Conviction Error Rate ====

  ## Define function(s)
  min_max_scale <- function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
  
  ## Generate Period variable
  wc_dat$Period <- as.numeric(gsub("[0-9]{4}-","", wc_dat$DataYear))
  wc_dat$Period[is.na(wc_dat$Period)] <- wc_dat$StudyYear[is.na(wc_dat$Period)]
  
  ## Rescale p_est, create weight variable(s)
  epsilon <- 1e-6
  wc_dat <- wc_dat %>%
    mutate(
      recency = min_max_scale(Period),
      pij = pmin(pmax(p_est, epsilon), 1 - epsilon),
      type_wt = ifelse(weight == "secondary", 1, 0.9),
      wij = recency * type_wt * (n / N),  # total weight
    )
  
  ## formula
  wc_formula <- brmsformula(
    pij ~ 1 + (1 | Ref),
    phi ~ 1 + wij + (1 | Ref)
  )
  
  ## priors
  wc_priors <- c(
    ### Grand mean prior
    prior(normal(logit(0.03), 0.75), class = "Intercept"),
    
    ### Precision fixed effects
    prior(normal(0, 1), class = "Intercept", dpar = "phi"),
    prior(normal(0, 1), class = "b", dpar = "phi"),
    
    ### Random-effect SDs
    prior(exponential(2), class = "sd"),
    prior(exponential(2), class = "sd", dpar = "phi")
  )
  
  ## family
  wc_family <- Beta(link = "logit",
                    link_phi = "log")
  
    ### prior predictive check
    prior_check <- brm(
      data    = wc_dat,
      formula = wc_formula,
      prior   = wc_priors,
      family  = wc_family,
      sample_prior = "only",
      chains  = 4,
      iter    = 2000,
      warmup  = 0,
      seed    = 32608
    )
    
    pp_prior <- posterior_predict(prior_check)
    
    ### visualize
    stack(as.data.frame(t(pp_prior[1:200,]))) %>%
      ggplot(aes(x = values, fill = ind)) +
      geom_histogram() + 
      theme(legend.position = "none")
    
  ## posterior sampling
  wc_model <- brm(
    data    = wc_dat,
    formula = wc_formula,
    prior   = wc_priors,
    family  = wc_family,
    chains  = 4,
    iter    = 251000,
    warmup = 1000,
    cores   = 8,
    seed    = 32608,
    control = list(adapt_delta = 0.995)
  )
  
    ### posterior predictive checks
    pp_post <- posterior_predict(wc_model)
    
      #### Overlay densities
      ppc_dens_overlay(
        y = wc_dat$pij,
        yrep = pp_post[1:200, ]
      )
  
  ## extract posterior distribution for conviction error rate
  post_mu_wc <- posterior_linpred(
    wc_model,
    transform = TRUE,
    re_formula = NULL
  )
  post_mu_wc <- as.vector(post_mu_wc)
  
    ### descriptives
    quantile(post_mu_wc, c(.01, .05, .5, .95, .99))
    mean(post_mu_wc)

#==============================================================================#

# 3.1.1.2. False confession prevalence within wrongful convictions ====
    
  ## Generate Period variable
  fc_dat$Period <- as.numeric(gsub("[0-9]{4}-","", fc_dat$DataYear))
  fc_dat$Period[is.na(fc_dat$Period)] <- fc_dat$StudyYear[is.na(fc_dat$Period)]
  
  ## Rescale p_est, introduce weight variable(s)
  fc_dat <- fc_dat %>%
    mutate(
      recency = min_max_scale(Period),
      pij = pmin(pmax(p_est, epsilon), 1 - epsilon),
      type_wt = ifelse(weight == "secondary", 1, 0.9),
      wj = recency * type_wt * (N / max(N)),  # total weight
    )
  
  ## Formula
  fc_formula <- brmsformula(
    pij ~ 1 + (1 | Ref),
    phi ~ 1 + wj + (1 | Ref)
  )
    
  ## Priors
  fc_priors <- c(
    ### Grand mean prior
    prior(normal(logit(0.15), 0.50), class = "Intercept"),
    
    ### Precision fixed effects
    prior(normal(0, 1), class = "Intercept", dpar = "phi"),
    prior(normal(0, 1), class = "b", dpar = "phi"),
    
    ### Random-effect SDs
    prior(exponential(2), class = "sd"),
    prior(exponential(2), class = "sd", dpar = "phi")
  )
  
  ## Family
  fc_family <- Beta(link = "logit",
                    link_phi = "log")
  
  ## Modeling false confession (abbreviated code)
  fc_model <- brm(
    formula = fc_formula,
    data = fc_dat,
    family = fc_family,
    prior = fc_priors,
    chains = 4,
    iter = 251000,
    warmup = 1000,
    cores = 8,
    control = list(adapt_delta = 0.995)
  )
  
  ## posterior predictive checks
  pp_post <- posterior_predict(fc_model)
  
    ### Overlay densities.
    ppc_dens_overlay(
      y = fc_dat$pij,
      yrep = pp_post[1:200, ]
    )

  ## extract posterior distribution for conviction error rate
  post_mu_fc <- posterior_linpred(
    fc_model,
    transform = TRUE,
    re_formula = NULL
  )
  post_mu_fc <- as.vector(post_mu_fc)
  
  quantile(post_mu_fc, c(.01, .05, .5, .95, .99))
  mean(post_mu_fc)

#==============================================================================#
  
# 3.1.1.3 Joint probability base rate
  
  draws_wc <- as_draws_df(wc_model)$b_Intercept %>%
    plogis() %>%
    sample(size = 1e6, replace = TRUE)
  
  draws_fc <- as_draws_df(fc_model)$b_Intercept %>%
    plogis() %>%
    sample(size = 1e6, replace = TRUE)
  
  ## Calculate FCWC base rate
  fcwc_base_rate <- draws_wc * draws_fc
  
  ## Descriptives for FCWC base rate
  mean(fcwc_base_rate)
  median(fcwc_base_rate)
  quantile(fcwc_base_rate, c(0.025, 0.975))
  
#==============================================================================#

# Critique 2 
  ## Mourtgos & Adams (2026) argue that their specificity / sensitivity
  ## Likelihood statistics (which are uniformly sampled for both FCWC and
  ## non-FCWC cases) are the probability of a "problematic" interrogation
  ## tactic, but if this were true they would have attempted to model
  ## the fact that these more problematic are more likely to be present in
  ## FCWC cases. Even if divorced from the intent of the law enforcement
  ## officers, the existing literature clearly shows that FCWC interrogations
  ## are, on average, much longer than non-FCWC interrogations (see: 
  ## Leo, 2004; Kassin et al., 2007). Longer interrogations necessarily
  ## provide greater opportunity (and arguably an impetus) for problematic
  ## interrogation tactics.

  ## This could be approached by modelling the number of interrogation tactics:
  ## Kaplan & Cutler (2021) Co-occurrences among interrogation tactics in actual
  ## criminal investigations. Psychology, Crime & Law, 28(1): 1-19.
    
  ## However, we instead opt to simulate the plausible range of latent 
  ## interrogation intensities for FCWC and non-FCWC separately based on
  ## previous research on the average interrogation length for each.
      
    ### We have already estimated the FCWC base rate
    pi_fcwc <- fcwc_base_rate
    
    ### plausible interrogation intensity parameters
    fcwc_studies <- data.frame(study = c("Leo_2004", "Redlich_2011"),
               mean = c(16.23, 3.07),
               sd = c(15.05, 4.21),
               n = c(44, 35))
    fcwc_studies <- fcwc_studies[which(fcwc_studies$study=="Redlich_2011"),]

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

    fcwc_mix <- mixture_dist(fcwc_studies, weighting = "sample_size")
    nonfcwc_mix <- mixture_dist(nonfcwc_studies, weighting = "sample_size")

    ### Monte Carlo to approximate likelihoods
      #### Goal: Define function to generate average likelihoods for each level of K
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
      
      # Simulate durations from study-level mixture
      rlnorm_mixture <- function(n, study_mix,
                                 lower = 0, upper = 48){
        
        component <- sample(seq_len(nrow(study_mix)),
          size = n,
          replace = TRUE,
          prob = study_mix$weight)
        
        sims <- numeric(n)
        
        for (j in seq_len(nrow(study_mix))) {
          ind <- which(component == j)
          if (length(ind) > 0) {
            sims[ind] <- rlnorm_trunc(n = length(ind),
              meanlog = study_mix$meanlog[j],
              sdlog = study_mix$sdlog[j],
              lower = lower,
              upper = upper)
          }
        }
        
        return(data.frame(duration = sims,
                          study = study_mix$study[component]))
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

      # Check estimates
      pooled_mean <- sum(nonfcwc_studies$n * nonfcwc_studies$mean) / sum(nonfcwc_studies$n)
      results <- posterior_fcwc_given_H_mixture(h = pooled_mean,
                                                pi = pi_fcwc,
                                                fcwc_mix = fcwc_mix,
                                                nonfcwc_mix = nonfcwc_mix,
                                                upper = 48)
      c(mean = mean(results$posterior),
        median = median(results$posterior),
        lower_95 = quantile(results$posterior, 0.025),
        upper_95 = quantile(results$posterior, 0.975))
      
      pooled_mean <- sum(fcwc_studies$n * fcwc_studies$mean) / sum(fcwc_studies$n)
      results <- posterior_fcwc_given_H_mixture(h = pooled_mean,
                                                pi = pi_fcwc,
                                                fcwc_mix = fcwc_mix,
                                                nonfcwc_mix = nonfcwc_mix,
                                                upper = 48)
      c(mean = mean(results$posterior),
        median = median(results$posterior),
        lower_95 = quantile(results$posterior, 0.025),
        upper_95 = quantile(results$posterior, 0.975))
      
      
      results <- posterior_fcwc_given_H_mixture(h = 6.277229,
                                                pi = pi_fcwc,
                                                fcwc_mix = fcwc_mix,
                                                nonfcwc_mix = nonfcwc_mix,
                                                upper = 48)
      c(mean = mean(results$posterior),
        median = median(results$posterior),
        lower_95 = quantile(results$posterior, 0.025),
        upper_95 = quantile(results$posterior, 0.975))
      
      results <- posterior_fcwc_given_H_mixture(h = 12.060216,
                                                pi = pi_fcwc,
                                                fcwc_mix = fcwc_mix,
                                                nonfcwc_mix = nonfcwc_mix,
                                                upper = 48)
      c(mean = mean(results$posterior),
        median = median(results$posterior),
        lower_95 = quantile(results$posterior, 0.025),
        upper_95 = quantile(results$posterior, 0.975))
      
      results <- posterior_fcwc_given_H_mixture(h = 33.963305,
                                                pi = pi_fcwc,
                                                fcwc_mix = fcwc_mix,
                                                nonfcwc_mix = nonfcwc_mix,
                                                upper = 48)
      c(mean = mean(results$posterior),
        median = median(results$posterior),
        lower_95 = quantile(results$posterior, 0.025),
        upper_95 = quantile(results$posterior, 0.975))

#==============================================================================#

    #### Visualize likelihood distributions for FCWC v. non-FCWC
      ##### Line Graph
        hours <- seq(0, 48, by = 0.5)
        
          L1_boot <- bootstrap_mixture_likelihood(
            hours = hours,
            study_mix = fcwc_mix,
            n_rep = 1000,
            upper = 48) %>%
            mutate(condition = "FCWC")
          
          L0_boot <- bootstrap_mixture_likelihood(
            hours = hours,
            study_mix = nonfcwc_mix,
            n_rep = 1000,
            upper = 48) %>%
            mutate(condition = "non-FCWC")
          
          likelihood_boot <- bind_rows(L1_boot, L0_boot)
          
          likelihood_hdi <- likelihood_boot %>%
            group_by(hour, condition) %>%
            summarise(
              hdi_lower = HDInterval::hdi(
                likelihood,
                credMass = 0.95)[1],
              hdi_upper = HDInterval::hdi(
                likelihood,
                credMass = 0.95)[2],
              .groups = "drop"
            )
        
          likelihood_point <- data.frame(
            hour = hours,
            L1 = sapply(hours,
                        mixture_likelihood_H,
                        study_mix = fcwc_mix,
                        upper = 48),
            L0 = sapply(hours,
                        mixture_likelihood_H,
                        study_mix = nonfcwc_mix,
                        upper = 48)) %>%
            pivot_longer(
              cols = c(L1, L0),
              names_to = "condition",
              values_to = "likelihood"
            ) %>%
            mutate(condition = recode(
              condition,
              L1 = "FCWC",
              L0 = "non-FCWC"))
          
          likelihood_plot <- likelihood_point %>%
            left_join(likelihood_hdi,
              by = c("hour", "condition"))

        ggplot(likelihood_plot, aes(x = hour, y = likelihood)) +
          geom_ribbon(aes(ymin = hdi_lower,
                          ymax = hdi_upper,
                          fill = condition),
                      alpha = 0.20,
                      color = NA) +
          geom_line(aes(linetype = condition,
                        color = condition),
            linewidth = 0.75) +
          scale_x_continuous(breaks = seq(0, 48, by = 2)) +
          scale_y_continuous(limits = c(0, NA),
                             expand = expansion(mult = c(0, 0.05))) +
          scale_color_manual(values = c("FCWC" = "#D81B60",
                                        "non-FCWC" = "#004D40")) +
          scale_fill_manual(values = c("FCWC" = "#D81B60",
                                       "non-FCWC" = "#004D40")) +
          labs(title = "Figure 2a.",
               subtitle = "Simulated likelihood functions for interrogation duration",
               caption = paste0("Lines depict mixtures of truncated log-normal distributions ",
                                "for FCWC and non-FCWC interrogation durations,\n",
                                "derived of multiple studies, weighted by sample size. Shaded regions ",
                                "represent 95% HDIs obtained through\nparametric bootstrap ",
                                "resampling of study-level duration distributions."),
               color = NULL,
               fill = NULL,
               linetype = NULL,
               y = expression(paste("Likelihood, " ~ italic(f),"(", h ~ "|" ~ Y, ")")),
               x = "Interrogation Hours") +
          guides(fill = "none") +
          theme_classic() +
          theme(text = element_text(size = 27, family = "serif"),
                axis.title.y = element_text(margin = margin(r = 10, unit = "pt")),
                axis.title.x = element_text(margin = margin(t = 10, unit = "pt")),
                plot.title = element_text(size = 27, margin = margin(b = 5, unit = "pt")),
            plot.subtitle = element_text(size = 23, margin = margin(b = 20, unit = "pt")),
            plot.caption = element_text(size = 17, hjust = 0, margin = margin(t = 20, unit = "pt")),
            legend.key.spacing.y = unit(1.0, "cm"))        

#==============================================================================#
      
      ## Posterior over time
      hours <- seq(0, 48, by = 0.5)
      
      posterior_summary <- do.call(rbind,
        lapply(hours, function(h) {
          tmp <- posterior_fcwc_given_H_mixture(
            h = h,
            pi = pi_fcwc,
            fcwc_mix = fcwc_mix,
            nonfcwc_mix = nonfcwc_mix,
            upper = 48)
          tmp$hour <- h
          tmp$iteration <- seq_len(nrow(tmp))
          tmp})) %>%
        filter(hour > 0, is.finite(posterior)) %>%
        group_by(hour) %>%
        summarise(posterior_median = median(posterior, na.rm = TRUE),
                  posterior_mean = mean(posterior, na.rm = TRUE),
                  hdi_lower = HDInterval::hdi(posterior, credMass = 0.95)[1],
                  hdi_upper = HDInterval::hdi(posterior, credMass = 0.95)[2],
                  lower_95 = quantile(posterior, probs = 0.025),
                  upper_95 = quantile(posterior, probs = 0.975),
                  .groups = "drop")
      
      thresholds <- c(0.10, 0.25, 0.50)
      find_crossing <- function(x, y, threshold){
        d <- y - threshold
        idx <- which(d[-length(d)] * d[-1] <= 0)
        if (length(idx) == 0){return(NA_real_)}
        i <- idx[1]
        approx(x = y[c(i, i + 1)],
               y = x[c(i, i + 1)],
               xout = threshold)$y}
      
      crossing_dat <- data.frame(
        threshold = thresholds,
        hour = sapply(
          thresholds,
          function(p){find_crossing(
              x = posterior_summary$hour,
              y = posterior_summary$posterior_median,
              threshold = p)}))

      ggplot(posterior_summary, aes(x = hour, y = posterior_median)) +
        geom_ribbon(aes(ymin = hdi_lower, ymax = hdi_upper),
          alpha = 0.25, fill = "#D81B60") +
        geom_line(linewidth = 1, color = "#D81B60") +
        geom_hline(data = crossing_dat,
                   aes(yintercept = threshold),
                   linetype = "dashed",
                   linewidth = 0.6) +
        geom_vline(data = crossing_dat,
                   aes(xintercept = hour),
                   linetype = "dotted",
                   linewidth = 0.6) +
        geom_point(data = crossing_dat,
                   aes(x = hour,
                       y = threshold),
                   inherit.aes = FALSE,
                   size = 3) +
        geom_text(data = crossing_dat,
                  aes(x = c(40,40,0),
                      y = threshold,
                      label = paste0("Pr(FCWC) = ", threshold)),
                  inherit.aes = FALSE,
                  hjust = 0,
                  vjust = -0.5,
                  family = "serif",
                  size = 5) +
        scale_x_continuous(breaks = seq(0, 48, by = 2),
                           limits = c(0, 48)) +
        scale_y_continuous(breaks = seq(0, 0.75, by = 0.15),
                           limits = c(0, 0.75)) +
        labs(x = "Interrogation Hours",
             y = expression(Pr(FCWC ~ "|" ~ H == h)),
             title = "Figure 2b.",
             subtitle = "Posterior probability of FCWC by interrogation duration",
             caption = "Line depicts posterior median; shaded region depicts the 95% HDI.") +
        theme_classic() + 
        theme(legend.position = "none",
              text = element_text(size = 27, family = "serif"),
              axis.text.x = element_text(size = 17),
              axis.title.y = element_text(margin = margin(r = 10, unit = "pt")),
              axis.title.x = element_text(margin = margin(t = 10, unit = "pt")),
              plot.title = element_text(size = 27,
                                        margin = margin(b = 5, unit = "pt")),
              plot.subtitle = element_text(size = 23, 
                                           margin = margin(b = 20, unit = "pt")),
              plot.caption = element_text(size = 17, hjust = 0, 
                                          margin = margin(t = 20, unit = "pt")),
              legend.key.spacing.y = unit(1.0, "cm"))
      
#==============================================================================#
  
      hours <- seq(1, 24, by = 1)
      bin_width <- 1
      
      ## Choose the size of the population/sample represented
      N_total <- 10000
      
      posterior_draws <- do.call(rbind, lapply(hours, function(h){
          tmp <- posterior_fcwc_given_H_mixture(h = h,
            pi = pi_fcwc,
            fcwc_mix = fcwc_mix,
            nonfcwc_mix = nonfcwc_mix,
            upper = 48)
          tmp$hour <- h
          tmp$iteration <- seq_len(nrow(tmp))
          tmp})) %>%
        filter(is.finite(L1),
               is.finite(L0),
               is.finite(pi),
               L1 >= 0,
               L0 >= 0,
               pi >= 0,
               pi <= 1)
      
      duration_count_dat <- posterior_draws %>%
        mutate(joint_fcwc = pi * L1 * bin_width,
               joint_nonfcwc = (1 - pi) * L0 * bin_width) %>%
        group_by(hour) %>%
        summarise(FCWC = mean(joint_fcwc, na.rm = TRUE),
                  `Non-FCWC` = mean(joint_nonfcwc, na.rm = TRUE),
                  .groups = "drop")
      
      duration_count_dat <- duration_count_dat %>%
        mutate(total_mass = FCWC + `Non-FCWC`,
               normalization_constant = sum(total_mass),
               FCWC = N_total * FCWC / normalization_constant,
               `Non-FCWC` = N_total * `Non-FCWC` / normalization_constant) %>%
        select(hour, FCWC, `Non-FCWC`) %>%
        pivot_longer(cols = c(FCWC, `Non-FCWC`),
                     names_to = "condition",
                     values_to = "expected_cases") %>%
        mutate(condition = factor(condition, levels = c("FCWC", "Non-FCWC")))
    
      facet_bar_dat <- duration_count_dat %>%
        group_by(hour) %>%
        mutate(proportion = expected_cases / sum(expected_cases)) %>%
        ungroup() %>%
        pivot_longer(cols = c(proportion, expected_cases),
                     names_to = "measure",
                     values_to = "value") %>%
        mutate(measure = factor(measure,
                                levels = c("expected_cases", "proportion"),
                                labels = c("i. Simulated expected count by interrogation hour (N = 10,000)",
                                           "ii. Simulated proportion within interrogation hour")),
          condition = factor(condition, levels = c("FCWC", "Non-FCWC")))
      
      expected_n_labels <- duration_count_dat %>%
        group_by(hour) %>%
        summarise(total_n = sum(expected_cases, na.rm = TRUE),
                  .groups = "drop") %>%
        mutate(measure = factor("i. Simulated expected count by interrogation hour (N = 10,000)",
                                levels = levels(facet_bar_dat$measure)),
               label = paste0("n = ", round(total_n)))
      label_offset <- 0.04 * max(expected_n_labels$total_n, na.rm = TRUE)
      
      
      ggplot(facet_bar_dat,
        aes(x = hour, y = value, fill = condition)) +
        geom_col(width = bin_width * 0.95) +
        geom_text(data = expected_n_labels,
                  aes(x = hour, 
                      y = total_n + label_offset, 
                      label = label),
                  inherit.aes = FALSE,
                  angle = 90,
                  hjust = 0,
                  vjust = 0.6,
                  size = 6,
                  family = "serif") +
        facet_wrap(~ measure,
                   ncol = 1,
                   scales = "free_y") +
        ggh4x::facetted_pos_scales(y = list(measure == "ii. Simulated proportion within interrogation hour" ~
              scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1),
                                 limits = c(0, 1.05),
                                 labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
                                 expand = expansion(mult = c(0, 0))))) +
        scale_x_continuous(breaks = seq(1, 24, by = 1),
                           limits = c(0, 25),
                           expand = expansion(mult = c(0, 0))) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
        scale_fill_manual(values = c("FCWC" = "#D81B60",
                                     "Non-FCWC" = "gray5"),
                          name = NULL) +
        labs(x = "Interrogation Hours",
             y = NULL, 
             caption = paste("The upper panel depicts the marginal frequency of each duration",
                             "in a simulated sample of 10,000 cases.\nThe lower panel shows the",
                             "FCWC/non-FCWC composition conditional on duration.")) +
        coord_cartesian(clip = "off") +
        theme_classic() +
        theme(text = element_text(size = 27, family = "serif"),
              axis.text.x = element_text(size = 17),
              axis.title.x = element_text(margin = margin(t = 10, unit = "pt")),
              strip.background = element_blank(),
              strip.text = element_text(size = 22, hjust = 0.5, margin = margin(b = 15)),
              strip.placement = "outside",
              plot.title = element_text(size = 27, margin = margin(b = 5, unit = "pt")),
              plot.subtitle = element_text(size = 23, margin = margin(b = 20, unit = "pt")),
              plot.caption = element_text(size = 17, hjust = 0, margin = margin(t = 20, unit = "pt")),
              legend.position = "bottom",
              panel.spacing = unit(1.2, "cm")) +
        ggh4x::force_panelsizes(rows = grid::unit(c(2.5, 1), "null"))      
      
#==============================================================================#
      
      ### Visualize posterior distributions at 1.6 and 16 hours across sensitivity and specificity
      library(interp)
      
      hours <- c(1.53, 3.07, 6.28, 12.06)
      n_boot <- 2000
      
      L1_boot <- bootstrap_mixture_likelihood(
        hours = hours,
        study_mix = fcwc_mix,
        n_rep = n_boot,
        lower = 0,
        upper = 48) %>%
        rename(L1 = "likelihood")
      
      L0_boot <- bootstrap_mixture_likelihood(
        hours = hours,
        study_mix = nonfcwc_mix,
        n_rep = n_boot,
        lower = 0,
        upper = 48) %>%
        rename(L0 = "likelihood")
      
      likelihood_boot <- L1_boot %>%
        inner_join(L0_boot, by = c("iteration", "hour")) %>%
        filter(is.finite(L1), is.finite(L0), L1 > 0, L0 > 0) %>%
        mutate(LR = L1 / L0, log_LR = log(LR),
               f_lab = paste0(hour, " hours"))
      
      fitted_likelihoods <- map_dfr(hours,
        function(h) {
          tmp <- posterior_fcwc_given_H_mixture(
            h = h,
            pi = pi_fcwc,
            fcwc_mix = fcwc_mix,
            nonfcwc_mix = nonfcwc_mix,
            lower = 0,
            upper = 48)
          data.frame(hour = h,
            L1 = tmp$L1[1],
            L0 = tmp$L0[1],
            LR = tmp$LR[1],
            posterior = mean(tmp$posterior,
                             na.rm = TRUE))}) %>%
        mutate(log_LR = log(LR),
          f_lab = paste0(hour, " hours"))
      
      likelihood_ranges <- likelihood_boot %>%
        group_by(hour) %>%
        summarise(L0_lower = quantile(L0, 0, na.rm = TRUE),
                  L0_upper = quantile(L0, 1, na.rm = TRUE),
                  L1_lower = quantile(L1, 0, na.rm = TRUE),
                  L1_upper = quantile(L1, 1, na.rm = TRUE),
                  .groups = "drop")
      
      grid_df <- likelihood_ranges %>%
        group_by(hour) %>%
        group_split() %>%
        map_dfr(function(x) {
            h <- x$hour[1]
            expand_grid(L0 = seq(x$L0_lower, x$L0_upper, length.out = 150),
                        L1 = seq(x$L1_lower, x$L1_upper, length.out = 150)) %>%
              mutate(hour = h,
                     LR = L1 / L0,
                     log_LR = log(LR),
                     f_lab = paste0(h, " hours"))})
      
      n_pi_eval <- 5000
      pi_eval <- quantile(pi_fcwc[is.finite(pi_fcwc) & pi_fcwc > 0 & pi_fcwc < 1],
        probs = (seq_len(n_pi_eval) - 0.5) / n_pi_eval,
        names = FALSE)
      
      log_LR_range <- range(
        c(grid_df$log_LR,
          likelihood_boot$log_LR,
          fitted_likelihoods$log_LR),
        finite = TRUE)
      
      log_LR_lookup <- seq(
        log_LR_range[1],
        log_LR_range[2],
        length.out = 5000)
      
      posterior_lookup <- vapply(
        log_LR_lookup,
        function(log_lr){
          LR <- exp(log_lr)
          posterior_draws <- (LR * pi_eval) / (LR * pi_eval + (1 - pi_eval))
          mean(posterior_draws, na.rm = TRUE)},
        FUN.VALUE = numeric(1))
      
      posterior_from_log_LR <- function(log_LR) {
        approx(x = log_LR_lookup,
               y = posterior_lookup,
               xout = log_LR,
               rule = 2)$y}
      
      grid_df <- grid_df %>%
        mutate(posterior = posterior_from_log_LR(log_LR),
               posterior_odds = posterior / (1 - posterior),
               posterior_log_odds = qlogis(posterior))
      
      likelihood_boot <- likelihood_boot %>%
        mutate(posterior = posterior_from_log_LR(log_LR),
               posterior_log_odds = qlogis(posterior))
      
      facet_levels <- paste0(hours, " hours")
      
      grid_df <- grid_df %>%
        mutate(f_lab = factor(f_lab, levels = facet_levels))
      
      likelihood_boot <- likelihood_boot %>%
        mutate(f_lab = factor(f_lab, levels = facet_levels))
      
      fitted_likelihoods <- fitted_likelihoods %>%
        mutate(f_lab = factor(f_lab, levels = facet_levels))
      
      posterior_breaks <- c(0.001, 0.005,
                            0.01, 0.025, 0.05, 0.075,
                            0.10, 0.25, 0.50, 0.75, 1.00)
      
      legend_labels <- paste0("(",head(posterior_breaks, -1),
                              ", ", tail(posterior_breaks, -1), "]")
      
      grid_df %>%
        group_by(f_lab) %>%
        summarise(min_posterior = min(posterior, na.rm = TRUE),
                  max_posterior = max(posterior, na.rm = TRUE))
      
      ggplot(grid_df, aes(x = L0, y = L1, z = posterior)) +
        geom_contour_filled(breaks = posterior_breaks, na.rm = TRUE) +
        geom_point(data = likelihood_boot,
                   aes(x = L0,
                       y = L1),
                   inherit.aes = FALSE,
                   alpha = 0.3,
                   size = 0.6,
                   color = "white") +
        geom_point(data = fitted_likelihoods,
                   aes(x = L0,
                       y = L1),
                   inherit.aes = FALSE,
                   shape = 21,
                   size = 4,
                   stroke = 1.1,
                   fill = "white",
                   color = "black") +
        labs(x = expression(L[0](H)*","~"\u00ACFCWC"),
             y = expression(L[1](H)*","~"FCWC"),
             fill = expression(Pr(FCWC ~ "|" ~ H == h))) +
        scale_fill_viridis_d(option = "magma",
                             direction = -1,
                             na.translate = FALSE,
                             na.value = "transparent",
                             labels = legend_labels,
                             drop = FALSE) +
        facet_wrap(~f_lab, scales = "free") +
        theme_classic() +
        theme(text = element_text(size = 27, family = "serif"),
              axis.title.y = element_text(margin = margin(r = 10, unit = "pt")),
              axis.title.x = element_text(margin = margin(t = 10, unit = "pt")),
              plot.title = element_text(size = 26,
                                        margin = margin(b = 5, unit = "pt")),
              plot.subtitle = element_text(size = 23, 
                                           margin = margin(b = 20, unit = "pt")),
              plot.caption = element_text(size = 17, hjust = 0, 
                                          margin = margin(t = 20, unit = "pt")),
              strip.text = element_text(size = 23,
                                        margin = margin(t = 10, b = 10, 
                                                        unit = "pt")),
              plot.margin = margin(30, 30, 30, 30, "pt"),
              legend.position = "bottom",
              legend.title = element_text(margin = margin(r = 20)))            
