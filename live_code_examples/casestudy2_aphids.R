# Load packages
library(mvgam)           # Fit, interrogate and forecast DGAMs
library(mgcv)            # Fit GAMs
library(tidyverse)       # Tidy and flexible data manipulation
library(ggplot2)         # Flexible plotting
library(gratia)          # Graceful plotting of GAM components

# Set up plotting environment
theme_set(theme_bw(base_size = 15,
                   base_family = 'serif'))
myhist = function(...){
  geom_histogram(col = 'white',
                 fill = '#B97C7C', ...)
}
hist_theme = function(){
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
}
par(family = "serif",
    las = 0,
    mar = c(4.2,
            4.4,
            2.2,
            2.2),
    mgp = c(2.2,
            0.5,
            0),
    bty = "l",
    cex.axis = 1.25,
    cex.lab = 1.4,
    cex.main = 1.5,
    xaxs = 'r',
    yaxs = 'r',
    pch = 16)

# Load the 'aphids' data from the ecostats package. In each of two fields 
# (one of oats, one of wheat) there were eight plots, four with plastic netting 
# to exclude birds, and four without. Aphid abundance was counted on seven 
# different occasions over the first 38 weeks following netting. The hypothesis
# was that aphid numbers would DECREASE when birds were excluded, because an important 
# food source to tree sparrows is aphid predators, hoverflies and ladybird beetles,
# so presence of birds may be limit the effectiveness of a biological control of aphids.
load(url('https://github.com/dwarton/ecostats/raw/main/data/aphids.RData'))

# Bind the two datasets (experimental observations of aphid abundances
# over time in oat and wheat crop plots, under two experimental treatments)
aphid_dat <- aphids$oat %>%
  dplyr::mutate(crop = 'oat',
                series = paste0('oat_plot_', Plot)) %>%
  dplyr::bind_rows(aphids$wheat %>%
                     dplyr::mutate(crop = 'wheat',
                                   series = paste0('wheat_plot_', Plot)))
# View the data structure
dplyr::glimpse(aphid_dat)

# Wrangle data to improve variable names and create a
# time_since_treat variable
aphid_dat %>%
  dplyr::mutate(series = as.factor(series),
                crop = as.factor(crop),
                birds_excluded = as.factor(ifelse(
                  Treatment == 'excluded', 'yes', 'no')),
                time_since_treat = Time) %>%
  janitor::clean_names() %>%
  dplyr::select(-plot, -logcount, -treatment) %>%
  dplyr::arrange(series, time) -> aphid_ts

# Plot the data
aphid_ts %>%
  ggplot(., aes(x = time_since_treat,
                y = counts,
                col = birds_excluded)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ series) +
  labs(x = 'Time since treatment',
       y = 'Counts of aphids')

# An mgcv model to get started
mod0 <- gam(
  # Observation formula
  formula = counts ~
    
    # Hierarchical (random) intercepts per individual plot
    s(series, bs = 're') +
    
    # Parametric interaction between crop type and bird exclusion,
    # used to centre the hierarchical smooth effects below
    crop * birds_excluded +
    
    # 'Average' nonlinear effect of time since treatment
    s(time_since_treat,
      k = 7,
      bs = 'bs') +
    
    # 'Deviation' nonlinear effects of time since treatment, where
    # every level of the interaction between crop type and bird exclusion
    # has a different smooth function that can deviate from the 'Average' 
    # function above. We use the 'sz' basis to ensure the deviations sum 
    # to zero, which is good for inference but requires that we centre the
    # effects using the parametric interaction effect above
    s(time_since_treat,
      birds_excluded, 
      crop,
      bs = 'sz', 
      k = 7),
  
  # The Aphid data in 'long' format
  data = aphid_ts,
  
  # A Negative Binomial observation family to capture excess
  # overdispersion
  family = nb()
)
summary(mod0)
draw(mod0)
appraise(mod0, method = 'simulate')

# Inspect fit against the observed data
plot_predictions(mod0,
                 by = c('time_since_treat',
                        'birds_excluded',
                        'series'),
                 points = 0.5)

# Average predictions for each treatment
plot_predictions(mod0,
                 by = c('time_since_treat',
                        'birds_excluded',
                        'crop'))

# Smoother predictions on a fine grid
plot_predictions(mod0,
                 by = c('time_since_treat',
                        'birds_excluded'),
                 newdata = datagrid(model = mod0,
                                    time_since_treat = 3:40,
                                    crop = 'oat',
                                    birds_excluded = unique)) +
  labs(title = 'Oat crops')

plot_predictions(mod0,
                 by = c('time_since_treat',
                        'birds_excluded'),
                 newdata = datagrid(model = mod0,
                                    time_since_treat = 30:70,
                                    crop = 'wheat',
                                    birds_excluded = unique))+
  labs(title = 'Wheat crops')

# It seems there might be some evidence to support the hypothesis that
# aphid numbers would DECREASE when birds were excluded (i.e. counts would be
# lower when birds_excluded == 'yes'). But this effect is not entirely clean. Is 
# our model accurately fitting the data well enough for us to make this conclusion?

# Inspect residuals for possible remaining autocorrelation
aphid_ts %>%
  dplyr::bind_cols(resids = residuals(mod0)) %>%
  ggplot(., aes(time_since_treat, 
                resids)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ series) +
  geom_hline(yintercept = 0, 
             linetype = 'dashed')

# There is clearly unmodelled residual autocorrelation for a lot of these plots.
# We will need an mvgam model, but how to deal with the non-discrete 'time' aspect?
head(aphid_ts, 12)
?mvgam::mvgam_trends
?mvgam::CAR

# As usual, a State-Space model is what we're after. But we will need to use a 
# Continuous-time autoregressive model (CAR model) to account for the irregular 
# temporal sampling
p <- get_mvgam_priors(
  # Observation formula, which only contains the plot-level 
  # hierarchical intercepts
  formula = counts ~ s(series, bs = 're'),
  
  # Process model formula, which contains the hierarchical smooth
  # functions and the crop x birds excluded interaction
  # effects, which will allow the deviation smooths to be zero-centred
  trend_formula = ~
    crop * birds_excluded +
    s(time_since_treat, 
      k = 7,
      bs = 'bs') +
    s(time_since_treat,
      birds_excluded, 
      crop,
      bs = 'sz', 
      k = 7) - 1,
  
  # A Continuous Time Autoregressive process for the latent dynamics
  trend_model = CAR(),
  
  # Aphid time series in 'long' format
  data = aphid_ts,
  
  # A Poisson observation model, which is simpler than the 
  # Negative Binomial and doesn't have an overdispersion parameter
  # that can compete with the CAR process
  family = poisson(),
)
View(p)

# The priors for the regression coefficients and the CAR variances
# are certainly too vague. But we can also use a more informative prior
# to reflect that we expect the AR1 params to be positive
mod1 <- mvgam(  
  # Observation formula, which only contains the plot-level 
  # hierarchical intercepts
  formula = counts ~ s(series, bs = 're'),
  
  # Process model formula, which contains the hierarchical smooth
  # functions and the crop x birds excluded interaction
  # effects, which will allow the deviation smooths to be zero-centred
  trend_formula = ~
    crop * birds_excluded +
    s(time_since_treat, 
      k = 7,
      bs = 'bs') +
    s(time_since_treat,
      birds_excluded, 
      crop,
      bs = 'sz', 
      k = 7) - 1,
  
  # A Continuous Time Autoregressive process for the latent dynamics,
  # estimated using the noncentred parameterization for improved sampler
  # efficiency
  trend_model = CAR(),
  noncentred = TRUE,
  
  # Updated priors using brms::prior()
  priors = c(prior(std_normal(),
                   class = b),
             prior(beta(1, 1),
                   class = ar1,
                   lb = 0,
                   ub = 1),
             prior(exponential(5),
                   class = sigma)),
  
  # Aphid time series in 'long' format
  data = aphid_ts,
  
  # A Poisson observation model, which is simpler than the 
  # Negative Binomial and doesn't have an overdispersion parameter
  # that can compete with the CAR process
  family = poisson(),
  control = list(max_treedepth = 11)
)

# Inspect the Stan code
stancode(mod1)

# Diagnostics
summary(mod1, include_betas = FALSE)
mcmc_plot(mod1, 
          type = 'rhat_hist')
mcmc_plot(mod1,
          variable = 'sigma',
          regex = TRUE,
          type = 'hist')
mcmc_plot(mod1,
          variable = 'ar1',
          regex = TRUE,
          type = 'hist')

# Inferences
draw(mod1)
draw(mod1, trend_effects = TRUE, parametric = TRUE)
conditional_effects(mod1)

# Plot unconditional predictions against the observed data to understand
# how the nonlinear effects have been estimated
plot_predictions(mod1,
                 by = c('time_since_treat',
                        'birds_excluded',
                        'series'),
                 points = 0.5)

# Unconditional posterior predictive checks
pp_check(mod1,
         type = "ecdf_overlay_grouped",
         group = "crop",
         ndraws = 50)
pp_check(mod1,
         type = "dens_overlay_grouped",
         group = "birds_excluded",
         ndraws = 50)

# Average predictions for each treatment
plot_predictions(mod1,
                 by = c('time_since_treat',
                        'birds_excluded',
                        'crop'),
                 type = 'expected')

# Smoother predictions on a fine grid
plot_predictions(mod1,
                 by = c('time_since_treat',
                        'birds_excluded'),
                 newdata = datagrid(model = mod1,
                                    time_since_treat = 3:40,
                                    crop = 'oat',
                                    birds_excluded = unique),
                 type = 'expected') +
  labs(title = 'Oat crops')

plot_predictions(mod1,
                 by = c('time_since_treat',
                        'birds_excluded'),
                 newdata = datagrid(model = mod1,
                                    time_since_treat = 30:70,
                                    crop = 'wheat',
                                    birds_excluded = unique),
                 type = 'expected') +
  labs(title = 'Wheat crops')

# Marginal effects at each timepoint
marginaleffects::avg_predictions(mod1,
                                 variables = c('time_since_treat',
                                               'crop', 'birds_excluded'),
                                 type = 'expected')

# Residuals
resids <- residuals(mod1)
aphid_ts %>%
  dplyr::bind_cols(med = resids[,1],
                   lower = resids[,3],
                   upper = resids[,4]) %>%
  ggplot(., aes(time_since_treat, med)) +
  geom_ribbon(aes(ymin = lower,
                  ymax = upper),
              alpha = 0.3) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ series) +
  geom_hline(yintercept = 0, linetype = 'dashed')

plot(mod1, type = 'forecast', series = 16)

# We can decompose the latent process model into the contributions from the
# hierarchical smooth functions and the contribution from the CAR(1)
# process by using unconditional predictions (that ignore the CAR(1)) 
# and comparing them to conditional hindcasts (that used the CAR(1))
smooth_preds <- predict(mod1, summary = FALSE,
                        type = 'link',
                        process_error = FALSE)
trend_preds <- do.call(cbind, 
                       hindcast(mod1, 
                                type = 'link')$hindcasts)

# This plot will show the unconditional predictions in grey (i.e. the overall
# effect of the hierarchical smooth functions) together with the conditional 
# predictions in red (which use the CAR(1) process as well)
ggplot(aphid_ts %>%
         dplyr::mutate(pred = apply(smooth_preds, 2, mean),
                       upper = apply(smooth_preds, 2, function(x) 
                         quantile(x, probs = 0.975)),
                       lower = apply(smooth_preds, 2, function(x) 
                         quantile(x, probs = 0.025)),
                       predcar = apply(trend_preds, 2, mean),
                       uppercar = apply(trend_preds, 2, function(x) 
                         quantile(x, probs = 0.975)),
                       lowercar = apply(trend_preds, 2, function(x) 
                         quantile(x, probs = 0.025))),
       aes(x = time, y = pred)) +
  geom_ribbon(aes(ymax = upper,
                  ymin = lower),
              alpha = 0.2) +
  geom_line() +
  geom_ribbon(aes(ymax = uppercar,
                  ymin = lowercar),
              fill = 'darkred',
              alpha = 0.1) +
  geom_line(aes(y = predcar),
            col = 'darkred',
            alpha = 0.5) +
  facet_wrap(~ series, scales = 'free') +
  labs(y = 'Posterior linear predictions (log scale)',
       x = '',
       title = 'Hierarchical smooths (black) + CAR(1)s (red)')
