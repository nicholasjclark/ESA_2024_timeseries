# Load packages
library(mvgam)           # Fit, interrogate and forecast DGAMs
library(tidyverse)       # Tidy and flexible data manipulation
library(ggplot2)         # Flexible plotting

# Set up plotting environment
theme_set(theme_classic(base_size = 15,
                        base_family = 'serif'))
myhist = function(...){
  geom_histogram(col = 'white',
                 fill = 'darkred', ...)
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

# Load the annual American kestrel, Falco sparverius, abundance 
# time series taken in British Columbia, Canada. These data have 
# been collected annually, corrected for changes in observer 
# coverage and detectability, and logged. They can be found in 
# the MARSS package
load(url('https://github.com/atsa-es/MARSS/raw/master/data/kestrel.rda'))
head(kestrel)


# Arrange the data into a long-format data.frame
regions <- c("BC",
             "Alb",
             "Sask")
model_data <- do.call(rbind,
                      lapply(
                        seq_along(regions),
                        function(x){
                          data.frame(year = kestrel[, 1],
                                     # Reverse the logging so that we deal 
                                     # directly with the detection-adjusted 
                                     # counts
                                     adj_count = exp(kestrel[, 1 + x]),
                                     region = regions[x])}
                        )
                      ) %>%
  # Add series and time indicators for mvgam modelling
  dplyr::mutate(yearfac = as.factor(year),
                region = as.factor(region),
                series = as.factor(region),
                time = year)

# Inspect modelling data structure
head(model_data)
dplyr::glimpse(model_data)
levels(model_data$series)

# Split the data into a training and a testing split
data_train <- model_data %>%
  dplyr::filter(year <= 2001)
data_test <- model_data %>%
  dplyr::filter(year > 2001)

# Plot all three time series together
plot_mvgam_series(data = data_train,
                  y = 'adj_count',
                  series = 'all') +
  theme_bw(base_size = 15,
           base_family = 'serif')

# Now plot features for just one series at a time
plot_mvgam_series(data = data_train,
                  newdata = data_test,
                  y = 'adj_count',
                  series = 1,
                  lines = FALSE) & # The & is used by patchwork 
                                   # to add themes to all sub-plots
  theme_bw(base_size = 15,
           base_family = 'serif')

plot_mvgam_series(data = data_train,
                  newdata = data_test,
                  y = 'adj_count',
                  series = 2) &
  theme_bw(base_size = 15,
           base_family = 'serif')

plot_mvgam_series(data = data_train,
                  newdata = data_test,
                  y = 'adj_count',
                  series = 3) &
  theme_bw(base_size = 15,
           base_family = 'serif')

# Look at the distribution of the outcome variable
ggplot(data_train,
       aes(x = adj_count)) +
  myhist() +
  labs(x = 'Adjusted count', y = '') +
  hist_theme()
summary(data_train$adj_count)

# Heavy-ish tail to the right; perhaps a Gamma distribution
# Inspect default priors for a simple model that only includes
# random intercepts for years, implying that all three series 
# share the same (random) year effects
?mgcv::random.effects
?mgcv::gam.models
?mvgam::get_mvgam_priors
p <- get_mvgam_priors(formula = adj_count ~ series +
                        s(yearfac, bs = 're'),
                      data = data_train,
                      family = Gamma())
View(p)

# prior() from brms can be used within mvgam()
# to change default prior distributions
?brms::prior

# Fit the model
mod1 <- mvgam(
  # Observation formula
  formula = adj_count ~ series +
                s(yearfac, bs = 're'),

  # Updated prior distributions using brms::prior()
  priors = c(prior(std_normal(), class = b),
             prior(exponential(2), class = sigma_raw)),

  # Training data in mvgam's long format
  data = data_train,

  # Testing data in mvgam's long format
  newdata = data_test,

  # Gamma observation model with shared shape parameter
  family = Gamma(),
  share_obs_params = TRUE,
  
  # Ensure all messages are reported for transparency
  silent = 1
  )

# Look at the structure of the object
str(mod1, max.level = 1)
?mvgam::`mvgam-class`
methods(class = "mvgam")

# Look at the Stan code to better understand the model
stancode(mod1)

# Generate a methods skeleton with references
how_to_cite(mod1)

# Diagnostics
summary(mod1)
summary(mod1,
        include_betas = FALSE,
        smooth_test = FALSE)
mcmc_plot(mod1,
          type = 'rhat_hist')
mcmc_plot(mod1,
          variable = 'obs_params',
          type = 'trace')
mcmc_plot(mod1,
          variable = c('mean', 'sd'),
          regex = TRUE,
          type = 'trace')
pairs(mod1,
      variable = c('mean', 'sd'),
      regex = TRUE)

# Inspect estimated effects on the outcome scale ...
conditional_effects(mod1)

# ... and on the link scale
conditional_effects(mod1,
                    type = 'link')

# Many other types of predictions and contrasts can be 
# made with marginaleffects
marginaleffects::avg_predictions(mod1, variable = 'series')

# Unconditional posterior predictive checks to 
# look at model fit
pp_check(mod1,
         type = "ecdf_overlay_grouped",
         group = "series",
         ndraws = 50)
pp_check(mod1,
         type = "dens_overlay_grouped",
         group = "series",
         ndraws = 50)

# Conditional posterior predictive checks
hcs <- hindcast(mod1)
class(hcs)
?mvgam::`mvgam_forecast-class`
methods(class = "mvgam_forecast")

layout(matrix(1:4, nrow = 2, byrow = TRUE))
plot(hcs, series = 1)
plot(hcs, series = 2)
plot(hcs, series = 3)
layout(1)

# Residual checks
plot(mod1, type = 'residuals', series = 1)
plot(mod1, type = 'residuals', series = 2)
plot(mod1, type = 'residuals', series = 3)

# Inspect forecasts, which were already computed by the
# model for the test data
fcs <- forecast(mod1)
class(fcs)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
plot(fcs, series = 1)
plot(fcs, series = 2)
plot(fcs, series = 3)
layout(1)

# Another way to look at forecasts for this example
plot_predictions(mod1, newdata = model_data,
                 by = c('yearfac', 'series', 'series'))

# As a quick aside, here is a model with splines of year 
# for each region, which are partially pooled toward
# a shared spline of year
?mgcv::factor.smooth.interaction

mod1.2 <- mvgam(
  # Observation formula containing region-level intercepts and
  # hierarchical splines of year for each region
  formula = adj_count ~ series +
    
    # Shared smooth of year for all series
    s(year, k = 30, bs = 'cr') +
    
    # Deviation smooths for each series
    s(year, series, k = 10, bs = 'sz'),

  # Updated prior distributions using brms::prior()
  priors = prior(std_normal(),
                 class = b),

  # Training and testing data in mvgam's long format
  data = data_train,
  newdata = data_test,

  # Gamma observation model
  family = Gamma(),
  share_obs_params = TRUE
)

# Draw the individual component smooths
gratia::draw(mod1.2)

# All in-sample, unconditional plots look excellent!
plot_predictions(mod1.2,
                 by = c('year', 'series', 'series'),
                 points = 0.5)
pp_check(mod1.2,
         type = "dens_overlay_grouped",
         group = "series",
         ndraws = 50)
pp_check(mod1.2,
         type = "pit_ecdf_grouped",
         group = "series",
         ndraws = 50)

# And this model is slightly favoured using loo in-sample
# model checks
loo_compare(mod1, mod1.2)

# But what about forecasts?!?
plot_predictions(mod1.2, 
                 newdata = model_data,
                 by = c('year', 'series', 'series'),
                 type = 'response') +
  geom_vline(xintercept = 2000,
             linetype = 'dashed')

fcs <- forecast(mod1.2)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
plot(fcs, series = 1)
plot(fcs, series = 2)
plot(fcs, series = 3)
layout(1)

# Expand to a State-Space model with more appropriate nonlinear 
# temporal effects; here we use Gaussian Processes for the 
# shared and deviation effects, which tend to extrapolate
# much better than splines do
?brms::gp
?mvgam::AR
?mvgam::mvgam_formulae
p <- get_mvgam_priors(formula = adj_count ~ series,
                      trend_formula = ~ 
                        gp(year, 
                           k = 32) + 
                        gp(year, 
                           by = trend,
                           k = 20) - 1,
                      trend_model = AR(cor = TRUE),
                      data = model_data,
                      family = Gamma())
View(p)

# Fit the model
mod2 <- mvgam(
  # Observation formula, only containing region-level intercepts
  formula = adj_count ~ series,
  
  # Process model formula, containing hierarchical GPs of time
  trend_formula = ~ 
    gp(year, 
       k = 32, 
       cov = 'exponential') + 
    gp(year, 
       by = trend,
       k = 20, 
       cov = 'exponential') - 1,
  
  # Additional autoregressive dynamics (using a correlated AR(1))
  trend_model = AR(cor = TRUE),
  
  # Updated prior distributions using brms::prior()
  priors = c(prior(beta(3, 10),
                   class = sigma,
                   lb = 0, 
                   ub = 1),
             prior(std_normal(),
                   class = `alpha_gp_trend(year)`),
             prior(std_normal(),
                   class = `alpha_gp_trend(year):trendtrend1`),
             prior(std_normal(),
                   class = `alpha_gp_trend(year):trendtrend2`),
             prior(std_normal(),
                   class = `alpha_gp_trend(year):trendtrend3`),
             prior(normal(0.5, 0.2),
                   class = ar1),
             prior(std_normal(),
                   class = b)),
  
  # Training and testing data in mvgam's long format
  data = data_train,
  newdata = data_test,
  
  # Gamma observation model with independent shape parameters
  family = Gamma(),
  share_obs_params = FALSE,
  
  # Stan MCMC control for slower but more precise sampling
  control = list(adapt_delta = 0.95)
)

# Inspect the Stan code
stancode(mod2)
how_to_cite(mod2)

# Diagnostics
summary(mod2,
        include_betas = FALSE,
        smooth_test = FALSE)
how_to_cite(mod2)
mcmc_plot(mod2,
          type = 'rhat_hist')
mcmc_plot(mod2,
          variable = c('sigma',
                       'ar1',
                       'shape'),
          regex = TRUE,
          type = 'trace')

# Unconditional posterior check
pp_check(mod2,
         type = "dens_overlay_grouped",
         group = "series",
         ndraws = 50)

# Inferences and unconditional predictions
plot_predictions(mod2,
                 condition = c('year', 'series'),
                 type = 'link',
                 conf_level = 0.5)
plot_predictions(mod2,
                 condition = c('year', 'series', 'series'),
                 points = 0.5)
marginaleffects::avg_predictions(mod2,
                                 variable = 'series')

# Rate of change, averaged over all regions
patchwork::wrap_plots(
  plot_predictions(
    mod2, 
    by = 'year',
    newdata = datagrid(model = mod2,
                       year = seq(min(model_data$year),
                                  max(model_data$year),
                                  length.out = 75),
                       series = unique),
    type = 'link'
  ) +
    labs(y = 'Linear predictor',
         x = ''),
  plot_slopes(
    mod2,
    variables = 'year',
    by = 'year',
    newdata = datagrid(model = mod2,
                       year = seq(min(model_data$year),
                                  max(model_data$year),
                                  length.out = 75),
                       series = unique),
    type = 'link'
  ) +
    labs(y = 'Slope of linear predictor',
         x = 'Year') +
    geom_hline(yintercept = 0, linetype = 'dashed'),
  nrow = 2
) + 
  patchwork::plot_annotation(
    title = 'Average conditional and marginal effects'
    ) 

# Inspect the AR1 variance-covariance parameters
Sigma_pars <- matrix(NA,
                     nrow = 3,
                     ncol = 3)
for(i in 1:3){
  for(j in 1:3){
    Sigma_pars[i, j] <- paste0('Sigma[', i, ',', j, ']')
  }
}
mcmc_plot(mod2,
          variable = as.vector(t(Sigma_pars)),
          type = 'hist') +
  geom_vline(xintercept = 0,
             col = 'white',
             linewidth = 2) +
  geom_vline(xintercept = 0,
             linewidth = 1)

# Compare models using in-sample fit metrics
loo_compare(mod1,
            mod2)

# Look at forecasts from each model and compare
fcs1 <- forecast(mod1)
fcs2 <- forecast(mod2)

layout(matrix(1:2,
              nrow = 2,
              byrow = TRUE))
for(x in 1:3){
  plot(fcs1,
       series = x)
  title('Random effects of year')
  plot(fcs2,
       series = x)
  title('GPs of year with AR1 dynamics')
}

# Score forecasts
?mvgam::score
score(fcs1,
      score = 'energy')
score(fcs1,
      score = 'energy')$all_series$score
score(fcs2,
      score = 'energy')$all_series$score

score(fcs1,
      score = 'variogram')$all_series$score
score(fcs2,
      score = 'variogram')$all_series$score

# Now a completely different model that uses a State-Space Vector
# Autoregression of order 1, with only the region-level intercepts
# as regression parameters
varmod <- mvgam(
  # Observation formula, empty to only consider the Gamma 
  # observation process
  formula = adj_count ~ -1,

  # Process model formula that includes regional intercepts
  trend_formula = ~ trend,

  # A VAR(1) dynamic process with fully parameterized covariance 
  # matrix Sigma
  trend_model = VAR(cor = TRUE),

  # Modified prior distributions using brms::prior()
  priors = c(prior(std_normal(),
                   class = Intercept_trend),
             prior(std_normal(),
                   class = b),
             prior(beta(3, 10),
                   class = sigma,
                   lb = 0, 
                   ub = 1)),

  # The time series data in 'long' format
  data = data_train,
  newdata = data_test,

  # Gamma observation model with independent shape parameters
  family = Gamma(),
  share_obs_params = FALSE,

  # Stan MCMC control for slower but more precise sampling
  control = list(adapt_delta = 0.95)
)
summary(varmod)
how_to_cite(varmod)
mcmc_plot(varmod,
          type = 'rhat_hist')

# Estimates of the autoregressive coefficients
A_pars <- matrix(NA,
                 nrow = 3,
                 ncol = 3)
for(i in 1:3){
  for(j in 1:3){
    A_pars[i, j] <- paste0('A[', i, ',', j, ']')
  }
}
mcmc_plot(varmod,
          variable = as.vector(t(A_pars)),
          type = 'hist') +
  geom_vline(xintercept = 0,
             col = 'white',
             linewidth = 2) +
  geom_vline(xintercept = 0,
             linewidth = 1)

# And of the variance-covariance parameters
Sigma_pars <- matrix(NA,
                     nrow = 3,
                     ncol = 3)
for(i in 1:3){
  for(j in 1:3){
    Sigma_pars[i, j] <- paste0('Sigma[', i, ',', j, ']')
  }
}
mcmc_plot(varmod,
          variable = as.vector(t(Sigma_pars)),
          type = 'hist') +
  geom_vline(xintercept = 0,
             col = 'white',
             linewidth = 2) +
  geom_vline(xintercept = 0,
             linewidth = 1)

# Impulse response functions
?mvgam::irf
irfs <- irf(varmod,
            h = 12,
            orthogonal = FALSE)
plot(irfs, series = 1)
plot(irfs, series = 2)
plot(irfs, series = 3)

# Forecast error variance decompositions
?mvgam::fevd
fevds <- fevd(varmod, h = 12)
plot(fevds) +
  scale_fill_manual(values = c("#DCBCBC",
                               "#A25050",
                               "#5C0000")) +
  labs(fill = 'Process')
?mvgam::stability

# In-sample comparisons don't suggest much difference with mod2
loo_compare(mod1,
            mod2,
            varmod)

# What about forecast comparisons?
fcsvar <- forecast(varmod)

layout(matrix(1:2,
              nrow = 2,
              byrow = TRUE))
for(x in 1:3){
  plot(fcs2,
       series = x)
  title('GPs of year with AR1 dynamics')
  plot(fcsvar,
       series = x)
  title('VAR1 dynamics')
}

# Compare energy scores
score(fcs2,
      score = 'energy')$all_series$score
score(fcsvar,
      score = 'energy')$all_series$score

# Compare variogram scores
score(fcs2,
      score = 'variogram')$all_series$score
score(fcsvar,
      score = 'variogram')$all_series$score

# What about an ensemble of these two?
ens <- ensemble(fcs2, fcsvar, ndraws = 2000)
sum(score(fcs2,
          score = 'energy')$all_series$score)
sum(score(fcsvar,
          score = 'energy')$all_series$score)
sum(score(ens,
          score = 'energy')$all_series$score)

# Perhaps the VAR1 isn't capturing the nonlinear trends as well
# as the hierarchical GPs; but we could easily combine the two!
