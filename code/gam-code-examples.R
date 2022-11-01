## GAM example model fitting code

## Packages
pkgs <- c("readr", "ggplot2", "mgcv", "dplyr", "readxl", "patchwork",
          "tibble", "tidyr", "rio", "purrr", "vegan")
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE)

##------------------------------------------------------------------------------
## WAIS Divide Ice Core Atmospheric CO2 record
uri <- 'https://www.ncei.noaa.gov/pub/data/paleo/icecore/antarctica/wais2021co2.txt'
ice_core <- read_tsv(uri,
                     na = c('', 'NA', 'NaN'),
                     comment = '#',
                     col_types = 'dddddc') %>%
  rename(depth = depth_m, # renames some variables
         age = age_calBP,
         co2_raw = CO2_blank_corrected,
         co2 = CO2_blank_gravity_corrected,
         std_err = CO2_se_1sigma,
         notes = notes) %>%
  mutate(neg_age = -age, # negative of age for getting time in right direction
         wts = 1/std_err, # we have a std. error for each observation
         wts = wts / mean(wts, na.rm = TRUE)) # convert to weights for gam()

## Plot
ice_core %>%
  ggplot(aes(x = neg_age, y = co2)) +
  geom_line() +
  labs(x = 'Calibrated years BP',
       y = expression(CO[2] ~ (ppm)))

# ?plotmath

## Model
ctrl <- gam.control(nthreads = 3, maxit = 500)

m1 <- gam(co2 ~ s(neg_age, k = 200), data = ice_core, method = 'REML',
          family = gaussian, weights = wts, control = ctrl)

## model summary
summary(m1)

## plot the fitted smooth
draw(m1, n = 500, residuals = TRUE)

plot(m1, n = 500)

## look at model diagnostics
appraise(m1, method = 'simulate')

layout(matrix(1:4, ncol = 2, byrow = TRUE))
gam.check(m1)
layout(1)

## check basis dimension
k.check(m1) # looks good

## heavy tails and variance isn't constant - try a Gaussian LS model
m2 <- gam(list(co2 ~ s(neg_age, k = 200), # model for the mean
               ~ s(wts) + s(neg_age, k = 30)), # model for variance
          data = ice_core, method = 'REML',
          family = gaulss(), control = ctrl)

draw(m2, n = 500)
appraise(m2)
k.check(m2)

## still some issues with data distribution
## technically CO2 can;t be negative so we can try a Gamma
m3 <- gam(co2 ~ s(neg_age, k = 300), data = ice_core, method = 'REML',
          family = Gamma(link = 'log'), weights = wts, control = ctrl)

draw(m3, n = 500)
appraise(m3, method = "simulate")
k.check(m3)

## is wiggliness the same on average over the series? try adaptive smooth
m4 <- gam(co2 ~ s(neg_age, k = 200, bs = 'ad', m = 10),
          data = ice_core, method = 'REML',
          family = Gamma(link = 'log'), weights = wts, control = ctrl)

AIC(m1, m2, m3, m4)

draw(m4, n = 500, residuals = TRUE)

appraise(m4, method = 'simulate')

k.check(m4)

## ------------------------------------------------------------------------------
## Lake 227 example

## Load data
## download file http://bit.ly/palaesig-227 save in working directory
tf <- tempfile(fileext = ".xlsx", tmpdir = "./")
tf
download.file("https://bit.ly/palaesig-227", destfile = tf)
lake227 <- read_excel(tf)

## Peter higlighted Fuco, Allox, Lutzeax, Pheo_b, Bcarot
## take only those variables and year
vars <- c('YEAR', 'FUCO', 'ALLOX', 'LUTZEAX', 'PHEO_B', 'BCAROT')
lake227 <- lake227 %>% select(all_of(vars))

## want nice names
names(lake227) <- c('Year', 'Fucoxanthin', 'Alloxanthin', 'LuteinZeaxanthin',
                    'Pheophytinb', 'BetaCarotene')

## take data from 1943 onward - replicate Cottingahm et al
lake227 <- lake227 %>% filter(Year >= 1943)

## to long format for modeling
lake227 <- lake227 %>%
  gather(key = Pigment, value = Concentration, - Year) %>%
  mutate(Pigment = factor(Pigment), cYear = (Year - mean(Year)) / 1000)

ctrl <- gam.control(nthreads = 2) # set up some control parameters

## fit the first model - intercept only for power
## single smooth for the scale
## Pigment specific smooths for for mean
mtwlss0 <- gam(list(Concentration ~ s(Year, Pigment, bs = 'fs', k = 10), # mean
                    ~ 1, # power, intercept only
                    ~ s(Year, k = 10)), # scale
               data = lake227,
               family = twlss,
               optimizer = 'efs', # twlss can be more stable with a the Extended Fellner Schall fit
               method = 'REML', # REML smoothness selection
               control = ctrl)

## see pigment specific smooths in scale also work, lower `k`
mtwlss <- gam(list(Concentration ~ s(Year, Pigment, bs = 'fs', k = 10),
                   ~ 1,
                   ~ s(Year, Pigment, bs = 'fs', k = 6)),
              data = lake227,
              family = twlss, optimizer = 'efs', method = 'REML',
              control = ctrl)

AIC(mtwlss0, mtwlss)

summary(mtwlss)

appraise(mtwlss)

draw(mtwlss)

k.check(mtwlss)

## ------------------------------------------------------------------------------
## Braya So
## load braya so data set
braya <- read_table(url("http://bit.ly/brayaso"), skip = 84,
  col_names = c("Depth", "DepthUpper", "DepthLower", "Year", "YearYoung",
    "YearOld", "UK37"), col_types = "ddddddd") |>
  mutate(sampleInterval = YearYoung - YearOld) # add a variable for the amount
# of time per sediment sample

## label for plotting
braya_ylabel <- expression(italic(U)[37]^{italic(k)})

## plot
ggplot(braya, aes(x = Year, y = UK37)) +
  geom_line(colour = "grey") +
  geom_point() +
  labs(y = braya_ylabel, x = "Year CE")

## fit the model with a continuous time AR1 --- needs optim as this is not a stable fit!
## also needs k setting lower than default
braya.car1 <- gamm(UK37 ~ s(Year, k = 5), data = braya,
                   correlation = corCAR1(form = ~ Year),
                   method = "REML",
                   control = list(niterEM = 0, optimMethod = "BFGS",
                                  opt = "optim"))
## fit model using GCV
braya.gcv <- gam(UK37 ~ s(Year, k = 30), data = braya)
braya.reml <- gam(UK37 ~ s(Year, k = 30), data = braya, method = "REML")

## CAR(1) parameter
brayaPhi <- intervals(braya.car1$lme)$corStruct
## fails - fit is non-positive definite

N <- 300
# number of points at which to evaluate the smooth
## data to predict at
newBraya <- data_slice(braya.reml, Year = evenly(Year, n = N))
# or
# newBraya <- with(braya, data.frame(Year = seq(min(Year), max(Year),
#                                               length.out = N)))
## add predictions from GAMM + CAR(1) model
newBraya <- cbind(newBraya,
                  data.frame(predict(braya.car1$gam, newBraya,
                                     se.fit = TRUE)))
crit.t <- qt(0.975, df = df.residual(braya.car1$gam))
newBraya <- transform(newBraya,
                      upper = fit + (crit.t * se.fit),
                      lower = fit - (crit.t * se.fit))
## add GAM GCV results
fit_gcv <- predict(braya.gcv, newdata = newBraya, se.fit = TRUE)
crit.t <- qt(0.975, df.residual(braya.gcv))
newGCV <- data.frame(Year
                     = newBraya[["Year"]],
                     fit
                     = fit_gcv$fit,
                     se.fit = fit_gcv$se.fit)
newGCV <- transform(newGCV,
                    upper = fit + (crit.t * se.fit),
                    lower = fit - (crit.t * se.fit))
# bind on GCV results
newBraya <- rbind(newBraya, newGCV)
## Add indicator variable for model
newBraya <- transform(newBraya,
                      Method = rep(c("GAMM (CAR(1))", "GAM (GCV)"),
                                   each = N))

## plot CAR(1) and GCV fits
braya_fitted <- ggplot(braya, aes(y = UK37, x = Year)) +
  geom_point() +
  geom_ribbon(data = newBraya,
              mapping = aes(x = Year, ymax = upper, ymin = lower,
                            fill = Method),
              alpha = 0.3, inherit.aes = FALSE) +
  geom_line(data = newBraya,
            mapping = aes(y = fit, x = Year, colour = Method)) +
  labs(y = braya_ylabel, x = "Year CE") +
  scale_color_manual(values = c("#5e3c99", "#e66101")) +
  scale_fill_manual(values = c("#5e3c99", "#e66101")) +
  theme(legend.position = "right")
braya_fitted

## lets fit the proper model
braya_reml <- gam(UK37 ~ s(Year, k = 40), data = braya,
                  method = "REML",
                  weights = sampleInterval / mean(sampleInterval))
summary(braya_reml)

draw(braya_reml, n = 500)

appraise(braya_reml, method = "simulate")

k.check(braya_reml)

## posterior simulation --- posterior smooths
post_sm <- smooth_samples(braya_reml, term = "s(Year)", n = 20, seed = 42,
                          unconditional = TRUE)

## plot
draw(post_sm, alpha = 0.3, colour = 'steelblue') +
  geom_line(data = smooth_estimates(braya_reml, n = 400),
            mapping = aes(x = Year, y = est, group = NULL),
            lwd = 1)

## derivatives of the smooth
dydt <- derivatives(braya_reml, term = "s(Year)", type = "central",
                    interval = "simultaneous")
dydt

draw(dydt)

## ----------------------------------------------------------------------------
# Gulf of Taranto in the Ionian Sea (Taricco et al., 2016)
#
# for more info see https://fromthebottomoftheheap.net/2016/12/16/pangaea-r-open-palaeo-data/
foram <- read_rds("data/taranto-foram-18-o.rds")
foram

ylabel <- expression(delta^{18} * O ~ "[‰ VPDB]")
xlabel <- "Age [ka BP]"

ggplot(foram, aes(y = d18O, x = Age_kaBP)) +
    geom_path() +
    scale_x_reverse(sec.axis = sec_axis( ~ 1950 - (. * 1000),
    name = "Age [AD]")) +
    scale_y_reverse() +
    labs(y = ylabel, x = xlabel)

# add an age variable going in the right direction
foram <- foram |> mutate(Age = -Age_kaBP)

# fit the GAM
m <- gam(d18O ~ s(Age, k = 100, bs = "ad"),
  data = foram, method = "REML")

appraise(m, metho = "simulate")

summary(m)

# data to predict at
ds <- data_slice(m, Age = evenly(Age, n = 200)) |>
  mutate(Age_kaBP = -Age)
fv <- fitted_values(m, data = ds, scale = "response")
fv

ggplot(foram, aes(y = d18O, x = Age_kaBP)) +
    geom_point() +
    geom_ribbon(data = fv,
      mapping = aes(x = Age_kaBP, ymin = lower, ymax = upper),
      fill = "grey", colour = NA, alpha = 0.7, inherit.aes  = FALSE) +
    geom_path(data = fv, mapping = aes(x = Age_kaBP, y = fitted),
      inherit.aes = FALSE, size = 1) +
    scale_x_reverse(sec.axis = sec_axis( ~ 1950 - (. * 1000),
      name = "Age [AD]")) +
    scale_y_reverse() +
    labs(y = ylabel, x = xlabel)

# fit the GAM
m_t <- gam(d18O ~ s(Age, k = 100, bs = "ad"),
  data = foram, method = "REML", family = scat())

fv_t <- fitted_values(m_t, data = ds, scale = "response")
fv_t

# stack the two sets of predictions / models
fv_both <- bind_rows(fv, fv_t) |>
  mutate(distribution = rep(c("Gaussian", "SCAT"), each = 200))

ggplot(foram, aes(y = d18O, x = Age_kaBP)) +
    geom_point() +
    geom_ribbon(data = fv_both,
      mapping = aes(x = Age_kaBP, ymin = lower, ymax = upper,
      fill = distribution), colour = NA, alpha = 0.2, inherit.aes  = FALSE) +
    geom_path(data = fv_both,
      mapping = aes(x = Age_kaBP, y = fitted, colour = distribution),
      inherit.aes = FALSE, size = 1) +
    scale_x_reverse(sec.axis = sec_axis( ~ 1950 - (. * 1000),
      name = "Age [AD]")) +
    scale_y_reverse() +
    labs(y = ylabel, x = xlabel)

# some support for the SCAT
AIC(m, m_t)

# Woodbridge et al

# Load data
# download file http://bit.ly/palaesig-227 save in working directory
tf <- tempfile(fileext = ".xlsx", tmpdir = "./")
tf
download.file("https://bit.ly/3fohMO6", destfile = tf)

source(url("https://bit.ly/3DRk0z1"))

allpoll_list <- rio::import_list(tf)

allpoll_nested <- tibble(Site = names(allpoll_list), polldata = allpoll_list)

allpoll_nested

# create character vector of non-pollen variables to remove
non_pollen <- c("Sample",
                "Radiocarbon years B.P.",
                "EPD default [yrs.BP.]",
                "EPD [yrs.BP.]",
                "Fossilva [yrs.BP.]",
                "Sum")

##################################################################
#
# Wrangling palaeodata with the Tidyverse
# Part 5: Script to apply modelling functions to all pollen sites
# Calculating palynological richness
#
##################################################################

# How does palynological richness vary with time?
# A simple measure is the number of taxa per level (T), but N depends on the
# total pollen count, which also varies among levels.  We therefore use
# rarefaction to calculate the expected number of taxa (t) in a
# sample of fixed size (n) taken from a larger sample (N) containing (T) taxa.

# The rarefy function in the vegan package takes a matrix or data frame of taxon
# counts and returns a vector of rarefied number of taxa.

# We wrap the rarefy function in a wrapper function that also removes non-pollen
# variables and renames columns of the original data, as before.

library("vegan") # for rarefy function
library("purrr")

fun_rare <- function(x) {
   x %>% select(!contains(non_pollen)) %>%
   rename("Depth" = `Depth (cm)`, "Age_BP" = `Cal. yr. BP`) %>%
   summarise(Depth, Age_BP, richness=rarefy(round(.[, -(1:2)]), 200)) %>%
   filter(Age_BP < 9000)
}

# apply fun_rare to our nested list using mutate / map
rich <- allpoll_nested |>
  mutate(polldata, rare = map(polldata, ~ fun_rare(.x))) |>
  unnest(rare) |>
  select(-polldata)
## Ignore warnings about count < 200

# examine our handywork
rich

# plot
ggplot(rich, aes(Age_BP, richness, col = Site)) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 20000, by = 2000)) +
  theme(legend.position = "none")

tf <- tempfile(fileext = ".xlsx", tmpdir = "./")
tf
download.file("https://bit.ly/3sPIXEM", destfile = tf)

site_info <- read_excel(tf, sheet="Site_Info")
site_info

rich2 <- rich %>%
  left_join(site_info, by = c("Site" = "Site_code"))

rich2

# Plot richness by region
# Add a regression line to identify trends
ggplot(rich2, aes(Age_BP, richness, col = Site)) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_continuous(breaks = seq(0, 20000, by = 2000)) +
  facet_wrap(~Region) +
  labs(x = "Years BP") +
  theme(legend.position = "none")

rich2 <- rich2 |>
  mutate(Site_name = factor(Site_name), neg_age = -Age_BP)

# Let's fit those as smooths via a GAM.
# GAM has a random effect for Site_name modelling average richness
# differences between sites, plus a simple smooth of time for each Site
# which we model as a random smooth
ctrl <- gam.control(trace = TRUE)
m_rich <- bam(richness ~ s(Site_name, bs = "re") +
    s(neg_age, Site_name, k = 5, bs = "fs"),
  data = rich2, family = tw(), method = "fREML", control = ctrl,
  nthreads = 4, discrete = TRUE)

# plot the estimated smooths
draw(m_rich)

# want to predict for all sites over all time!
# so we create a data slice from the model, add on an Age variable
# do a join with the site info to get the meta data we need based on
# Site_name, and finally we rename the Site variable
ds_rich <- data_slice(m_rich, neg_age = evenly(neg_age, n = 200),
    Site_name = evenly(Site_name)) |>
  mutate(Age_BP = -neg_age) |>
  left_join(site_info |> select(Site_name, Site_no, Region)) |>
  rename(Site = Site_no)

# get model estimated values of richness on the response scale
fv_rich <- fitted_values(m_rich, data = ds_rich, scale = "response")

# plot it
ggplot(rich2, aes(Age_BP, richness, colour = Site_name)) +
  geom_line(alpha = 0.5) +
  geom_line(data = fv_rich, aes(x = Age_BP, y = fitted)) +
  geom_ribbon(data = fv_rich, aes(x = Age_BP, ymin = lower, ymax = upper,
  fill = Site_name),
    alpha = 0.1, inherit.aes = FALSE) +
  scale_x_continuous(breaks = seq(0, 20000, by = 2000)) +
  facet_wrap(~Region) +
  labs(x = "Years BP") +
  theme(legend.position = "none")

# Boom! 👊🎤⏬