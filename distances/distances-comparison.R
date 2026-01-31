library(tidyverse)
library(brms)
library(spdep)
library(sf)

d <- read_csv('./data/LV_freq_data.csv')

#africa_sf <- spData::world %>% filter(continent == "Africa")

## dist weights: IDW ##
dsf <- st_as_sf(d, coords = c("long", "lat"), crs = 4326)

# get the spatial extent of the languages with frequency data #
ids <- which(dsf$LVpct != 0)
dc <- st_crop(dsf, st_bbox(dsf[ids,]))
dc$LVnorm <- scale(dc$LVpct)
dc$long <- st_coordinates(dc)[,1]
dc$lat <- st_coordinates(dc)[,2]

dnear <- dnearneigh(dc, 0, 500) 
dnear2 <- dnearneigh(dc, 0, 1000) 

listexp <- nb2listwdist(dnear, dc, type="exp", style="W", 
                        alpha = 1, zero.policy=TRUE)
listexp2 <- nb2listwdist(dnear2, dc, type="exp", style="W", 
                        alpha = 1, zero.policy=TRUE)

W <- listw2mat(listexp)
W2 <- listw2mat(listexp2)

# brms models #
fit_gp <- brm(bf(
    LVnorm ~ 1 + 
    (1 | family_wals) +
    gp(long, lat, gr = TRUE)
  ),
  data = dc,
  family = student(),
  prior = c(
    prior(normal(0, 2), "Intercept")
    #prior(normal(0, 1), class = "sdgp"),
    #prior(inv_gamma(5, 5), class = "lscale")
  ),
  chains = 4,
  cores = 4,
  iter = 3000,
  warmup = 1500,
  seed = 42,
  control = list(adapt_delta = 0.9)
)

write_rds(fit_gp, './results/fit_gp_euclidean.rds')


fit_sar500 <- brm(
  bf(LVnorm ~ 1 + 
       (1|family_wals) +
       sar(W, type = "lag")
  ),
  data = dc,
  data2 = list(W = W),
  family = student(),
  prior = c(
    prior(normal(0, 1), "Intercept"),
    prior(normal(0, 0.5), class = "lagsar")
  ),
  chains = 4,
  cores = 4,
  iter = 3000,
  warmup = 1500,
  seed = 42,
  control = list(adapt_delta = 0.9)
)

write_rds(fit_sar500, './results/fit_exp_sar500.rds')

fit_sar1000 <- brm(
  bf(LVnorm ~ 1 + 
       (1|family_wals) +
       sar(W2, type = "lag")
  ),
  data = dc,
  data2 = list(W2 = W2),
  family = student(),
  prior = c(
    prior(normal(0, 2), "Intercept"),
    prior(normal(0, 0.5), class = "lagsar")
  ),
  chains = 4,
  cores = 4,
  iter = 3000,
  warmup = 1500,
  seed = 42,
  control = list(adapt_delta = 0.9)
)

fit_sar1000 

write_rds(fit_sar1000, './results/fit_exp_sar1000.rds')

cdist <- st_distance(dc)
clean_units <- function(x){
  attr(x,"units") <- NULL
  class(x) <- setdiff(class(x),"units")
  x
}
D <- clean_units(cdist)

# stan GP models
library(rstan)

group_ids <- as.integer(factor(dc$family_wals))
group_ids

model_gp <- rstan::stan_model('./stan/gp_re_model.stan')
y <- as.vector(dc$LVnorm)

gp_data <- list(N = nrow(dc),
                Y = y, 
                Dist = D,
                N_1 = length(unique(group_ids)),
                M_1 = 1, 
                J_1 = group_ids)

fit <- sampling(model_gp,
              gp_data,
              chains = 4,
              cores = 4,
              iter = 3000,
              warmup = 1500,
              seed = 42)
write_rds(fit, './fit_gp_circle.rds')
