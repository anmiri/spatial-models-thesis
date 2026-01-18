library(rstan)
library(sf)
library(tidyverse)
library(posterior)
library(loo)
library(sf)

## prep the original data ##
all_langs = read_rds("./data/df-wals.rds")

all_langs$Latitude = st_coordinates(all_langs$geometry)[,2]
all_langs$Longitude = st_coordinates(all_langs$geometry)[,1]

afr_langs <- all_langs %>% filter(Macroarea == "Africa")
eur_langs <- all_langs %>% filter(Macroarea == "Eurasia")
nam_langs <- all_langs %>% filter(Macroarea == "North America")
sam_langs <- all_langs %>% filter(Macroarea == "South America")
pap_langs <- all_langs %>% filter(Macroarea == "Multinesia")
aus_langs <- all_langs %>% filter(Macroarea == "Australia")

ids_afr <- which(all_langs$Glottocode %in% afr_langs$Glottocode)
ids_eur <- which(all_langs$Glottocode %in% eur_langs$Glottocode)
ids_nam <- which(all_langs$Glottocode %in% nam_langs$Glottocode)
ids_sam <- which(all_langs$Glottocode %in% sam_langs$Glottocode)
ids_pap <- which(all_langs$Glottocode %in% pap_langs$Glottocode)
ids_aus <- which(all_langs$Glottocode %in% aus_langs$Glottocode)

loc_afr = scale(cbind(afr_langs$Longitude, afr_langs$Latitude))
nAfr <- length(loc_afr[,1])
loc_eur = scale(cbind(eur_langs$Longitude, eur_langs$Latitude))
nEur = length(loc_eur[,1])
loc_nam = scale(cbind(nam_langs$Longitude, nam_langs$Latitude))
nNam <- length(loc_nam[,1])
loc_sam = scale(cbind(sam_langs$Longitude, sam_langs$Latitude))
nSam = length(loc_sam[,1])
loc_pap = scale(cbind(pap_langs$Longitude, pap_langs$Latitude))
nPap <- length(loc_pap[,1])
loc_aus = scale(cbind(aus_langs$Longitude, aus_langs$Latitude))
nAus = length(loc_aus[,1])

## load the fitted model and extract parameters ##
fit_gp_phylo = read_rds("./results/fit_hurdle_approx_gp_phy_3.rds")

hurdle_draws = as_draws(fit_gp_phylo)
alpha_draws <- subset_draws(hurdle_draws, variable = c("alpha"))
alpha_df = summarise_draws(alpha_draws)

alpha_afr <- alpha_df[ids_afr,]$mean
alpha_eur <- alpha_df[ids_eur,]$mean
alpha_nam <- alpha_df[ids_nam,]$mean
alpha_sam <- alpha_df[ids_sam,]$mean
alpha_pap <- alpha_df[ids_pap,]$mean
alpha_aus <- alpha_df[ids_aus,]$mean

lscale_draws = subset_draws(hurdle_draws,variable=c("lscale"), regex=FALSE)
lscale_df = summarise_draws(lscale_draws)

lscale_afr = lscale_df$median[1]
lscale_eur = lscale_df$median[2]
lscale_nam = lscale_df$median[3]
lscale_sam = lscale_df$median[4]
lscale_pap = lscale_df$median[5]
lscale_aus = lscale_df$median[6]

hurdle_lscale_draws = subset_draws(hurdle_draws,variable=c("hurdle_lscale"), regex=FALSE)
hurdle_lscale_df = summarise_draws(hurdle_lscale_draws)

hurdle_lscale_afr = hurdle_lscale_df$median[1]
hurdle_lscale_eur = hurdle_lscale_df$median[2]
hurdle_lscale_nam = hurdle_lscale_df$median[3]
hurdle_lscale_sam = hurdle_lscale_df$median[4]
hurdle_lscale_aus = hurdle_lscale_df$median[5]
hurdle_lscale_pap = hurdle_lscale_df$median[6]

sdgp_draws = subset_draws(hurdle_draws,variable=c("sdgp"), regex=FALSE)
sdgp_df = summarise_draws(sdgp_draws)
sdgp_afr = sdgp_df$median[1]
sdgp_eur = sdgp_df$median[2]
sdgp_nam = sdgp_df$median[3]
sdgp_sam = sdgp_df$median[4]
sdgp_pap = sdgp_df$median[5]
sdgp_aus = sdgp_df$median[6]

hurdle_sdgp_draws = subset_draws(hurdle_draws,variable=c("hurdle_sdgp"), regex=FALSE)
hurdle_sdgp_df = summarise_draws(hurdle_sdgp_draws)
hurdle_sdgp_afr = hurdle_sdgp_df$median[1]
hurdle_sdgp_eur = hurdle_sdgp_df$median[2]
hurdle_sdgp_nam = hurdle_sdgp_df$median[3]
hurdle_sdgp_sam = hurdle_sdgp_df$median[4]
hurdle_sdgp_pap = hurdle_sdgp_df$median[5]
hurdle_sdgp_aus = hurdle_sdgp_df$median[6]

shape_draws = subset_draws(hurdle_draws,variable=c("shape"), regex=TRUE)
shape_df = summarise_draws(shape_draws)
shape = shape_df$median

beta_afr_draws = subset_draws(hurdle_draws,variable=c("beta_afr"), regex=FALSE)
beta_afr_df = summarise_draws(beta_afr_draws)
beta_afr = beta_afr_df$median

beta_eur_draws = subset_draws(hurdle_draws,variable=c("beta_eur"), regex=FALSE)
beta_eur_df = summarise_draws(beta_eur_draws)
beta_eur = beta_eur_df$median

beta_nam_draws = subset_draws(hurdle_draws,variable=c("beta_nam"), regex=FALSE)
beta_nam_df = summarise_draws(beta_nam_draws)
beta_nam = beta_nam_df$median

beta_sam_draws = subset_draws(hurdle_draws,variable=c("beta_sam"), regex=FALSE)
beta_sam_df = summarise_draws(beta_sam_draws)
beta_sam = beta_sam_df$median

beta_pap_draws = subset_draws(hurdle_draws,variable=c("beta_pap"), regex=FALSE)
beta_pap_df = summarise_draws(beta_pap_draws)
beta_pap = beta_pap_df$median

beta_aus_draws = subset_draws(hurdle_draws,variable=c("beta_aus"), regex=FALSE)
beta_aus_df = summarise_draws(beta_aus_draws)
beta_aus = beta_aus_df$median

hurdle_beta_afr_draws = subset_draws(hurdle_draws,variable=c("hurdle_beta_afr"), regex=FALSE)
hurdle_beta_afr_df = summarise_draws(hurdle_beta_afr_draws)
hurdle_beta_afr = hurdle_beta_afr_df$median

hurdle_beta_eur_draws = subset_draws(hurdle_draws,variable=c("hurdle_beta_eur"), regex=FALSE)
hurdle_beta_eur_df = summarise_draws(hurdle_beta_eur_draws)
hurdle_beta_eur = hurdle_beta_eur_df$median

hurdle_beta_nam_draws = subset_draws(hurdle_draws,variable=c("hurdle_beta_nam"), regex=FALSE)
hurdle_beta_nam_df = summarise_draws(hurdle_beta_nam_draws)
hurdle_beta_nam = hurdle_beta_nam_df$median

hurdle_beta_sam_draws = subset_draws(hurdle_draws,variable=c("hurdle_beta_sam"), regex=FALSE)
hurdle_beta_sam_df = summarise_draws(hurdle_beta_sam_draws)
hurdle_beta_sam = hurdle_beta_sam_df$median

hurdle_beta_sam_draws = subset_draws(hurdle_draws,variable=c("hurdle_beta_sam"), regex=FALSE)
hurdle_beta_sam_df = summarise_draws(hurdle_beta_sam_draws)
hurdle_beta_sam = hurdle_beta_sam_df$median

hurdle_beta_pap_draws = subset_draws(hurdle_draws,variable=c("hurdle_beta_pap"), regex=FALSE)
hurdle_beta_pap_df = summarise_draws(hurdle_beta_pap_draws)
hurdle_beta_pap = hurdle_beta_pap_df$median

hurdle_beta_aus_draws = subset_draws(hurdle_draws,variable=c("hurdle_beta_aus"), regex=FALSE)
hurdle_beta_aus_df = summarise_draws(hurdle_beta_aus_draws)
hurdle_beta_aus = hurdle_beta_aus_df$median

hurdle_intercept_draws = subset_draws(hurdle_draws, variable=c("hurdle_intercept"),regex=TRUE)
hurdle_intercept_df = summarise_draws(hurdle_intercept_draws)
hurdle_intercept = hurdle_intercept_df$median

intercept_draws = subset_draws(hurdle_draws, variable=c("Intercept"),regex=TRUE)
intercept_df = summarise_draws(intercept_draws)
intercept = intercept_df$median

mu_area_draws = subset_draws(hurdle_draws, variable=c("mu_area"),regex=FALSE)
mu_area_df = summarise_draws(mu_area_draws)

mu_area_afr = mu_area_df$median[1]
mu_area_eur = mu_area_df$median[2]
mu_area_nam = mu_area_df$median[3]
mu_area_sam = mu_area_df$median[4]
mu_area_pap = mu_area_df$median[5]
mu_area_aus = mu_area_df$median[6]

hu_mu_area_draws = subset_draws(hurdle_draws, variable=c("hurdle_mu_area"),regex=FALSE)
hu_mu_area_df = summarise_draws(hu_mu_area_draws)
hurdle_mu_area_afr = hu_mu_area_df$median[1]
hurdle_mu_area_eur = hu_mu_area_df$median[2]
hurdle_mu_area_nam = hu_mu_area_df$median[3]
hurdle_mu_area_sam = hu_mu_area_df$median[4]
hurdle_mu_area_pap = hu_mu_area_df$median[5]
hurdle_mu_area_aus = hu_mu_area_df$median[6]

sigma_intercept_draws = subset_draws(hurdle_draws, variable=c("sigma_intercept"),regex=FALSE)
sigma_intercept_df = summarise_draws(sigma_intercept_draws)
sigma_intercept = sigma_intercept_df$median[1]

hurdle_sigma_intercept_draws = subset_draws(hurdle_draws, variable=c("hurdle_sigma_intercept"),regex=TRUE)
hurdle_sigma_intercept_df = summarise_draws(hurdle_sigma_intercept_draws)
hurdle_sigma_intercept = hurdle_sigma_intercept_df$median[1]

# basis functions
m2 <- 22
m1 <- 22

m1*m2

indices <- matrix(NA, m1*m2, 2)

mm=0;
for (r in 1:m1){
  for (s in 1:m2){
    mm = mm+1
    indices[mm,] = c(r, s)
  }
}
indices

M_nD= m1*m2

## generate posterior predictions using original data ##

og_data_afr = list(N = nAfr,
                lonlat = loc_afr,
                max_lon = max(loc_afr[,1]),
                max_lat = max(loc_afr[,2]),
                M = c(m1,m2),
                M_nD = M_nD,
                indices = indices,
                c1 = 4,
                c2 = 4,
                lscale = lscale_afr,
                sdgp = sdgp_afr,
                shape = shape,
                hurdle_intercept = hurdle_intercept,
                hurdle_lscale=hurdle_lscale_afr,
                hurdle_sdgp=hurdle_sdgp_afr,
                Intercept = intercept,
                mu_area = mu_area_afr,
                hurdle_mu_area=hurdle_mu_area_afr,
                sigma_intercept = sigma_intercept,
                hurdle_sigma_intercept=hurdle_sigma_intercept,
                beta = beta_afr,
                hurdle_beta = hurdle_beta_afr,
                alpha = alpha_afr
)

og_data_eur = list(N = nEur,
               lonlat = loc_eur,
               max_lon = max(loc_eur[,1]),
               max_lat = max(loc_eur[,2]),
               M = c(m1,m2),
               M_nD = M_nD,
               indices = indices,
               c1 = 4,
               c2 = 4,
               lscale = lscale_eur,
               sdgp = sdgp_eur,
               shape = shape,
               hurdle_intercept = hurdle_intercept,
               hurdle_lscale=hurdle_lscale_eur,
               hurdle_sdgp=hurdle_sdgp_eur,
               Intercept = intercept,
               mu_area = mu_area_eur,
               hurdle_mu_area=hurdle_mu_area_eur,
               sigma_intercept = sigma_intercept,
               hurdle_sigma_intercept=hurdle_sigma_intercept,
               beta = beta_eur,
               hurdle_beta = hurdle_beta_eur,
               alpha = alpha_eur
)

og_data_nam = list(N = nNam,
                   lonlat = loc_nam,
                   max_lon = max(loc_nam[,1]),
                   max_lat = max(loc_nam[,2]),
                   M = c(m1,m2),
                   M_nD = M_nD,
                   indices = indices,
                   c1 = 4,
                   c2 = 4,
                   lscale = lscale_nam,
                   sdgp = sdgp_nam,
                   shape = shape,
                   hurdle_intercept = hurdle_intercept,
                   hurdle_lscale=hurdle_lscale_nam,
                   hurdle_sdgp=hurdle_sdgp_nam,
                   Intercept = intercept,
                   mu_area = mu_area_nam,
                   hurdle_mu_area=hurdle_mu_area_nam,
                   sigma_intercept = sigma_intercept,
                   hurdle_sigma_intercept=hurdle_sigma_intercept,
                   beta = beta_nam,
                   hurdle_beta = hurdle_beta_nam,
                   alpha = alpha_nam
)

og_data_sam = list(N = nSam,
                   lonlat = loc_sam,
                   max_lon = max(loc_sam[,1]),
                   max_lat = max(loc_sam[,2]),
                   M = c(m1,m2),
                   M_nD = M_nD,
                   indices = indices,
                   c1 = 4,
                   c2 = 4,
                   lscale = lscale_sam,
                   sdgp = sdgp_sam,
                   shape = shape,
                   hurdle_intercept = hurdle_intercept,
                   hurdle_lscale=hurdle_lscale_sam,
                   hurdle_sdgp=hurdle_sdgp_sam,
                   Intercept = intercept,
                   mu_area = mu_area_sam,
                   hurdle_mu_area=hurdle_mu_area_sam,
                   sigma_intercept = sigma_intercept,
                   hurdle_sigma_intercept=hurdle_sigma_intercept,
                   beta = beta_sam,
                   hurdle_beta = hurdle_beta_sam,
                   alpha=alpha_sam
)

og_data_pap = list(N = nPap,
                   lonlat = loc_pap,
                   max_lon = max(loc_pap[,1]),
                   max_lat = max(loc_pap[,2]),
                   M = c(m1,m2),
                   M_nD = M_nD,
                   indices = indices,
                   c1 = 4,
                   c2 = 4,
                   lscale = lscale_pap,
                   sdgp = sdgp_pap,
                   shape = shape,
                   hurdle_intercept = hurdle_intercept,
                   hurdle_lscale=hurdle_lscale_pap,
                   hurdle_sdgp=hurdle_sdgp_pap,
                   Intercept = intercept,
                   mu_area = mu_area_pap,
                   hurdle_mu_area=hurdle_mu_area_pap,
                   sigma_intercept = sigma_intercept,
                   hurdle_sigma_intercept=hurdle_sigma_intercept,
                   beta = beta_pap,
                   hurdle_beta = hurdle_beta_pap,
                   alpha=alpha_pap
)

og_data_aus = list(N = nAus,
                   lonlat = loc_aus,
                   max_lon = max(loc_aus[,1]),
                   max_lat = max(loc_aus[,2]),
                   M = c(m1,m2),
                   M_nD = M_nD,
                   indices = indices,
                   c1 = 4,
                   c2 = 4,
                   lscale = lscale_aus,
                   sdgp = sdgp_aus,
                   shape = shape,
                   hurdle_intercept = hurdle_intercept,
                   hurdle_lscale=hurdle_lscale_aus,
                   hurdle_sdgp=hurdle_sdgp_aus,
                   Intercept = intercept,
                   mu_area = mu_area_aus,
                   hurdle_mu_area=hurdle_mu_area_aus,
                   sigma_intercept = sigma_intercept,
                   hurdle_sigma_intercept=hurdle_sigma_intercept,
                   beta = beta_aus,
                   hurdle_beta = hurdle_beta_aus,
                   alpha=alpha_aus
)

gm_hu_gp_phy = rstan::stan_model("./stan/approx-gp-phylo-preds.stan") 

yrep_afr <- sampling(gm_hu_gp_phy, 
                       data   = og_data_afr, 
                       chains = 1,
                       cores  = 1,
                       warmup = 0,
                       iter   = 500,
                       seed   = 42,
                       pars   = c("preds"),
                       algorithm="Fixed_param")

yrep_afr_hu <- sampling(gm_hu_gp_phy, 
                     data   = og_data_afr, 
                     chains = 1,
                     cores  = 1,
                     warmup = 0,
                     iter   = 500,
                     seed   = 42,
                     pars   = c("hurdle_preds"),
                     algorithm="Fixed_param")


yrep_eur <- sampling(gm_hu_gp_phy, 
                     data   = og_data_eur, 
                     chains = 1,
                     cores  = 1,
                     warmup = 0,
                     iter   = 500,
                     seed   = 42,
                     pars   = c("preds"),
                     algorithm="Fixed_param")

yrep_eur_hu <- sampling(gm_hu_gp_phy, 
                     data   = og_data_eur, 
                     chains = 1,
                     cores  = 1,
                     warmup = 0,
                     iter   = 500,
                     seed   = 42,
                     pars   = c("hurdle_preds"),
                     algorithm="Fixed_param")

yrep_nam <- sampling(gm_hu_gp_phy, 
                     data   = og_data_nam, 
                     chains = 1,
                     cores  = 1,
                     warmup = 0,
                     iter   = 500,
                     seed   = 42,
                     pars   = c("preds"),
                     algorithm="Fixed_param")

yrep_nam_hu <- sampling(gm_hu_gp_phy, 
                     data   = og_data_nam, 
                     chains = 1,
                     cores  = 1,
                     warmup = 0,
                     iter   = 500,
                     seed   = 42,
                     pars   = c("hurdle_preds"),
                     algorithm="Fixed_param")

yrep_sam <- sampling(gm_hu_gp_phy, 
                     data   = og_data_sam, 
                     chains = 1,
                     cores  = 1,
                     warmup = 0,
                     iter   = 500,
                     seed   = 42,
                     pars   = c("preds"),
                     algorithm="Fixed_param")

yrep_sam_hu <- sampling(gm_hu_gp_phy,
                     data   = og_data_sam, 
                     chains = 1,
                     cores  = 1,
                     warmup = 0,
                     iter   = 500,
                     seed   = 42,
                     pars   = c("hurdle_preds"),
                     algorithm="Fixed_param")

yrep_pap <- sampling(gm_hu_gp_phy,
                     data   = og_data_pap, 
                     chains = 1,
                     cores  = 1,
                     warmup = 0,
                     iter   = 500,
                     seed   = 42,
                     pars   = c("preds"),
                     algorithm="Fixed_param")

yrep_pap_hu <- sampling(gm_hu_gp_phy,
                     data   = og_data_pap, 
                     chains = 1,
                     cores  = 1,
                     warmup = 0,
                     iter   = 500,
                     seed   = 42,
                     pars   = c("hurdle_preds"),
                     algorithm="Fixed_param")

yrep_aus <- sampling(gm_hu_gp_phy, 
                     data   = og_data_aus, 
                     chains = 1,
                     cores  = 1,
                     warmup = 0,
                     iter   = 500,
                     seed   = 42,
                     pars   = c("preds"),
                     algorithm="Fixed_param")

yrep_aus_hu <- sampling(gm_hu_gp_phy, 
                     data   = og_data_aus, 
                     chains = 1,
                     cores  = 1,
                     warmup = 0,
                     iter   = 500,
                     seed   = 42,
                     pars   = c("hurdle_preds"),
                     algorithm="Fixed_param")

# PP checks
y_afr = afr_langs$counts_wals
y_afr_hu = ifelse(y_afr > 0, 1, 0)
y_eur = eur_langs$counts_wals
y_eur_hu = ifelse(y_eur > 0, 1, 0)
y_nam = nam_langs$counts_wals
y_nam_hu = ifelse(y_nam > 0, 1, 0)
y_sam = sam_langs$counts_wals
y_sam_hu = ifelse(y_sam > 0, 1, 0)
y_pap = pap_langs$counts_wals
y_pap_hu = ifelse(y_pap > 0, 1, 0)
y_aus = aus_langs$counts_wals
y_aus_hu = ifelse(y_aus > 0, 1, 0)

yrep_draws_afr = as_draws_matrix(yrep_afr)
yrep_draws_afr = yrep_draws_afr[,-2339] # get rid of lp__ column

yrep_draws_eur_hu = as_draws_matrix(yrep_hu_eur)
yrep_draws_eur_hu = yrep_draws_eur_hu[,-1944] # get rid of lp__ column

yrep_draws_hu_nam = as_draws_matrix(yrep_nam_hu)
dim(yrep_draws_hu_nam)
yrep_draws_hu_nam = yrep_draws_hu_nam[,-742] # get rid of lp__ column

yrep_draws_hu_sam = as_draws_matrix(yrep_sam_hu)
dim(yrep_draws_hu_sam)
yrep_draws_hu_sam = yrep_draws_hu_sam[,-679] # get rid of lp__ column

yrep_draws_hu_pap = as_draws_matrix(yrep_pap_hu)
dim(yrep_draws_hu_pap)
yrep_draws_hu_pap = yrep_draws_hu_pap[,-2199] # get rid of lp__ column
yrep_draws_hu_pap

yrep_draws_hu_aus = as_draws_matrix(yrep_aus_hu)
dim(yrep_draws_hu_aus)
yrep_draws_hu_aus = yrep_draws_hu_aus[,-385] # get rid of lp__ column

ll <- extract_log_lik(fit_gp_phylo)

ll_afr <- ll[7500:7999,ids_afr]
ll_eur <- ll[7500:7999,ids_eur]
ll_nam <- ll[7500:7999,ids_nam]
ll_sam <- ll[7500:7999,ids_sam]
ll_aus <- ll[7500:7999,ids_aus]
ll_pap <- ll[7500:7999,ids_pap]

afr_bacc <- loo_predictive_metric(yrep_draws_hu_afr, y_afr_hu, 
                      ll_afr, metric = c("balanced_acc"))

eur_bacc <- loo_predictive_metric(yrep_draws_eur_hu, y_eur_hu, 
                                  ll_eur, metric = c("balanced_acc"))

nam_bacc <- loo_predictive_metric(yrep_draws_hu_nam, y_nam_hu, 
                                  ll_nam, metric = c("balanced_acc"))

sam_bacc <- loo_predictive_metric(yrep_draws_hu_sam, y_sam_hu, 
                                  ll_sam, metric = c("balanced_acc"))

pap_bacc <- loo_predictive_metric(yrep_draws_hu_pap, y_pap_hu, 
                                  ll_pap, metric = c("balanced_acc"))

aus_bacc <- loo_predictive_metric(yrep_draws_hu_aus, y_aus_hu, 
                                  ll_aus, metric = c("balanced_acc"))

all_bacc <- c(afr_bacc$estimate,
     eur_bacc$estimate,
     nam_bacc$estimate,
     sam_bacc$estimate,
     pap_bacc$estimate,
     aus_bacc$estimate)

yrep_draws_eur = as_draws_matrix(yrep_eur)
yrep_draws_eur = yrep_draws_eur[,-1944] # get rid of lp__ column

yrep_draws_nam = as_draws_matrix(yrep_nam)
yrep_draws_nam = yrep_draws_nam[,-742] # get rid of lp__ column

yrep_draws_sam = as_draws_matrix(yrep_sam)
yrep_draws_sam = yrep_draws_sam[,-679] # get rid of lp__ column

yrep_draws_pap = as_draws_matrix(yrep_pap)
yrep_draws_pap = yrep_draws_pap[,-2199] # get rid of lp__ column

yrep_draws_aus = as_draws_matrix(yrep_aus)
yrep_draws_aus = yrep_draws_aus[,-385] # get rid of lp__ column

ll <- extract_log_lik(fit_gp_phylo)

ll_afr <- ll[7500:7999,ids_afr]
ll_eur <- ll[7500:7999,ids_eur]
ll_nam <- ll[7500:7999,ids_nam]
ll_sam <- ll[7500:7999,ids_sam]
ll_aus <- ll[7500:7999,ids_aus]
ll_pap <- ll[7500:7999,ids_pap]

afr_rmse <- loo_predictive_metric(yrep_draws_afr, y_afr, 
                                  ll_afr, metric = c("rmse"))

eur_rmse <- loo_predictive_metric(yrep_draws_eur, y_eur, 
                                  ll_eur, metric = c("rmse"))

nam_rmse <- loo_predictive_metric(yrep_draws_nam, y_nam, 
                                  ll_nam, metric = c("rmse"))

sam_rmse <- loo_predictive_metric(yrep_draws_sam, y_sam, 
                                  ll_sam, metric = c("rmse"))

pap_rmse <- loo_predictive_metric(yrep_draws_pap, y_pap, 
                                  ll_pap, metric = c("rmse"))

aus_rmse <- loo_predictive_metric(yrep_draws_aus, y_aus, 
                                  ll_aus, metric = c("rmse"))

all_rmse <- c(afr_rmse$estimate,
              eur_rmse$estimate,
              nam_rmse$estimate,
              sam_rmse$estimate,
              pap_rmse$estimate,
              aus_rmse$estimate)
mean(all_rmse)


## draw spatial predictions using new locations ##

# sample new lon/lat coordinates along a grid for spatial interpolation
newcoords_afr <- tibble(
  longitude = rep(seq(from = -20, to = 60, length.out = 100), each = 100)
  , latitude = rep(seq(from = -35, to = 40, length.out = 100), 100)
)

newcoords_eur <- tibble(
  longitude = rep(seq(from = -180, to = 180, length.out = 100), each = 100)
  , latitude = rep(seq(from = -10, to = 81, length.out = 100), 100)
)

newcoords_sam <- tibble(
  longitude = rep(seq(from = -91, to = -34, length.out = 100), each = 100)
  , latitude = rep(seq(from = -56, to = 12, length.out = 100), 100)
)

newcoords_nam <- tibble(
  longitude = rep(seq(from = -180, to = 30, length.out = 100), each = 100)
  , latitude = rep(seq(from = 0, to = 80, length.out = 100), 100)
)

newcoords_aus <- tibble(
  longitude = rep(seq(from = 110, to = 160, length.out = 100), each = 100)
  , latitude = rep(seq(from = 0, to = -60, length.out = 100), 100)
)

newcoords_pap <- tibble(
  longitude = rep(seq(from = 100, to = 160, length.out = 100), each = 100)
  , latitude = rep(seq(from = -15, to = 15, length.out = 100), 100)
)

newlonlat_afr = as.matrix(unname(newcoords_afr))
newlonlat_afr[,1] = scale(newlonlat_afr[,1])
newlonlat_afr[,2] = scale(newlonlat_afr[,2])

newlonlat_eur = as.matrix(unname(newcoords_eur))
newlonlat_eur[,1] = scale(newlonlat_eur[,1])
newlonlat_eur[,2] = scale(newlonlat_eur[,2])

newlonlat_nam = as.matrix(unname(newcoords_nam))
newlonlat_nam[,1] = scale(newlonlat_nam[,1])
newlonlat_nam[,2] = scale(newlonlat_nam[,2])

newlonlat_sam = as.matrix(unname(newcoords_sam))
newlonlat_sam[,1] = scale(newlonlat_sam[,1])
newlonlat_sam[,2] = scale(newlonlat_sam[,2])

newlonlat_aus = as.matrix(unname(newcoords_aus))
newlonlat_aus[,1] = scale(newlonlat_aus[,1])
newlonlat_aus[,2] = scale(newlonlat_aus[,2])

newlonlat_pap = as.matrix(unname(newcoords_pap))
newlonlat_pap[,1] = scale(newlonlat_pap[,1])
newlonlat_pap[,2] = scale(newlonlat_pap[,2])

new_data_afr = list(N = nrow(newlonlat_afr),
                lonlat = newlonlat_afr,
                max_lon = max(newlonlat_afr[,1]),
                max_lat = max(newlonlat_afr[,2]),
                M = c(m1,m2),
                M_nD = M_nD,
                indices = indices,
                c1 = 4,
                c2 = 4,
                lscale = lscale_afr,
                sdgp = sdgp_afr,
                shape = shape,
                hurdle_intercept = hurdle_intercept,
                hurdle_lscale=hurdle_lscale_afr,
                hurdle_sdgp=hurdle_sdgp_afr,
                Intercept = intercept,
                mu_area = mu_area_afr,
                hurdle_mu_area=hurdle_mu_area_afr,
                sigma_intercept = sigma_intercept,
                hurdle_sigma_intercept=hurdle_sigma_intercept,
                beta = beta_afr,
                hurdle_beta = hurdle_beta_afr
)

new_data_eur = list(N = nrow(newlonlat_eur),
                    lonlat = newlonlat_eur,
                    max_lon = max(newlonlat_eur[,1]),
                    max_lat = max(newlonlat_eur[,2]),
                    M = c(m1,m2),
                    M_nD = M_nD,
                    indices = indices,
                    c1 = 4,
                    c2 = 4,
                    lscale = lscale_eur,
                    sdgp = sdgp_eur,
                    shape = shape,
                    hurdle_intercept = hurdle_intercept,
                    hurdle_lscale=hurdle_lscale_eur,
                    hurdle_sdgp=hurdle_sdgp_eur,
                    Intercept = intercept,
                    mu_area = mu_area_eur,
                    hurdle_mu_area=hurdle_mu_area_eur,
                    sigma_intercept = sigma_intercept,
                    hurdle_sigma_intercept=hurdle_sigma_intercept,
                    beta = beta_eur,
                    hurdle_beta = hurdle_beta_eur
)

new_data_nam = list(N = nrow(newlonlat_nam),
                    lonlat = newlonlat_nam,
                    max_lon = max(newlonlat_nam[,1]),
                    max_lat = max(newlonlat_nam[,2]),
                    M = c(m1,m2),
                    M_nD = M_nD,
                    indices = indices,
                    c1 = 4,
                    c2 = 4,
                    lscale = lscale_nam,
                    sdgp = sdgp_nam,
                    shape = shape,
                    hurdle_intercept = hurdle_intercept,
                    hurdle_lscale=hurdle_lscale_nam,
                    hurdle_sdgp=hurdle_sdgp_nam,
                    Intercept = intercept,
                    mu_area = mu_area_nam,
                    hurdle_mu_area=hurdle_mu_area_nam,
                    sigma_intercept = sigma_intercept,
                    hurdle_sigma_intercept=hurdle_sigma_intercept,
                    beta = beta_nam,
                    hurdle_beta = hurdle_beta_nam
)

new_data_sam = list(N = nrow(newlonlat_sam),
                    lonlat = newlonlat_sam,
                    max_lon = max(newlonlat_sam[,1]),
                    max_lat = max(newlonlat_sam[,2]),
                    M = c(m1,m2),
                    M_nD = M_nD,
                    indices = indices,
                    c1 = 4,
                    c2 = 4,
                    lscale = lscale_sam,
                    sdgp = sdgp_sam,
                    shape = shape,
                    hurdle_intercept = hurdle_intercept,
                    hurdle_lscale=hurdle_lscale_sam,
                    hurdle_sdgp=hurdle_sdgp_sam,
                    Intercept = intercept,
                    mu_area = mu_area_sam,
                    hurdle_mu_area=hurdle_mu_area_sam,
                    sigma_intercept = sigma_intercept,
                    hurdle_sigma_intercept=hurdle_sigma_intercept,
                    beta = beta_sam,
                    hurdle_beta = hurdle_beta_sam
)

new_data_pap = list(N = nrow(newlonlat_pap),
                    lonlat = newlonlat_pap,
                    max_lon = max(newlonlat_pap[,1]),
                    max_lat = max(newlonlat_pap[,2]),
                    M = c(m1,m2),
                    M_nD = M_nD,
                    indices = indices,
                    c1 = 4,
                    c2 = 4,
                    lscale = lscale_pap,
                    sdgp = sdgp_pap,
                    shape = shape,
                    hurdle_intercept = hurdle_intercept,
                    hurdle_lscale=hurdle_lscale_pap,
                    hurdle_sdgp=hurdle_sdgp_pap,
                    Intercept = intercept,
                    mu_area = mu_area_pap,
                    hurdle_mu_area=hurdle_mu_area_pap,
                    sigma_intercept = sigma_intercept,
                    hurdle_sigma_intercept=hurdle_sigma_intercept,
                    beta = beta_pap,
                    hurdle_beta = hurdle_beta_pap
)

new_data_aus = list(N = nrow(newlonlat_aus),
                    lonlat = newlonlat_aus,
                    max_lon = max(newlonlat_aus[,1]),
                    max_lat = max(newlonlat_aus[,2]),
                    M = c(m1,m2),
                    M_nD = M_nD,
                    indices = indices,
                    c1 = 4,
                    c2 = 4,
                    lscale = lscale_aus,
                    sdgp = sdgp_aus,
                    shape = shape,
                    hurdle_intercept = hurdle_intercept,
                    hurdle_lscale=hurdle_lscale_aus,
                    hurdle_sdgp=hurdle_sdgp_aus,
                    Intercept = intercept,
                    mu_area = mu_area_aus,
                    hurdle_mu_area=hurdle_mu_area_aus,
                    sigma_intercept = sigma_intercept,
                    hurdle_sigma_intercept=hurdle_sigma_intercept,
                    beta = beta_aus,
                    hurdle_beta = hurdle_beta_aus
)


gm_hurdle_gp = rstan::stan_model("./stan/approx-gp-preds.stan")

gp_preds_afr <- sampling(gm_hurdle_gp,
                           data   = new_data_afr, 
                           chains = 1,
                           cores  = 1,
                           warmup = 0,
                           iter   = 500,
                           seed   = 42,
                           pars = c("preds", "hurdle_preds"),
                           algorithm="Fixed_param")

gp_preds_eur <- sampling(gm_hurdle_gp, 
                         data   = new_data_eur, 
                         chains = 1,
                         cores  = 1,
                         warmup = 0,
                         iter   = 500,
                         seed   = 42,
                         pars = c("preds", "hurdle_preds"),
                         algorithm="Fixed_param")

gp_preds_nam <- sampling(gm_hurdle_gp, 
                         data   = new_data_nam, 
                         chains = 1,
                         cores  = 1,
                         warmup = 0,
                         iter   = 500,
                         seed   = 42,
                         pars = c("preds", "hurdle_preds"),
                         algorithm="Fixed_param")

gp_preds_sam <- sampling(gm_hurdle_gp, 
                         data   = new_data_sam, 
                         chains = 1,
                         cores  = 1,
                         warmup = 0,
                         iter   = 500,
                         seed   = 42,
                         pars = c("preds", "hurdle_preds"),
                         algorithm="Fixed_param")

gp_preds_aus <- sampling(gm_hurdle_gp, 
                         data   = new_data_aus, 
                         chains = 1,
                         cores  = 1,
                         warmup = 0,
                         iter   = 500,
                         seed   = 42,
                         pars = c("preds", "hurdle_preds"),
                         algorithm="Fixed_param")

gp_preds_pap <- sampling(gm_hurdle_gp, 
                         data   = new_data_pap, 
                         chains = 1,
                         cores  = 1,
                         warmup = 0,
                         iter   = 500,
                         seed   = 42,
                         pars = c("preds", "hurdle_preds"),
                         algorithm="Fixed_param")

## retrieve summarised predictions ##

afr_preds_df <- summarise_draws(subset_draws(as_draws(gp_preds_afr), variable = c("preds")))
afr_hu_df <- summarise_draws(subset_draws(as_draws(gp_preds_afr), variable = c("hurdle_preds")))
eur_preds_df <- summarise_draws(subset_draws(as_draws(gp_preds_eur), variable = c("preds")))
eur_hu_df <- summarise_draws(subset_draws(as_draws(gp_preds_eur), variable = c("hurdle_preds")))
nam_preds_df <- summarise_draws(subset_draws(as_draws(gp_preds_nam), variable = c("preds")))
nam_hu_df <- summarise_draws(subset_draws(as_draws(gp_preds_nam), variable = c("hurdle_preds")))
sam_preds_df <- summarise_draws(subset_draws(as_draws(gp_preds_sam), variable = c("preds")))
sam_hu_df <- summarise_draws(subset_draws(as_draws(gp_preds_sam), variable = c("hurdle_preds")))
pap_preds_df <- summarise_draws(subset_draws(as_draws(gp_preds_pap), variable = c("preds")))
pap_hu_df <- summarise_draws(subset_draws(as_draws(gp_preds_pap), variable = c("hurdle_preds")))
aus_preds_df <- summarise_draws(subset_draws(as_draws(gp_preds_aus), variable = c("preds")))
aus_hu_df <- summarise_draws(subset_draws(as_draws(gp_preds_aus), variable = c("hurdle_preds")))

hu_afr_mean <- afr_hu_df$mean
hu_eur_mean <- eur_hu_df$mean
hu_nam_mean <- nam_hu_df$mean
hu_sam_mean <- sam_hu_df$mean
hu_pap_mean <- pap_hu_df$mean
hu_aus_mean <- aus_hu_df$mean

counts_afr_mean <- round(afr_preds_df$mean)
counts_eur_mean <- round(eur_preds_df$mean)
counts_nam_mean <- round(nam_preds_df$mean)
counts_sam_mean <- round(sam_preds_df$mean)
counts_pap_mean <- round(pap_preds_df$mean)
counts_aus_mean <- round(aus_preds_df$mean)

## plot the predictions on a map ##

# fixed scale using the max mean predicted counts for all areas (found in Eurasia)
grandmin <- 0
grandmax <- max(counts_eur_mean)

humin <- 0
humax <- 1

# define the number of breaks
mybreaks <- seq(0, 34, length.out = 11)
hubreaks <- seq(0, 1, length.out = 11)

# add predictions to grid coordinates
newcoords_afr$hu_mean <- hu_afr_mean
newcoords_afr$count_mean <- counts_afr_mean
newcoords_eur$hu_mean <- hu_eur_mean
newcoords_eur$count_mean <- counts_eur_mean
newcoords_nam$hu_mean <- hu_nam_mean
newcoords_nam$count_mean <- counts_nam_mean
newcoords_sam$hu_mean <- hu_sam_mean
newcoords_sam$count_mean <- counts_sam_mean
newcoords_pap$hu_mean <- hu_pap_mean
newcoords_pap$count_mean <- counts_pap_mean
newcoords_aus$hu_mean <- hu_aus_mean
newcoords_aus$count_mean <- counts_aus_mean

myviridis = c("#fde725", "#b5de2b", "#6ece58", "#35b779", "#1f9e89", 
              "#26828e", "#31688e", "#3e4989", "#482878", "#440154")

## Africa ##

# if this throws an error, do all_langs <- st_drop_geometry(all_langs) to turn it into a non-sf dataframe
afr_langs <- filter(all_langs, Macroarea == "Africa") %>%
  mutate(has_0 = as.numeric(counts_wals == 0))
eur_langs <- filter(all_langs, Macroarea == "Eurasia") %>%
  mutate(has_0 = as.numeric(counts_wals == 0))
nam_langs <- filter(all_langs, Macroarea == "North America") %>%
  mutate(has_0 = as.numeric(counts_wals == 0))
sam_langs <- filter(all_langs, Macroarea == "South America") %>%
  mutate(has_0 = as.numeric(counts_wals == 0))
pap_langs <- filter(all_langs, Macroarea == "Multinesia") %>%
  mutate(has_0 = as.numeric(counts_wals == 0))
aus_langs <- filter(all_langs, Macroarea == "Australia") %>%
  mutate(has_0 = as.numeric(counts_wals == 0))

plot_hurdle_afr = ggplot() +
  geom_contour_filled(data = newcoords_afr
                      , aes(x = longitude
                            , y = latitude, z = hu_mean),
                      breaks = hubreaks) +
  geom_sf(data = africa, aes(geometry = geometry), alpha = 0, color = "black") +
  scale_fill_manual(name = "Prob. non-missing", 
                    values = rev(myviridis), 
                    drop=FALSE) +
  geom_point(data = afr_langs, alpha = 0.7,
             aes(x = Longitude, y = Latitude, shape = as.factor(has_0))) +
  theme(legend.position = "bottom"
        , legend.box = "vertical"
        , legend.key = element_rect(fill = "transparent"
                                    , colour = "transparent")) +
  labs(shape = "is missing", fill = "GP") +
  ylab("latitude") +
  xlab("longitude")

plot_counts_afr = ggplot() +
  geom_contour_filled(data = newcoords_afr
                      , aes(x = longitude
                            , y = latitude, z = count_mean),
                      breaks = mybreaks) +
  geom_sf(data = africa, aes(geometry = geometry), alpha = 0, color = "black") +
  scale_fill_manual(name = "Predicted counts", 
                    values = rev(myviridis), 
                    drop=FALSE) +
  geom_point(data = afr_langs, alpha = 0.7,
             aes(x = Longitude, y = Latitude, shape = as.factor(has_0))) +
  theme(legend.position = "bottom"
        , legend.box = "vertical"
        , legend.key = element_rect(fill = "transparent"
                                    , colour = "transparent")) +
  labs(shape = "is missing") +
  ylab("latitude") +
  xlab("longitude")

## Eurasia - needs to be projected before plotting ##

proj = paste('PROJCS["ProjWiz_Custom_Patterson",
 GEOGCS["GCS_WGS_1984",
  DATUM["D_WGS_1984",
   SPHEROID["WGS_1984",6378137.0,298.257223563]],
  PRIMEM["Greenwich",0.0],
  UNIT["Degree",0.0174532925199433]],
 PROJECTION["Patterson"],
 PARAMETER["False_Easting",0.0],
 PARAMETER["False_Northing",0.0],
 PARAMETER["Central_Meridian",20],
 UNIT["Meter",1.0]]', sep=(","))

eurasia = st_transform(eurasia, crs=proj)
preds_eur = st_as_sf(newcoords_eur, coords=c('longitude', 'latitude'), crs = 4326)
hu_preds_eur = st_as_sf(newcoords_eur, coords=c('longitude', 'latitude'), crs = 4326)
preds_eur_proj = st_transform(preds_eur, crs = proj)
hu_preds_eur_proj = st_transform(hu_preds_eur, crs = proj)
preds_eur_proj$longitude = st_coordinates(preds_eur_proj)[,1]
preds_eur_proj$latitude = st_coordinates(preds_eur_proj)[,2]
hu_preds_eur_proj$longitude = st_coordinates(hu_preds_eur_proj)[,1]
hu_preds_eur_proj$latitude = st_coordinates(hu_preds_eur_proj)[,2]
eur_langs = st_as_sf(eur_langs, coords=c('Longitude', 'Latitude'), crs = 4326)
eur_langs = st_transform(eur_langs, crs=proj)
eur_langs$proj_long = st_coordinates(eur_langs)[,1]
eur_langs$proj_lat = st_coordinates(eur_langs)[,2]

preds_eur_proj = st_crop(preds_eur_proj, xmin = -5000000, ymin=min(preds_eur_proj$latitude),
                         xmax = max(preds_eur_proj$longitude), ymax=max(preds_eur_proj$latitude))

hu_preds_eur_proj = st_crop(hu_preds_eur_proj, xmin = -5000000, ymin=min(hu_preds_eur_proj$latitude),
                            xmax = max(preds_eur_proj$longitude), ymax=max(hu_preds_eur_proj$latitude))

preds_eur_df = preds_eur_proj %>% 
  select('latitude', 'longitude', 'hu_mean', 'count_mean') %>% 
  as_tibble()

plot_hurdle_eur = ggplot() +
  geom_contour_filled(data = newcoords_eur
                      , aes(x = longitude
                            , y = latitude, z = hu_mean),
                      breaks = hubreaks) +
  geom_sf(data = continents, aes(geometry = geometry), alpha = 0, color = "black") +
  coord_sf(xlim = c(-20, 
                    max(newcoords_eur$longitude)), 
           ylim=c(0, 
                  80),
           expand=FALSE) + 
  scale_fill_manual(name = "Prob. non-missing", 
                    values = rev(myviridis), 
                    drop=FALSE) +
  geom_point(data = eur_langs, alpha = 0.7,
             aes(x = Longitude, y = Latitude, shape = as.factor(has_0))) +
  theme(legend.position = "bottom"
        , legend.box = "vertical"
        , legend.key = element_rect(fill = "transparent"
                                    , colour = "transparent")) +
  labs(shape = "is missing") +
  ylab("latitude") +
  xlab("longitude")

plot_counts_eur = ggplot() +
  geom_contour_filled(data = newcoords_eur
                      , aes(x = longitude
                            , y = latitude, z = count_mean),
                      breaks = mybreaks) +
  geom_sf(data = continents, aes(geometry = geometry), alpha = 0, color = "black") +
  coord_sf(xlim = c(-20, 
                    max(newcoords_eur$longitude)), 
           ylim=c(0, 
                  80),
           expand=FALSE) + 
  scale_fill_manual(name = "Predicted counts", 
                    values = rev(myviridis), 
                    drop=FALSE) +
  geom_point(data = eur_langs, alpha = 0.7,
             aes(x = Longitude, y = Latitude, shape = as.factor(has_0))) +
  theme(legend.position = "bottom"
        , legend.box = "vertical"
        , legend.key = element_rect(fill = "transparent"
                                    , colour = "transparent")) +
  labs(shape = "is missing") +
  ylab("latitude") +
  xlab("longitude")

## North America ##

plot_hurdle_nam = ggplot() +
  geom_contour_filled(data = newcoords_nam
                      , aes(x = longitude
                            , y = latitude, z = hu_mean),
                      breaks = hubreaks) +
  geom_sf(data = north_america, aes(geometry = geometry), alpha = 0, color = "black") +
  coord_sf(xlim = c(min(newcoords_nam$longitude), 0), 
           ylim=c(min(newcoords_nam$latitude), max(newcoords_nam$latitude)), 
           expand=FALSE) +
  scale_fill_manual(name = "Prob. non-missing", 
                    values = rev(myviridis), 
                    drop=FALSE) +
  theme(legend.position = "bottom"
        , legend.box = "vertical"
        , legend.key = element_rect(fill = "transparent"
                                    , colour = "transparent")) +
  ylab("latitude") +
  xlab("longitude")

plot_counts_nam = ggplot() +
  geom_contour_filled(data = newcoords_nam
                      , aes(x = longitude
                            , y = latitude, z = count_mean),
                      breaks = mybreaks) +
  geom_sf(data = north_america, aes(geometry = geometry), alpha = 0, color = "black") +
  coord_sf(xlim = c(min(newcoords_nam$longitude), 0), 
           ylim=c(min(newcoords_nam$latitude), max(newcoords_nam$latitude)),
           expand=FALSE) + ## necessary to not include Hawaii
  scale_fill_manual(name = "Predicted counts", 
                    values = rev(myviridis), 
                    drop=FALSE) +
  geom_point(data = nam_langs, alpha = 0.7,
             aes(x = Longitude, y = Latitude, shape = as.factor(has_0))) +
  theme(legend.position = "bottom"
        , legend.box = "vertical"
        , legend.key = element_rect(fill = "transparent"
                                    , colour = "transparent")) +
  labs(shape = "is missing") + 
  ylab("latitude") +
  xlab("longitude")

## South America ##

plot_hurdle_sam = ggplot() +
  geom_contour_filled(data = newcoords_sam
                      , aes(x = longitude
                            , y = latitude, z = hu_mean),
                      breaks = hubreaks) +
  geom_sf(data = south_america, aes(geometry = geometry), alpha = 0, color = "black") +
  geom_point(data = sam_langs, alpha = 0.7,
             aes(x = Longitude, y = Latitude, shape = as.factor(has_0))) +
  scale_fill_manual(name = "Prob. non-missing", 
                    values = rev(myviridis), 
                    drop=FALSE) +
  theme(legend.position = "bottom"
        , legend.box = "vertical"
        , legend.key = element_rect(fill = "transparent"
                                    , colour = "transparent")) +
  labs(shape = "is missing") +
  ylab("latitude") +
  xlab("longitude")

plot_counts_sam = ggplot() +
  geom_contour_filled(data = newcoords_sam
                      , aes(x = longitude
                            , y = latitude, z = count_mean),
                      breaks = mybreaks) +
  geom_sf(data = south_america, aes(geometry = geometry), alpha = 0, color = "black") +
  scale_fill_manual(name = "Predicted counts", 
                    values = rev(myviridis), 
                    #values=breaklabel_hu(10), 
                    drop=FALSE) +
  geom_point(data = sam_langs, alpha = 0.7,
             aes(x = Longitude, y = Latitude, shape = as.factor(has_0))) +
  theme(legend.position = "bottom"
        , legend.box = "vertical"
        , legend.key = element_rect(fill = "transparent"
                                    , colour = "transparent")) +
  labs(shape = "is missing") + 
  ylab("latitude") +
  xlab("longitude")


## Australia ##

plot_hurdle_aus = ggplot() +
  geom_contour_filled(data = newcoords_aus
                      , aes(x = longitude
                            , y = latitude, z = hu_mean),
                      breaks = hubreaks) +
  geom_sf(data = australia, aes(geometry = geometry), alpha = 0, color = "black") +
  coord_sf(xlim = c(min(newcoords_aus$longitude), max(newcoords_aus$longitude)), 
           ylim=c(-45, -8), 
           expand=FALSE) +
  geom_point(data = aus_langs, alpha = 0.7,
             aes(x = Longitude, y = Latitude, shape = as.factor(has_0))) +
  scale_fill_manual(name = "Prob. non-missing", 
                    values = rev(myviridis), 
                    #values=breaklabel_hu(10), 
                    drop=FALSE) +
  theme(legend.position = "bottom"
        , legend.box = "vertical"
        , legend.key = element_rect(fill = "transparent"
                                    , colour = "transparent")) +
  labs(shape = "is missing") +
  ylab("latitude") +
  xlab("longitude")


plot_counts_aus = ggplot() +
  geom_contour_filled(data = newcoords_aus
                      , aes(x = longitude
                            , y = latitude, z = count_mean),
                      breaks = mybreaks) +
  geom_sf(data = australia, aes(geometry = geometry), alpha = 0, color = "black") +
  coord_sf(xlim = c(min(newcoords_aus$longitude), max(newcoords_aus$longitude)), 
           ylim=c(-45, -8), 
           expand=FALSE) +
  geom_point(data = aus_langs, alpha = 0.7,
             aes(x = Longitude, y = Latitude, shape = as.factor(has_0))) +
  scale_fill_manual(name = "Predicted counts", 
                    values = rev(myviridis), 
                    #values=breaklabel_hu(10), 
                    drop=FALSE) +
  theme(legend.position = "bottom"
        , legend.box = "vertical"
        , legend.key = element_rect(fill = "transparent"
                                    , colour = "transparent")) +
  labs(shape = "is missing") +
  ylab("latitude") +
  xlab("longitude")

## Multinesia ##

plot_hurdle_pap = ggplot() +
  geom_contour_filled(data = newcoords_pap
                      , aes(x = longitude
                            , y = latitude, z = hu_mean),
                      breaks = hubreaks) +
  geom_sf(data = continents, aes(geometry = geometry), alpha = 0, color = "black") +
  coord_sf(xlim = c(min(newcoords_pap$longitude), max(newcoords_pap$longitude)), 
           ylim=c(min(newcoords_pap$latitude), max(newcoords_pap$latitude)), 
           expand=FALSE) +
  geom_point(data = pap_langs, alpha = 0.7,
             aes(x = Longitude, y = Latitude, shape = as.factor(has_0))) +
  scale_fill_manual(name = "Prob. non-missing", 
                    values = rev(myviridis), 
                    #values=breaklabel_hu(10), 
                    drop=FALSE) +
  theme(legend.position = "bottom"
        , legend.box = "vertical"
        , legend.key = element_rect(fill = "transparent"
                                    , colour = "transparent")) +
  labs(shape = "is missing") +
  ylab("latitude") +
  xlab("longitude")

plot_counts_pap = ggplot() +
  geom_contour_filled(data = newcoords_pap
                      , aes(x = longitude
                            , y = latitude, z = count_mean),
                      breaks = mybreaks) +
  geom_sf(data = continents, aes(geometry = geometry), alpha = 0, color = "black") +
  coord_sf(xlim = c(min(newcoords_pap$longitude), max(newcoords_pap$longitude)), 
           ylim=c(min(newcoords_pap$latitude), max(newcoords_pap$latitude)), 
           expand=FALSE) +
  geom_point(data = pap_langs, alpha = 0.7,
             aes(x = Longitude, y = Latitude, shape = as.factor(has_0))) +
  scale_fill_manual(name = "Predicted counts", 
                    values = rev(myviridis), 
                    #values=breaklabel_hu(10), 
                    drop=FALSE) +
  theme(legend.position = "bottom"
        , legend.box = "vertical"
        , legend.key = element_rect(fill = "transparent"
                                    , colour = "transparent")) +
  labs(shape = "is missing") +
  ylab("latitude") +
  xlab("longitude")
