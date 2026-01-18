library(rstan)
library(tidyverse)
library(loo)

all_langs = read_csv("./data/df-wals.csv")

n_areas <- length(unique(all_langs$Macroarea))

ids_afr <- which(all_langs$Macroarea == "Africa")
ids_eur <- which(all_langs$Macroarea == "Eurasia")
ids_nam <- which(all_langs$Macroarea == "North America")
ids_sam <- which(all_langs$Macroarea == "South America")
ids_pap <- which(all_langs$Macroarea == "Multinesia")
ids_aus <- which(all_langs$Macroarea == "Australia")

nAfr <- length(ids_afr)
nEur <- length(ids_eur)
nNam <- length(ids_nam)
nSam <- length(ids_sam)
nPap <- length(ids_pap)
nAus <- length(ids_aus)

nGr <- c(nAfr,nEur,nNam,nSam,nPap,nAus)

coords <- cbind(all_langs$Longitude, all_langs$Latitude)

loc_afr <- coords[ids_afr,]
loc_eur <- coords[ids_eur,]
loc_nam <- coords[ids_nam,]
loc_sam <- coords[ids_sam,]
loc_pap <- coords[ids_pap,]
loc_aus <- coords[ids_aus,]

# scale long and lat in order to estimate lscale on a not-ridiculous scale
loc_afr[,1] <- scale(loc_afr[,1])[,1]
loc_eur[,1] <- scale(loc_eur[,1])[,1]
loc_nam[,1] <- scale(loc_nam[,1])[,1]
loc_sam[,1] <- scale(loc_sam[,1])[,1]
loc_pap[,1] <- scale(loc_pap[,1])[,1]
loc_aus[,1] <- scale(loc_aus[,1])[,1]

loc_afr[,2] <- scale(loc_afr[,2])[,1]
loc_eur[,2] <- scale(loc_eur[,2])[,1]
loc_nam[,2] <- scale(loc_nam[,2])[,1]
loc_sam[,2] <- scale(loc_sam[,2])[,1]
loc_pap[,2] <- scale(loc_pap[,2])[,1]
loc_aus[,2] <- scale(loc_aus[,2])[,1]

## m1 and m2 are the number of basis functions
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

M_nD = m1*m2

## model with only an approximate Gaussian Process ##
hurdle_gp <- rstan::stan_model("./stan/hurdle/hurdle-negbin-gp.stan")

bfdata_hurdle <- list(N = nrow(all_langs)
                    , nAreas = n_areas
                    , y = all_langs$counts_wals
                    , ngr = nGr 
                    , ids_afr = ids_afr 
                    , ids_eur = ids_eur 
                    , ids_nam = ids_nam
                    , ids_sam = ids_sam
                    , ids_pap = ids_pap
                    , ids_aus = ids_aus 
                    , lonlat_afr = loc_afr
                    , lonlat_eur = loc_eur
                    , lonlat_nam = loc_nam
                    , lonlat_sam = loc_sam 
                    , lonlat_pap = loc_pap
                    , lonlat_aus = loc_aus
                    , max_lon_afr = max(loc_afr[,1])
                    , max_lat_afr = max(loc_afr[,2])
                    , max_lon_eur = max(loc_eur[,1])
                    , max_lat_eur = max(loc_eur[,2])
                    , max_lon_nam = max(loc_nam[,1])
                    , max_lat_nam = max(loc_nam[,2])
                    , max_lon_sam = max(loc_sam[,1])
                    , max_lat_sam = max(loc_sam[,2])
                    , max_lon_pap = max(loc_pap[,1])
                    , max_lat_pap = max(loc_pap[,2])
                    , max_lon_aus = max(loc_aus[,1])
                    , max_lat_aus = max(loc_aus[,2])
                    , c1 = 4 # boundary condition 
                    , c2 = 4
                    , M= c(m1, m2)
                    , M_nD= m1*m2
                    , indices= indices
                      )

fit_gp <- sampling(hurdle_gp,
                    data   = bfdata_hurdle, 
                    chains = 4,
                    cores  = 4,
                    warmup = 2000,
                    iter   = 4000,
                    seed   = 42,
                    control = list(adapt_delta = 0.9))

## model with approx. GP and phylogenetic regression ##
hurdle_gp_phy = rstan::stan_model("./stan/hurdle/hurdle-negbin-gp-phylo.stan")

library(ape)

## retrieve the Glottolog phylogeny ##

phylo <- read_rds("./data/phylo-mis.rds")
dphylo <- ape::vcv.phylo(phylo, corr = TRUE)

## this is the covariance matrix we use in the model ##

chol_phylo <- t(chol(dphylo[all_langs$Glottocode, 
                            all_langs$Glottocode]))

mdata_phylo <- list(N = nrow(all_langs)
                      , nAreas = n_areas
                      , y = all_langs$counts_wals
                      , ngr = nGr 
                      , Lcov = chol_phylo
                      , ids_afr = ids_afr 
                      , ids_eur = ids_eur 
                      , ids_nam = ids_nam
                      , ids_sam = ids_sam
                      , ids_pap = ids_pap
                      , ids_aus = ids_aus 
                      , lonlat_afr = loc_afr
                      , lonlat_eur = loc_eur
                      , lonlat_nam = loc_nam
                      , lonlat_sam = loc_sam 
                      , lonlat_pap = loc_pap
                      , lonlat_aus = loc_aus
                      , max_lon_afr = max(loc_afr[,1])
                      , max_lat_afr = max(loc_afr[,2])
                      , max_lon_eur = max(loc_eur[,1])
                      , max_lat_eur = max(loc_eur[,2])
                      , max_lon_nam = max(loc_nam[,1])
                      , max_lat_nam = max(loc_nam[,2])
                      , max_lon_sam = max(loc_sam[,1])
                      , max_lat_sam = max(loc_sam[,2])
                      , max_lon_pap = max(loc_pap[,1])
                      , max_lat_pap = max(loc_pap[,2])
                      , max_lon_aus = max(loc_aus[,1])
                      , max_lat_aus = max(loc_aus[,2])
                      , c1 = 4  
                      , c2 = 4
                      , M= c(m1, m2)
                      , M_nD= m1*m2
                      , indices= indices
)

fit_gp_phylo <- sampling(hurdle_gp_phy,
                            data   = mdata_phylo, 
                            chains = 4,
                            cores  = 4,
                            warmup = 2000,
                            iter   = 4000,
                            seed   = 42,
                            control = list(adapt_delta = 0.9))

## model with only phylogenetic regression, no GP ##
hurdle_phy = rstan::stan_model("./stan/hurdle-negbin-phylo.stan")

mdata_phy <- list(N = nrow(all_langs)
                    , nAreas = n_areas
                    , y = all_langs$counts_wals
                    , ngr = nGr 
                    , ids_afr = ids_afr 
                    , ids_eur = ids_eur 
                    , ids_nam = ids_nam
                    , ids_sam = ids_sam
                    , ids_pap = ids_pap
                    , ids_aus = ids_aus 
                    , Lcov = chol_phylo
)

fit_phylo <- sampling(hurdle_phy,
                         data   = mdata_phy, 
                         chains = 4,
                         cores  = 4,
                         warmup = 2000,
                         iter   = 4000,
                         seed   = 42,
                         control = list(adapt_delta = 0.9))

## hurdle negative binomial model without any controls ##
hurdle_base = rstan::stan_model("./stan/hurdle-negbin.stan")

mdata_base <- list(N = nrow(all_langs)
                  , nAreas = n_areas
                  , y = all_langs$counts_wals
                  , ngr = nGr 
                  , ids_afr = ids_afr 
                  , ids_eur = ids_eur 
                  , ids_nam = ids_nam
                  , ids_sam = ids_sam
                  , ids_pap = ids_pap
                  , ids_aus = ids_aus 
)

fit_base <- sampling(hurdle_base,
                      data   = mdata_base, 
                      chains = 4,
                      cores  = 4,
                      warmup = 2000,
                      iter   = 4000,
                      seed   = 42,
                      control = list(adapt_delta = 0.9))


## loo CV ##
library(loo)
loo_phylo <- loo(fit_phylo, save_psis = TRUE, cores = 10)
loo_gp <- loo(fit_gp, save_psis = TRUE, cores = 10)
loo_base <- loo(fit_base, save_psis = TRUE, cores = 10)
loo_gp_phylo <- loo(fit_gp_phylo, save_psis = TRUE, cores = 10)

loo_compare(loo_gp_phylo, 
            loo_gp,
            loo_phylo, 
            loo_base)
