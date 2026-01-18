-- Environment

* requirements.txt: you can use this create a Conda environment with all the required packages for running the models.

Alternatively, you can manually install the following R packages. We used R version 4.3.2:
 - tidyverse 
 - rstan
 - loo
 - posterior 
 - sf
 - ggplot2

-- R code:

* generate-preds.R : This file contains the code for generating posterior predictions and spatial predictions for new locations along a grid.
It also contains the code for doing the posterior predictive checks shown in the Supplementary Materials, and for plotting the spatial predictions on a map of each macroarea.
The file assumes you have the results of the model in a folder called `results`.

* run.R : Code for running the Stan models including LOO-CV 

-- Stan code:

* stan/hurdle-negbin.stan: hurdle negative binomial model without any controls
* stan/hurdle-negbin-gp.stan: hurdle negative binomial model with a latent approximate GP as an areal control
* stan/hurdle-negbin-phylo.stan: hurdle negative binomial model with phylogenetic controls
* stan/hurdle-negbin-gp-phylo.stan: multivaraite probit model with both areal and phylogenetic controls
* stan/approx-gp-preds.stan: script for drawing spatial predictions for new locations given model parameters
* stan/approx-gp-preds.stan: script for drawing posterior predictions given model parameters

-- Data:

* phylo-mis.rds : phylogeny based on Glottolog 4.8

* df-wals: count data for WALS and Glottolog languages.
  -- Columns:
  - ID: language code from WALS
  - id2: Glottocode
  - Longitude: longitude information
  - Latitude: latitude information
  - Macroarea: macroarea
  - counts_wals: counts for languages in WALS, 0 if language is not in WALS

* World_Continents.shp: The shapefiles (continents without national borders) used for plotting the maps.