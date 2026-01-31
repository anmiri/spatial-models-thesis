library(tidyverse)
library(sf)

noLV <- read_tsv('./data/NoLV.tsv') # languages in this df have no LVs
LVfreq <- read_tsv('./data/LVFreq.tsv') 
LVnf <- read_tsv('./data/LVNoFreq.tsv') # languages in this df have LVs
LVswadesh <- read_tsv('./data/LVSwadesh200.tsv')
LVswadesh$LVFreq
LVswadesh$LVFreq200

LVfreq$LVFreq

plot(density(LVfreq$LVTypes))

## create a df with all the data ##
LVfreq$has_lv <- 1
LVnf$has_lv <- 1
LVnf$LVTypes <- -1 # absent
noLV$has_lv <- 0
noLV$LVTypes <- 0
noLV$LVFreq <- 0
LVnf$LVFreq <- -1
noLV <- noLV %>% select(code_glot, has_lv, LVTypes, LVFreq, lat, long, family_wals, genus)
LVfreq <- LVfreq %>% select(code_glot, has_lv, LVTypes, LVFreq, lat, long, family_wals, genus)
LVnf <- LVnf %>% select(code_glot, has_lv, LVTypes, LVFreq, lat, long, family_wals, genus)

d <- rbind(LVfreq, LVnf, noLV) %>%
  arrange(code_glot)

write_csv(d, './data/LV_data.csv')

d$LVFreq <- str_replace(d$LVFreq, ",", ".")
d$LVFreq <- as.numeric(d$LVFreq)

median(d$LVFreq[d$LVFreq != 0 & d$LVFreq != -1])
d$LVFreq[d$LVFreq == -1] <- NA
d$LVpct <- d$LVFreq/100
plot(density(d$LVpct))
standardize(d$LVFreq)

dfreq <- d %>% drop_na(LVFreq)
dsf <- dfreq %>% sf::st_as_sf(coords = c('long', 'lat'), crs = 4326)
write_csv(dfreq, './data/LV_freq_data.csv')
