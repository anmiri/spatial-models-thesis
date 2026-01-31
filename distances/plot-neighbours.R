library(tidyverse)
library(geostan)
library(spdep)
library(sf)

d <- read_csv('./data/LV_freq_data.csv')

africa_sf <- spData::world %>% filter(continent == "Africa")

## dist weights: IDW ##
d <- read_csv('./data/LV_freq_data.csv')
dsf <- st_as_sf(d, coords = c("long", "lat"), crs = 4326)

# get the spatial extent of the languages with frequency data #
plot(density(d$LVpct))
plot(density(d$LVFreq))

ids <- which(dsf$LVpct != 0)
st_bbox(dsf[ids,])

dc <- st_crop(dsf, st_bbox(dsf[ids,]))
plot(density(dc$LVFreq))
dc$LVnorm <- scale(dc$LVFreq)
plot(density(dc$LVnorm))

## get the tree ##

library(ape)

df_tree <- read_rds('./data/phylogeny.rds')
phylo <- ape::vcv.phylo(df_tree, corr=TRUE)
dphy <- dc %>% filter(code_glot %in% df_tree$tip.label) %>%
  arrange(code_glot)
dphy <- dphy[!duplicated(dphy$code_glot),]
phymat <- phylo[dphy$code_glot, dphy$code_glot]
all(dphy$code_glot == rownames(phymat))

### neighbour matrix  w/ distance weights ###

dnear <- dnearneigh(dphy, 0, 500) # symmetric

listidw <- nb2listwdist(dnear, dphy, type="idw", style="raw", 
                        alpha = 1, zero.policy=TRUE)

listexp <- nb2listwdist(dnear, dphy, type="exp", style="raw", 
                        alpha = 1, zero.policy=TRUE)


idw_W <- listw2mat(listidw)
idw_Ws <- row_standardize(idw_W)
exp_W <- listw2mat(listexp)
exp_Ws <- row_standardize(exp_W)
sar_exp <- prep_sar_data(exp_Ws)
sar_idw <- prep_sar_data(idw_Ws)

W <- sar_idw$W
rownames(W) <- dphy$code_glot
colnames(W) <- dphy$code_glot
W <- as.data.frame(W)
W$code_glot <- dphy$code_glot

dphy$longitude <- st_coordinates(dphy$geometry)[,1]
dphy$latitude <- st_coordinates(dphy$geometry)[,2]

wdf <- W %>% pivot_longer(-code_glot) %>%
  mutate(group = str_c(pmin(code_glot, name), pmax(code_glot, name))) %>%
  left_join(dphy, by = "code_glot") %>%
  distinct()

expp <- ggplot() +
  geom_sf(data = africa_sf, fill = 'antiquewhite1') +
  geom_line(data = wdf
            , aes(x = longitude, y = latitude, group = group, alpha = value)
            , color = "darkblue", linewidth = 0.8) +
  geom_point(data = dphy, aes(x = longitude, y = latitude)
          , size = 1
          , color = "purple4") +
  scale_alpha_continuous(name = "exponential distance weights, r = 500"
                         , range = c(0, 1)) +
  ggtitle("Spatial weights") +
  theme(plot.title = element_text(size = 20, face = "bold"), 
        panel.background = element_rect(fill = "aliceblue"),
        legend.position='bottom') +
  coord_sf(xlim = c(-17, 32),
           ylim = c(-6, 13))

ggsave('./plots/expdist-map.pdf', expp, height = 6, width = 11)


idwp <- ggplot() +
  geom_sf(data = africa_sf, fill = 'antiquewhite1') +
  geom_line(data = wdf
            , aes(x = longitude, y = latitude, group = group, alpha = value)
            , color = "darkblue", linewidth = 0.8) +
  geom_point(data = dphy, aes(x = longitude, y = latitude)
             , size = 1
             , color = "purple4") +
  scale_alpha_continuous(name = "inverse distance weights, r = 500"
                         , range = c(0, 1)) +
  ggtitle("Spatial weights") +
  theme(plot.title = element_text(size = 20, face = "bold"), 
        panel.background = element_rect(fill = "aliceblue"),
        legend.position='bottom') +
  coord_sf(xlim = c(-17, 32),
           ylim = c(-6, 13))

ggsave('./plots/idwdist-map.pdf', idwp, height = 6, width = 11)

