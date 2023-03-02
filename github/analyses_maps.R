# Objective: make maps showing sites within ecoregions for abundance and diversity
# datasets

# load package
library(tidyverse)
library(sf)
library(mapproj)

# read in datasets
# drop ecoregion-year combos for which SC is not between 0.7 and 0.8 from region
# level data and site level data

# diversity
dat_div_all <- read.csv("./data/dat_div_NRSA_USGS_4Jan.csv", header = T, 
                        row.names = 1) %>%
  filter(sc >= 0.7 & sc <= 0.8) %>%
  unite("EcoRegion_Year", c(Ecoregion_NARS, Year), remove = F)

# site from which region level diversity metrics are calculated
dat_alpha_all <- read.csv("./data/dat_alpha_all_4Jan.csv", header = T,
                          row.names = 1) %>%
  unite("EcoRegion_Year", c(Ecoregion_NARS, Year), remove = F) %>%
  #exclude ER-YR combos that don't meet the SC range
  filter(EcoRegion_Year %in% dat_div_all$EcoRegion_Year) 

# abundance
dat_abun_all <- read.csv("./data/dat_abun_all.csv", header = T,
                          row.names = 1) 

###
### Make a make showing sites for abundance and diversity from EPA and USGS
###

# diversity sites

# collapse dataset into sites only
dat_site_div<- dat_alpha_all %>%
  group_by(SiteNumber) %>%
  slice(1)

# read in shapefiles for ecoregions
ecor <- st_read("./ecoregion_shapefile/Aggr_Ecoregions_2015.shp") 

# change projection to NAD84 to match sites
ecor <- st_transform(ecor, crs = 4269)
# reduce size of shapefile to make plot file smaller
ecor <- rmapshaper::ms_simplify(ecor, keep = 0.01, keep_shapes = T)

# create a lab for panel plot
divlab <- c(EPA = "EPA - diversity",
            USGS = "USGS - diversity")

# two panel plot showing diversity sites from EPA and USGS
f1 <- ggplot() +
  geom_sf(data = ecor,
          aes(fill = WSA9),
          alpha = 0.4) +  
  geom_point(data = dat_site_div,
             aes(x = Longitude_dd,
                 y = Latitude_dd),
             shape = 21,
             alpha = 0.4) +
  facet_wrap(~Agency,
             labeller = labeller(Agency = divlab)) +
  xlab("") +
  ylab("") +
  scale_fill_brewer(palette="Set1",
                    name = "Ecoregions",
                    labels = c("Coastal Plains (CPL)",
                               "Northern Appalachians (NAP)",
                               "Northern Plains (NPL)",
                               "Southern Appalachians (SAP)",
                               "Southern Plains (SPL)",
                               "Temperate Plains (TPL)",
                               "Upper Midwest (UMW)",
                               "Western Mountains (WMT)",
                               "Xeric (XER)")) +
  theme_bw() +
  theme(legend.position="top") 

# abundance sites

# collapse dataset into sites only
dat_site_abun <- dat_abun_all %>%
  group_by(SiteNumber) %>%
  slice(1)

# create a lab for panel plot
denlab <- c(EPA = "EPA - density",
            USGS = "USGS - density")

# two panel plot showing diversity sites from EPA and USGS
f2 <- ggplot() +
  geom_sf(data = ecor,
          aes(fill = WSA9),
          alpha = 0.4) +  
  geom_point(data = dat_site_abun,
             aes(x = Longitude_dd,
                 y = Latitude_dd),
             shape = 21,
             alpha = 0.4) +
  facet_wrap(~Agency,
             labeller = labeller(Agency = denlab)) +
  xlab("") +
  ylab("") +
  scale_fill_brewer(palette="Set1",
                    name = "Ecoregions",
                    labels = c("Coastal Plains (CPL)",
                               "Northern Appalachians (NAP)",
                               "Northern Plains (NPL)",
                               "Southern Appalachians (SAP)",
                               "Southern Plains (SPL)",
                               "Temperate Plains (TPL)",
                               "Upper Midwest (UMW)",
                               "Western Mountains (WMT)",
                               "Xeric (XER)")) +
  theme_bw() +
  theme(legend.position="top") 

# join panels together
A = f2 + theme(legend.position = "none")
B = f1 + theme(legend.position = "none")
legend = cowplot::get_legend(f1)

f2 <- cowplot::plot_grid(A, B, labels = c("A", "B"),
                         label_size = 20,
                         ncol = 1,
                         align = "v")
f2_legend <- cowplot::plot_grid(legend, f2, ncol = 1, rel_heights = c(2,15))

ggsave("./figures_supplement/map_abun_div.png", f2_legend, dpi = 200, height = 8, 
       width = 10)

### Make one map with sites

dat_site_all <- unique(rbind(dat_site_abun %>% select(Longitude_dd,
                                         Latitude_dd,
                                         SiteNumber,
                                         Agency),
                      dat_site_div %>% select(Longitude_dd,
                                               Latitude_dd,
                                               SiteNumber,
                                               Agency)
                      
))


fig_main <- ggplot() +
  geom_sf(data = ecor,
          aes(fill = WSA9),
          alpha = 0.4) +  
  geom_point(data = dat_site_all,
             aes(x = Longitude_dd,
                 y = Latitude_dd),
             shape = 21,
             alpha = 0.4) +
  facet_wrap(~Agency) +
  xlab("") +
  ylab("") +
  scale_fill_brewer(palette="Set1",
                    name = "Ecoregions",
                    labels = c("Coastal Plains (CPL)",
                               "Northern Appalachians (NAP)",
                               "Northern Plains (NPL)",
                               "Southern Appalachians (SAP)",
                               "Southern Plains (SPL)",
                               "Temperate Plains (TPL)",
                               "Upper Midwest (UMW)",
                               "Western Mountains (WMT)",
                               "Xeric (XER)")) +
  theme_bw() +
  theme(legend.position="top") 

fig_main

ggsave("./figures_maintext/map_all.png", fig_main, dpi = 200, height = 8, 
       width = 10)
