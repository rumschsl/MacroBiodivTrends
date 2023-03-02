# 30 September 2021
# Samantha Rumschlag
# Objective: Examine correlations that exist among total density, area sampled,
# proportion sample ID'd, and time

# Update: added panels for USGS and EPA data

# Checked by SLR on Feb 23, 2022

# load packages
library(tidyverse)

# read in data
dat_abun_all <- read.csv("./data/dat_abun_all.csv", header = T,
                         row.names = 1) 

# AREA SAMPLED VS TIME
area_time_fig <- dat_abun_all %>%
  ggplot(aes(x = YearCont + 1992, y = AreaSampTot_m2, color = Ecoregion_NARS,
             fill = Ecoregion_NARS)) +
  geom_point(alpha = 0.3) +
  facet_wrap(Ecoregion_NARS~Agency, ncol=6) +
  scale_x_continuous(breaks=c(1993,1998,2003,2008,2013,2018))+
  scale_color_brewer(palette="Set1", name = "Ecoregion") +
  scale_fill_brewer(palette="Set1", name = "Ecoregion") +
  geom_smooth(method = "lm") +
  ylab(expression("Total area sampled (m"^2*")")) +
  xlab("Year") +
  theme_bw() +
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1))

area_time_fig

# AREA SAMPLED VS DENSITY
area_dens_fig <- dat_abun_all %>%
  ggplot(aes(x = AreaSampTot_m2, y = tot_abun, color = Ecoregion_NARS,
             fill = Ecoregion_NARS)) +
  geom_point(alpha = 0.3) +
  facet_wrap(Ecoregion_NARS~Agency, ncol=6) +
  scale_color_brewer(palette="Set1", name = "Ecoregion") +
  scale_fill_brewer(palette="Set1", name = "Ecoregion") +
  geom_smooth(method = "lm") +
  ylab(expression("Total density (ind. / m"^2*")")) +
  xlab(expression("Total area sampled (m"^2*")")) +
  scale_y_continuous(trans = "log", breaks =c(10, 150, 3000, 60000) )+
  theme_bw() +
  theme(text = element_text(size = 18))

area_dens_fig

# join panels together
A = area_time_fig + theme(legend.position = "top") +
  guides(colour = guide_legend(nrow = 1))
B = area_dens_fig + theme(legend.position = "none")

f2 <- cowplot::plot_grid(A, B, labels = c("A", "B"),
                         label_size = 20,
                         ncol = 1,
                         align = "v")
ggsave("./figures_supplement/area_time_dens.png", f2, dpi = 200, height = 16, width = 12)


# Proportion total sample identified VS TIME

propID_time_fig <- dat_abun_all %>%
  ggplot(aes(x = YearCont + 1992, y = PropID, 
             color = Ecoregion_NARS,
             fill = Ecoregion_NARS)) +
  geom_point(alpha = 0.3) +
  facet_wrap(Ecoregion_NARS~Agency, ncol=6) +
  scale_x_continuous(breaks=c(1993,1998,2003,2008,2013,2018))+
  scale_color_brewer(palette="Set1", name = "Ecoregion") +
  scale_fill_brewer(palette="Set1", name = "Ecoregion") +
  geom_smooth(method = "lm") +
  ylab("Proportion total sample identified") +
  xlab("Year") +
  theme_bw() +
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1))

propID_time_fig

# Proportion total sample identified

propID_dens_fig <- dat_abun_all %>%
  ggplot(aes(x = PropID, 
             y = tot_abun, color = Ecoregion_NARS,
             fill = Ecoregion_NARS)) +
  geom_point(alpha = 0.3) +
  facet_wrap(Ecoregion_NARS~Agency, ncol=6) +
  scale_color_brewer(palette="Set1", name = "Ecoregion") +
  scale_fill_brewer(palette="Set1", name = "Ecoregion") +
  geom_smooth(method = "lm") +
  ylab(expression("Total density (ind. / m"^2*")")) +
  xlab("Proportion total sample identified") +
  scale_y_continuous(trans = "log", breaks =c(10, 150, 3000, 60000) )+
  scale_x_continuous(breaks =c(0, 0.25, 0.5, 0.75, 1.00), limits = c(0,1) )+
  theme_bw() +
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1))

propID_dens_fig

# join panels together
A = propID_time_fig + theme(legend.position = "top") +
  guides(colour = guide_legend(nrow = 1))
B = propID_dens_fig + theme(legend.position = "none")

f2 <- cowplot::plot_grid(A, B, labels = c("A", "B"),
                         label_size = 20,
                         ncol = 1,
                         align = "v")

ggsave("./figures_supplement/propID_time_dens.png", f2, dpi = 200, height = 16, width = 12)

