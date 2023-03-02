# Objective: Diversity analyses according to landuse and including random effect 
# of agency with linear models with insects only

# load packages
library(tidyverse)
library(lme4)
library(ggeffects)
library(emmeans)
library(glmmTMB)

# read in density data
dat_abun <- read.csv("./insects_only/data/dat_abun_all.csv", row.names = 1, header =T) 

# read in diversity data
dat_div_all_LU_Ag <- read.csv("./insects_only/data/dat_div_all_LU_Ag.csv", 
                              header = T, 
                              row.names = 1) %>%
  # drop ER_YR_LU_AG combos not estimated within 0.7-0.8 sample coverage
  filter(sc >= 0.7 & sc <= 0.8)

###         ###
### DENSITY ###
###         ###

# this model fits
m2 <- lmer(log(tot_abun+1) ~ YearCont + landuse  + Ecoregion_NARS +
             PropID + AreaSampTot_m2 + Agency + Gen_ID_Prop +
             Ecoregion_NARS:Agency +
             YearCont:landuse + 
             YearCont:Ecoregion_NARS +
             (1|SiteNumber),
           data = dat_abun,
           REML = TRUE)

performance::r2(m2)
# R2 for Mixed Models
# 
# Conditional R2: 0.793
#    Marginal R2: 0.737

summary(m2)

car::Anova(m2, type = "III", test = "F") 
# Analysis of Deviance Table (Type III Wald F tests with Kenward-Roger df)
# 
# Response: log(tot_abun + 1)
#                                  F Df Df.res    Pr(>F)    
# (Intercept)              3741.4417  1 7668.5 < 2.2e-16 ***
# YearCont                    0.1176  1 7741.5 0.7317015    
# landuse                     2.2015  3 6699.5 0.0857386 .  
# Ecoregion_NARS              4.7806  8 7376.8 7.029e-06 ***
# PropID                  14236.5621  1 8002.2 < 2.2e-16 ***
# AreaSampTot_m2           1597.7208  1 7488.2 < 2.2e-16 ***
# Agency                     12.4548  1 5032.0 0.0004207 ***
# Gen_ID_Prop              1356.6655  1 8090.8 < 2.2e-16 ***
# Ecoregion_NARS:Agency      16.1184  8 4372.0 < 2.2e-16 ***
# YearCont:landuse            0.8323  3 7286.8 0.4758740    
# YearCont:Ecoregion_NARS     0.7700  8 7625.3 0.6293273         

dat1 <- data.frame(ggemmeans(m2, terms = c("YearCont[1:27]","landuse","Agency"), 
                             type = "fixed",back.transform=TRUE))

# filter so effects on span the length of which we have data
min_max <- dat_abun %>%
  group_by(Agency, landuse) %>%
  summarize(min = min(Year),
            max = max(Year))

dat1 <- dat1 %>% 
  left_join(min_max, by = c("group" = "landuse",
                            "facet" = "Agency")) %>%
  filter(x >= min,
         x <= max)

m2_LU_fig <- ggplot(data = dat1, 
       aes(x = x + 1992, y = predicted)) +
  facet_wrap(~facet) +
  geom_line(aes(color = group), size = 1) +
  geom_ribbon(aes(x=x+ 1992, ymin = conf.low , ymax = conf.high, fill = group), 
              alpha=0.1) +
  scale_color_brewer(palette="Dark2", name = "Dominant landuse",
                     limits = c("Forest_Wetland","Agriculture","Urban","Other"),
                     breaks = c("Forest_Wetland","Other","Agriculture","Urban"),
                     labels = c("Forest/wetland","Grassland/shrub","Agriculture","Urban")) +
  scale_fill_brewer(palette="Dark2", name = "Dominant landuse",
                    limits = c("Forest_Wetland","Agriculture","Urban","Other"),
                    breaks = c("Forest_Wetland","Other","Agriculture","Urban"),
                    labels = c("Forest/wetland","Grassland/shrub","Agriculture","Urban")) +
  scale_x_continuous(breaks=c(1993,1998,2003,2008,2013, 2018))+
  ylab(expression("Total density (ind. / m"^2*")")) +
  xlab("Year") +
  theme_bw() +
  theme(text = element_text(size = 18))

###                            ###
### SITE-LEVEL ALPHA DIVERSITY ###    
###                            ###

### rarefied to 300 at site level
dat_alpha_all <- read.csv("./insects_only/data/dat_alpha_all.csv", row.names=1, header=T)

### doesn't converge but leaving so model same form as all macroinvertebrates
m_LU_a <- glmmTMB::glmmTMB(alpha ~ Year + landuse  + Ecoregion_NARS +
             PropID + AreaSampTot_m2 + Agency + Gen_ID_Prop +
             Ecoregion_NARS:Agency +
             Year:landuse + 
             Year:Ecoregion_NARS +
             (1|SiteNumber),
             family = nbinom2,
           data = dat_alpha_all)

performance::r2(m_LU_a)
# R2 for Mixed Models
# 
# Conditional R2: 0.634
# Marginal R2: 0.335

summary(m_LU_a)

car::Anova(m_LU_a, type = "III") 
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: alpha
#                           Chisq Df Pr(>Chisq)    
# (Intercept)           2018.8264  1  < 2.2e-16 ***
# Year                     0.0229  1  0.8797078    
# landuse                 94.4901  3  < 2.2e-16 ***
# Ecoregion_NARS         143.5091  8  < 2.2e-16 ***
# PropID                   8.1925  1  0.0042064 ** 
# AreaSampTot_m2          39.1201  1  3.985e-10 ***
# Agency                   2.5557  1  0.1098967    
# Gen_ID_Prop            194.8131  1  < 2.2e-16 ***
# Ecoregion_NARS:Agency   69.5664  8  5.995e-12 ***
# Year:landuse             6.6773  3  0.0829258 .  
# Year:Ecoregion_NARS     27.0006  8  0.0007068 ***

# pair wise differences of landuse
pairs(emmeans(m_LU_a, "landuse"), adjust = "none")
# contrast                     estimate     SE   df t.ratio p.value
# Agriculture - Forest_Wetland  -0.1224 0.0131 7630  -9.372  <.0001
# Agriculture - Other           -0.0482 0.0163 7630  -2.950  0.0032
# Agriculture - Urban            0.1437 0.0176 7630   8.182  <.0001
# Forest_Wetland - Other         0.0743 0.0144 7630   5.166  <.0001
# Forest_Wetland - Urban         0.2661 0.0159 7630  16.766  <.0001
# Other - Urban                  0.1918 0.0188 7630  10.192  <.0001
#
# make a figure
dat1 <- data.frame(ggemmeans(m_LU_a, terms = c("Year[1:27]","landuse","Agency"), 
                             type = "fixed",back.transform=TRUE))

# filter so effects on span the length of which we have data
min_max <- dat_alpha_all %>%
  group_by(Agency, landuse) %>%
  summarize(min = min(Year),
            max = max(Year))

dat1 <- dat1 %>% 
  left_join(min_max, by = c("group" = "landuse",
                            "facet" = "Agency")) %>%
  filter(x >= min,
         x <= max)

m_LU_a_fig <- ggplot(data = dat1, 
                    aes(x = x + 1992, y = predicted)) +
  facet_wrap(~facet) +
  geom_line(aes(color = group), size = 1) +
  geom_ribbon(aes(x=x+ 1992, ymin = conf.low , ymax = conf.high, fill = group), 
              alpha=0.1) +
  scale_color_brewer(palette="Dark2", name = "Dominant landuse",
                     limits = c("Forest_Wetland","Agriculture","Urban","Other"),
                     breaks = c("Forest_Wetland","Other","Agriculture","Urban"),
                     labels = c("Forest/wetland","Grassland/shrub","Agriculture","Urban")) +
  scale_fill_brewer(palette="Dark2", name = "Dominant landuse",
                    limits = c("Forest_Wetland","Agriculture","Urban","Other"),
                    breaks = c("Forest_Wetland","Other","Agriculture","Urban"),
                    labels = c("Forest/wetland","Grassland/shrub","Agriculture","Urban")) +
  scale_x_continuous(breaks=c(1993,1998,2003,2008,2013, 2018))+
  ylab(expression(atop(""*alpha*" diversity","(richness of genera at sites, rarefied)"))) +
  xlab("Year") +
  theme_bw() +
  theme(text = element_text(size = 18))

m_LU_a_fig


###       ###
### GAMMA ###
###       ###

# divide one by CI of estimate of gamma
dat_div_all_LU_Ag$wt.sc.CL <- 1 / dat_div_all_LU_Ag$sc.CL
# standardize to mean to avoid small range of values
dat_div_all_LU_Ag$wt.std <- dat_div_all_LU_Ag$wt.sc.CL/ mean(dat_div_all_LU_Ag$wt.sc.CL)


m2g_LU <- lmer(sc.gamma ~ Year + landuse + Agency + Gen_ID_Prop +
                   Year:landuse + 
                   (1 | Ecoregion_NARS:landuse) + 
                   (Agency + 0 |Ecoregion_NARS),
                 weights = wt.std,
                 data = dat_div_all_LU_Ag,
                 control = lmerControl(optimizer = "nlminbwrap",
                                     optCtrl = list(maxfun = 2e7)),
                 REML = TRUE)

summary(m2g_LU)

slm2g <- emtrends(m2g_LU, specs =  ~ landuse, var = "Year")
slm2g
#  landuse        Year.trend    SE  df lower.CL upper.CL
# Agriculture       -0.0650 0.127 530   -0.314   0.1837
# Forest_Wetland     0.0386 0.114 580   -0.184   0.2616
# Other              0.1282 0.166 372   -0.199   0.4549
# Urban             -0.3842 0.147 383   -0.674  -0.0945

performance::r2(m2g_LU) 
# R2 for Mixed Models
# 
# Conditional R2: 0.696
#    Marginal R2: 0.614

summary(m2g_LU)

car::Anova(m2g_LU, test = "F", type = "III")
# Analysis of Deviance Table (Type III Wald F tests with Kenward-Roger df)
# 
# Response: sc.gamma
#                    F Df Df.res    Pr(>F)    
# (Intercept)  68.3707  1  58.35 2.084e-11 ***
# Year          0.2639  1 530.35   0.60770    
# landuse       3.7387  3  47.18   0.01720 *  
# Agency       56.8017  1   8.29 5.553e-05 ***
# Gen_ID_Prop   2.1581  1 373.02   0.14266    
# Year:landuse  2.8275  3 447.20   0.03821 *  

pairs(emtrends(m2g_LU, "landuse", var = "Year"), adjust = "none")
# contrast                     estimate    SE  df t.ratio p.value
# Agriculture - Forest_Wetland  -0.1037 0.150 578  -0.692  0.4891
# Agriculture - Other           -0.1933 0.196 413  -0.988  0.3240
# Agriculture - Urban            0.3191 0.174 443   1.829  0.0681
# Forest_Wetland - Other        -0.0896 0.188 402  -0.476  0.6345
# Forest_Wetland - Urban         0.4228 0.163 465   2.600  0.0096 #
# Other - Urban                  0.5124 0.209 364   2.453  0.0146 #

# make a figure
dat1 <- data.frame(ggpredict(m2g_LU, terms = c("Year[1:27]","landuse","Agency"), 
                             type = "fixed"))

# filter dat1, so effects on span the length of which we have data
min_max <- dat_div_all_LU_Ag %>%
  group_by(landuse, Agency) %>%
  summarize(min = min(Year),
            max = max(Year))

dat1 <- dat1 %>% 
  left_join(min_max, by = c("group" = "landuse",
                            "facet" = "Agency")) %>%
  filter(x >= min,
         x <= max) 

m2g_LU_fig <-ggplot(data = dat1, 
       aes(x = x + 1992, y = predicted)) +
  facet_wrap(~facet) +
  geom_line(aes(color = group), size = 1) +
  geom_ribbon(aes(x=x+ 1992, ymin = conf.low , ymax = conf.high, fill = group), 
              alpha=0.1) +
  scale_color_brewer(palette="Dark2", name = "Dominant landuse",
                     limits = c("Forest_Wetland","Agriculture","Urban","Other"),
                     breaks = c("Forest_Wetland","Other","Agriculture","Urban"),
                     labels = c("Forest/wetland","Grassland/shrub","Agriculture","Urban")) +
  scale_fill_brewer(palette="Dark2", name = "Dominant landuse",
                    limits = c("Forest_Wetland","Agriculture","Urban","Other"),
                    breaks = c("Forest_Wetland","Other","Agriculture","Urban"),
                    labels = c("Forest/wetland","Grassland/shrub","Agriculture","Urban")) +
  scale_x_continuous(breaks=c(1993,1998,2003,2008,2013,2018))+
  ylab(expression(atop(gamma[est]*" diversity",
                       "(total no. of genera within ecoregions)")))+
  xlab("Year") +
  theme_bw() +
  theme(text = element_text(size = 18)) 

m2g_LU_fig

###       ###
### BETA  ###
###       ###

#### recalulate beta
dat_div_all_LU_Ag$prop.beta.sc <- dat_div_all_LU_Ag$alpha.bar / dat_div_all_LU_Ag$sc.gamma

# make a figure
m2b_LU <- lmer(prop.beta.sc ~ Year + landuse + Agency + Gen_ID_Prop +
                   Year:landuse + 
                   (1 | Ecoregion_NARS:landuse) +
                   (Agency + 0 |Ecoregion_NARS),
                 weights = wt.std,
                 data = dat_div_all_LU_Ag,
                 control = lmerControl(optimizer = "nlminbwrap",
                                     optCtrl = list(maxfun = 2e5)),
                 REML = TRUE)

performance::r2(m2b_LU)
# R2 for Mixed Models
# 
# Conditional R2: 0.463
#    Marginal R2: 0.291

car::Anova(m2b_LU, test = "F", type = "III")
# Analysis of Deviance Table (Type III Wald F tests with Kenward-Roger df)
# 
# Response: prop.beta.sc
#                    F Df Df.res    Pr(>F)    
# (Intercept)  45.0038  1 216.44 1.692e-10 ***
# Year          0.1772  1 529.21    0.6740    
# landuse       1.7533  3  61.62    0.1655    
# Agency       52.2837  1   8.69 5.885e-05 ***
# Gen_ID_Prop   1.1991  1 389.53    0.2742    
# Year:landuse  1.4625  3 454.33    0.2241        

slm2b <- emtrends(m2b_LU, specs =  ~ landuse, var = "Year")
slm2b
# landuse        Year.trend      SE  df lower.CL upper.CL
# Agriculture     -0.000586 0.00139 529 -0.00332  0.00215
# Forest_Wetland   0.001046 0.00124 597 -0.00139  0.00349
# Other           -0.002191 0.00182 384 -0.00576  0.00138
# Urban            0.002067 0.00161 387 -0.00111  0.00524

dat1 <- data.frame(ggpredict(m2b_LU, terms = c("Year[1:27]","landuse","Agency"), 
                             type = "fixed"))

# filter dat1, so effects on span the length of which we have data
dat1 <- dat1 %>% 
  left_join(min_max, by = c("group" = "landuse",
                            "facet" = "Agency")) %>%
  filter(x >= min,
         x <= max) 

m2b_LU_fig <- ggplot(data = dat1, 
       aes(x = x + 1992, y = predicted)) +
  facet_wrap(~facet) +
  geom_line(aes(color = group), size = 1) +  
  geom_ribbon(aes(x=x+ 1992, ymin = conf.low , ymax = conf.high, fill = group), 
              alpha=0.1) +
  scale_color_brewer(palette="Dark2", name = "Dominant landuse",
                     limits = c("Forest_Wetland","Agriculture","Urban","Other"),
                     breaks = c("Forest_Wetland","Other","Agriculture","Urban"),
                     labels = c("Forest/wetland","Grassland/shrub","Agriculture","Urban")) +
  scale_fill_brewer(palette="Dark2", name = "Dominant landuse",
                    limits = c("Forest_Wetland","Agriculture","Urban","Other"),
                    breaks = c("Forest_Wetland","Other","Agriculture","Urban"),
                    labels = c("Forest/wetland","Grassland/shrub","Agriculture","Urban")) +
  scale_x_continuous(breaks=c(1993,1998,2003,2008,2013,2018))+
  ylab(expression(atop(""*beta*" diversity","(spatial turnover within ecoregions)"))) +
  xlab("Year") +
  theme_bw() +
  theme(text = element_text(size = 18)) 

m2b_LU_fig

# join panels together
A = m2_LU_fig + theme(legend.position = "none")
B = m_LU_a_fig + theme(legend.position = "none")
C = m2g_LU_fig + theme(legend.position = "none")
D = m2b_LU_fig + theme(legend.position = "none")
legend = cowplot::get_legend(m2b_LU_fig + theme(legend.position = "top",
                                   text = element_text(size = 25)) +
                               guides(colour = guide_legend(nrow = 1)))

f2 <- cowplot::plot_grid(A, B, C, D, labels = c("A", "B", "C", "D"),
                         label_size = 20,
                         ncol = 2,
                         nrow = 2,
                         align = "v")

f2_legend <- cowplot::plot_grid(legend, f2, ncol = 1, rel_heights = c(1,15))

ggsave("./insects_only/figures/div_time_landuse.png", f2_legend, dpi = 200, 
       height = 10.5, width = 16)
