# Objective: Diversity analyses including random effect of agency with linear 
# models - using only insects

# load packages
library(tidyverse)
library(lme4)
library(ggeffects)
library(emmeans)

# read in datasets

# abubndance/density
dat_abun <- read.csv("./insects_only/data/dat_abun_all.csv", row.names = 1, header =T) 

# diversity metrics
# filter out observations with in which sample coverage gamma was not 
# estimated between 0.7 and 0.8 
dat_div_NRSA_USGS <- read.csv("./insects_only/data/dat_div_NRSA_USGS.csv", header = T,
                              row.names = 1) %>%
  filter(sc >= 0.7 & sc <= 0.8)

# create a variable for model weights for gamma and beta models
# the inverse of the CI estimate of gamma, standardized to the mean

# divide one by CI of estimate of gamma
dat_div_NRSA_USGS$wt.sc.CL <- 1 / dat_div_NRSA_USGS$sc.CL
# standardize to mean to avoid small range of values
dat_div_NRSA_USGS$wt.std <- dat_div_NRSA_USGS$wt.sc.CL/ mean(dat_div_NRSA_USGS$wt.sc.CL)


###         ###
### DENSITY ###
###         ###


# create a model
m2 <- lmer(log(tot_abun+1) ~ YearCont + Ecoregion_NARS + YearCont:Ecoregion_NARS +
             PropID + AreaSampTot_m2 + Agency + Gen_ID_Prop +
             Ecoregion_NARS:Agency + 
             (1|SiteNumber),
           data = dat_abun,
           REML = TRUE)

sl2 <- emtrends(m2, specs = ~ YearCont, var = 'YearCont', 
                  pbkrtest.limit = 8236,
                  lmerTest.limit = 8236)
summary(sl2)
# YearCont YearCont.trend     SE   df lower.CL upper.CL
#     16.6       -0.00986 0.00152 7837  -0.0128 -0.00687


#YearCont trend is exp(-0.00986) for year 1 year increase
(exp(-0.00986)-1) * 100
# So, every year abundance decreases by 0.98%

(exp(27*-0.00981155))-1
# 23.3% reductions in total density over 27 years

car::Anova(m2, type = "II", test = "F")
# Analysis of Deviance Table (Type II Wald F tests with Kenward-Roger df)
# 
# Response: log(tot_abun + 1)
#                                  F Df Df.res    Pr(>F)    
# YearCont                   49.6104  1 7705.9 2.037e-12 ***
# Ecoregion_NARS             22.5739  8 4491.4 < 2.2e-16 ***
# PropID                  14242.7363  1 7999.4 < 2.2e-16 ***
# AreaSampTot_m2           1589.8427  1 7518.1 < 2.2e-16 ***
# Agency                     17.2307  1 4657.5 3.370e-05 ***
# Gen_ID_Prop              1381.5556  1 8093.9 < 2.2e-16 ***
# YearCont:Ecoregion_NARS     0.7541  8 7672.8    0.6435    
# Ecoregion_NARS:Agency      16.3199  8 4383.4 < 2.2e-16 ***

summary(m2)

performance::r2(m2)
#Conditional R2: 0.793
#   Marginal R2: 0.737

# get overall effect of YearCont
dat1 <- data.frame(ggemmeans(m2, terms = "YearCont[1:27]", type = "fixed", back.transform = TRUE))
# get effects of YearCont across ERs and Agencies 
dat2 <- data.frame(ggemmeans(m2, terms = c("YearCont[1:27]", 
                                           "Ecoregion_NARS","Agency"), 
                             type = "fixed", back.transform = TRUE))
# filter dat2, so effects on span the length of which we have data
min_max <- dat_abun %>%
  group_by(Ecoregion_NARS, Agency) %>%
  summarize(min = min(Year),
            max = max(Year))

dat2 <- dat2 %>% 
  left_join(min_max, by = c("group" = "Ecoregion_NARS",
                            "facet" = "Agency")) %>%
  filter(x >= min,
         x <= max)

m2_fig <- ggplot() +
  geom_line(data=dat1, 
            aes(x = x + 1992, y = predicted), 
            color = "black", size=2) +
  geom_ribbon(data=dat1, 
              aes(x=x + 1992, ymin = conf.low, ymax = conf.high),
              #aes(x=x + 1992, ymin = predicted - std.error, ymax = predicted + std.error),
              alpha=0.3) +
  geom_line(data = dat2, 
            aes(x = x + 1992, y = predicted, color = group),
            linetype = "dashed", size = 1) +
  facet_wrap(~facet) +
  scale_color_brewer(palette="Set1", name = "Ecoregion") +
  scale_x_continuous(breaks=c(1993,1998,2003,2008,2013,2018))+
  scale_y_continuous(breaks=c(250,500,750,1000,1250))+
  ylab(expression("Total density (ind. / m"^2*")")) +
  xlab("Year") +
  theme_bw() +
  theme(text = element_text(size = 18)) 

m2_fig

###                            ###
### SITE-LEVEL ALPHA DIVERSITY ###    
###                            ###

### rarefied to 300 at site level
dat_alpha_all_Sn <- read.csv("./insects_only/data/dat_alpha_all.csv", row.names=1, header=T)

m_Sn <- glmmTMB::glmmTMB(alpha ~ Year + Ecoregion_NARS + Year:Ecoregion_NARS +
                PropID + AreaSampTot_m2 + Agency + Gen_ID_Prop +
                Ecoregion_NARS:Agency + 
             (1|SiteNumber),
             family = poisson,
             data = dat_alpha_all_Sn)

sl1 <- emtrends(m_Sn, specs = ~ Year, var = 'Year',
                pbkrtest.limit = 8240,
                lmerTest.limit = 8240)
summary(sl1)
#Year Year.trend       SE   df lower.CL upper.CL
#15.3   -0.00263 0.000679 7637 -0.00396 -0.00129

(exp(-0.00263)-1) * 100
# [1] -0.2626545
# So, every year richness decreases by 0.26%

((exp(27*-0.002626545))-1) *100
# [1] -6.846053
# So, across 27 year richness decreased 6.8%

car::Anova(m_Sn, type = "II") 
# Response: alpha
#                          Chisq Df Pr(>Chisq)    
# Year                    17.774  1  2.487e-05 ***
# Ecoregion_NARS        1034.425  8  < 2.2e-16 ***
# PropID                  13.446  1  0.0002455 ***
# AreaSampTot_m2          41.461  1  1.203e-10 ***
# Agency                 448.082  1  < 2.2e-16 ***
# Gen_ID_Prop            209.655  1  < 2.2e-16 ***
# Year:Ecoregion_NARS     25.040  8  0.0015303 ** 
# Ecoregion_NARS:Agency   85.610  8  3.606e-15 ***

performance::r2(m_Sn)
# R2 for Mixed Models
# 
# Conditional R2: 0.672
#    Marginal R2: 0.293

# get overall effect of YearCont
dat1 <- data.frame(ggemmeans(m_Sn, terms = "Year[1:27]", type = "fixed", back.transform = TRUE))
# get effects of YearCont across ERs and Agencies 
dat2 <- data.frame(ggemmeans(m_Sn, terms = c("Year[1:27]", 
                                           "Ecoregion_NARS","Agency"), 
                             type = "fixed", back.transform = TRUE))
# filter dat2, so effects on span the length of which we have data
min_max <- dat_alpha_all_Sn %>%
  group_by(Ecoregion_NARS, Agency) %>%
  summarize(min = min(Year),
            max = max(Year))

dat2 <- dat2 %>% 
  left_join(min_max, by = c("group" = "Ecoregion_NARS",
                            "facet" = "Agency")) %>%
  filter(x >= min,
         x <= max)

m3_fig_a <- ggplot() +
  geom_line(data=dat1, 
            aes(x = x + 1992, y = predicted), 
            color = "black", size=2) +
  geom_ribbon(data=dat1, 
              aes(x=x + 1992, ymin = conf.low, ymax = conf.high),
              #aes(x=x + 1992, ymin = predicted - std.error, ymax = predicted + std.error),
              alpha=0.3) +
  geom_line(data = dat2, 
            aes(x = x + 1992, y = predicted, color = group),
            linetype = "dashed", size = 1) +
  facet_wrap(~facet) +
  scale_color_brewer(palette="Set1", name = "Ecoregion") +
  scale_x_continuous(breaks=c(1993,1998,2003,2008,2013,2018))+
  scale_y_continuous(breaks=c(15,20,25,30,35,40), limits = c(15,40))+
  ylab(expression(atop(""*alpha*" diversity","(richness of genera at sites, rarefied)"))) +
  xlab("Year") +
  theme_bw() +
  theme(text = element_text(size = 18)) 

m3_fig_a


###                 ###
### GAMMA DIVERSITY ###
###                 ###

m1g <- lmer(sc.gamma ~ Year + Gen_ID_Prop +
              (Year + 0 |Ecoregion_NARS) +
              (1|Ecoregion_NARS:Agency), 
            weights = wt.std,
            data = dat_div_NRSA_USGS,
            control=lmerControl(optimizer="nlminbwrap",
                                optCtrl=list(maxfun=2e5)),
            REML = T) 

summary(m1g)
# Fixed effects:
#             Estimate Std. Error t value
# (Intercept)  32.6745     7.5111   4.350
# Year         -0.2194     0.1326  -1.655
# Gen_ID_Prop  36.0081     9.8626   3.651
# slope of year is -0.2194, so for each year we loose 0.22 genera

car::Anova(m1g, test = "F", type = "II") 
# Analysis of Deviance Table (Type II Wald F tests with Kenward-Roger df)
# 
# Response: sc.gamma
#                   F Df  Df.res    Pr(>F)    
# Year         2.7385  1  17.911 0.1153740    
# Gen_ID_Prop 13.1230  1 171.883 0.0003838 ***

performance::r2(m1g)
# R2 for Mixed Models
# 
# Conditional R2: 0.797
#    Marginal R2: 0.019

dat1 <- data.frame(ggpredict(m1g, terms = "Year [1:27]", type = "fixed"))
dat2 <- data.frame(ggpredict(m1g, terms = c("Year [1:27]", "Ecoregion_NARS",
                                            "Agency"), 
                             type = "random"))

dat2 <- dat2 %>% 
  left_join(min_max, by = c("group" = "Ecoregion_NARS",
                            "facet" = "Agency")) %>%
  filter(x >= min,
         x <= max) 

m1g_fig <- ggplot() +
  geom_line(data=dat1, aes(x = x + 1992, y = predicted), color = "black", size=2) +
  geom_ribbon(data=dat1, aes(x=x+ 1992, ymin = conf.low , ymax = conf.high),
              alpha=0.3) +
  geom_line(data = dat2, 
            aes(x = x + 1992, y = predicted, color = group),
            linetype = "dashed", size = 1) +
  facet_wrap(~facet) +
  scale_color_brewer(palette="Set1", name = "Ecoregion") +
  scale_x_continuous(breaks=c(1993,1998,2003,2008,2013,2018))+
  ylab(expression(atop(gamma[est]*" diversity",
                       "(total no. of genera within ecoregions)")))+
  xlab("Year") +
  theme_bw() +
  theme(text = element_text(size = 18))

m1g_fig

###                ###
### BETA DIVERSITY ###   
###                ###

#### recalulate beta
dat_div_NRSA_USGS$prop.beta.sc <- dat_div_NRSA_USGS$alpha.bar / dat_div_NRSA_USGS$sc.gamma

m1b <- lmer(prop.beta.sc ~ Year + Gen_ID_Prop +
                (Year + 0 |Ecoregion_NARS) +
                (1|Ecoregion_NARS:Agency), 
              weights = wt.std,
              data = dat_div_NRSA_USGS,
            control=lmerControl(optimizer="nlminbwrap",
                                optCtrl=list(maxfun=2e5)),
             REML = TRUE)   

car::Anova(m1b, test = "F",type = "II") 
# Analysis of Deviance Table (Type II Wald F tests with Kenward-Roger df)
# 
# Response: prop.beta.sc
#                  F Df Df.res   Pr(>F)   
# Year        1.4155  1  12.56 0.256148   
# Gen_ID_Prop 9.4269  1 175.18 0.002479 **

performance::r2(m1b)
# R2 for Mixed Models
# 
# Conditional R2: 0.713
#    Marginal R2: 0.018

dat1 <- data.frame(ggpredict(m1b, terms = "Year[1:27]", type = "fixed"))
dat2 <- data.frame(ggpredict(m1b, terms = c("Year[1:27]", "Ecoregion_NARS","Agency"), 
                             type = "random"))

dat2 <- dat2 %>% 
  left_join(min_max, by = c("group" = "Ecoregion_NARS",
                            "facet" = "Agency")) %>%
  filter(x >= min,
         x <= max) 

m1b_fig <- ggplot() +
  geom_line(data=dat1, aes(x = x + 1992, y = predicted), color = "black", size=2) +
  geom_ribbon(data=dat1, aes(x = x+ 1992, ymin = conf.low , ymax = conf.high),
              alpha=0.3) +
  geom_line(data = dat2, 
            aes(x = x + 1992, y = predicted, color = group),
            linetype = "dashed", size = 1) +
  facet_wrap(~facet) +
  scale_color_brewer(palette="Set1", name = "Ecoregion") +
  scale_x_continuous(breaks=c(1993,1998,2003,2008,2013,2018))+
  #scale_y_continuous(breaks=c(0.3,0.4,0.5))+
  ylab(expression(atop(""*beta*" diversity","(spatial turnover within ecoregions)"))) +
  xlab("Year") +
  theme_bw() +
  theme(text = element_text(size = 18))

m1b_fig

# join panels together - main text Figure 2
A = m2_fig + theme(legend.position = "none")
B = m3_fig_a + theme(legend.position = "none")
C = m1g_fig + theme(legend.position = "none")
D = m1b_fig + theme(legend.position = "none")
legend = cowplot::get_legend(m1b_fig + 
                               theme(legend.position = "top",
                                     text = element_text(size = 25)) +
                               guides(colour = guide_legend(nrow = 1)))

f2 <- cowplot::plot_grid(A, B, C, D, labels = c("A", "B", "C", "D"),
                         label_size = 20,
                         ncol = 2,
                         nrow = 2,
                         align = "v")
f2_legend <- cowplot::plot_grid(legend, f2, ncol = 1, rel_heights = c(1,15))

ggsave("./insects_only/figures/div_time.png", f2_legend, dpi = 200, 
       height = 10.5, width = 16)
