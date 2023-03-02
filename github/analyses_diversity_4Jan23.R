# Objective: Diversity analyses including random effect of agency with linear 
# models

# load packages
library(tidyverse)
library(lme4)
library(ggeffects)
library(emmeans)

# read in datasets

# abundance/density of all Id'd to phylum
dat_phy_abun <- read.csv("./data/dat_phy_abun_all.csv", row.names = 1, header =T)

# abubndance/density
dat_abun <- read.csv("./data/dat_abun_all.csv", row.names = 1, header =T) 

# diversity metrics
# filter out observations with in which sample coverage gamma was not 
# estimated between 0.7 and 0.8 
dat_div_NRSA_USGS <- read.csv("./data/dat_div_NRSA_USGS_4Jan.csv", header = T,
                              row.names = 1) %>%
  filter(sc >= 0.7 & sc <= 0.8)

# create a variable for model weights for gamma and beta models
# the inverse of the CI estimate of gamma, standardized to the mean

# divide one by CI of estimate of gamma
dat_div_NRSA_USGS$wt.sc.CL <- 1 / dat_div_NRSA_USGS$sc.CL
# standardize to mean to avoid small range of values
dat_div_NRSA_USGS$wt.std <- dat_div_NRSA_USGS$wt.sc.CL/ mean(dat_div_NRSA_USGS$wt.sc.CL)

###                 ###
### DENSITY FOR ALL ###
###                 ###

# create a model
m1 <- lmer(log(tot_abun) ~ YearCont + Ecoregion_NARS + YearCont:Ecoregion_NARS +
             PropID + AreaSampTot_m2 + Agency + 
             Ecoregion_NARS:Agency + 
             (1|SiteNumber),
           data = dat_phy_abun,
           REML = TRUE)

sl1 <- emtrends(m1, specs = ~ YearCont, var = 'YearCont',
                 pbkrtest.limit = 8240,
                 lmerTest.limit = 8240)
summary(sl1) 
# YearCont YearCont.trend      SE   df lower.CL upper.CL
#      16.6       -0.00179 0.0013 7597 -0.00433  0.00076    

#YearCont trend is exp(-0.00179) for year 1 year increase
(exp(-0.00179)-1) * 100
# So, every year abundance decreases by 0.18%

(exp(27*-0.001788399))-1
# [1] -0.04713951
# 4.7% reductions in total density over 27 years

car::Anova(m1, type = "II", test = "F")
# Analysis of Deviance Table (Type II Wald F tests with Kenward-Roger df)
# 
# Response: log(tot_abun)
#                                  F Df Df.res    Pr(>F)    
# YearCont                    0.8543  1 7463.3    0.3554    
# Ecoregion_NARS              9.9749  8 4322.0 7.328e-14 ***
# PropID                  19426.4840  1 7903.3 < 2.2e-16 ***
# AreaSampTot_m2           2253.2057  1 7007.4 < 2.2e-16 ***
# Agency                      5.8299  1 4258.1    0.0158 *  
# YearCont:Ecoregion_NARS     1.3005  8 7481.6    0.2380    
# Ecoregion_NARS:Agency      15.3048  8 3963.9 < 2.2e-16 ***

summary(m1)

performance::r2(m1)
# R2 for Mixed Models
# 
#  Conditional R2: 0.812
#     Marginal R2: 0.774

# get overall effect of YearCont
dat1 <- data.frame(ggemmeans(m1, terms = "YearCont[1:27]", type = "fixed", back.transform = TRUE))
# get effects of YearCont across ERs and Agencies 
dat2 <- data.frame(ggemmeans(m1, terms = c("YearCont[1:27]", 
                                           "Ecoregion_NARS","Agency"), 
                             type = "fixed", back.transform = TRUE))
# filter dat2, so effects on span the length of which we have data
min_max <- dat_phy_abun %>%
  group_by(Ecoregion_NARS, Agency) %>%
  summarize(min = min(YearCont),
            max = max(YearCont))

dat2 <- dat2 %>% 
  left_join(min_max, by = c("group" = "Ecoregion_NARS",
                            "facet" = "Agency")) %>%
  filter(x >= min,
         x <= max)

m1_fig <- ggplot() +
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
  scale_y_continuous(breaks=c(750,1000,1250,1500,1750), limits=c(670,1750))+
  ylab(expression("Total density (ind. / m"^2*")")) +
  xlab("Year") +
  theme_bw() +
  theme(text = element_text(size = 18)) 

m1_fig


###         ###
### DENSITY ###
###         ###


# create a model
m2 <- lmer(log(tot_abun) ~ YearCont + Ecoregion_NARS + YearCont:Ecoregion_NARS +
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
#     16.6       -0.00432 0.00136 7713 -0.00699 -0.00165

#YearCont trend is exp(-0.0043) for year 1 year increase
(exp(-0.00432)-1) * 100
# So, every year abundance decreases by 0.43%

(exp(27*-0.004310682))-1
# [1] -0.1098706
# 11% reductions in total density over 27 years

car::Anova(m2, type = "II", test = "F") 
# Analysis of Deviance Table (Type II Wald F tests with Kenward-Roger df)
# 
# Response: log(tot_abun)
#                                  F Df Df.res    Pr(>F)    
# YearCont                    8.4181  1 7585.9  0.003726 ** 
# Ecoregion_NARS             11.2323  8 4361.7 7.399e-16 ***
# PropID                  18764.8864  1 7914.7 < 2.2e-16 ***
# AreaSampTot_m2           2169.7149  1 7173.3 < 2.2e-16 ***
# Agency                      0.2311  1 4368.1  0.630721    
# Gen_ID_Prop              1310.3042  1 8022.0 < 2.2e-16 ***
# YearCont:Ecoregion_NARS     1.0119  8 7530.5  0.424323    
# Ecoregion_NARS:Agency      14.6496  8 4068.7 < 2.2e-16 ***

summary(m2)

performance::r2(m2)
# R2 for Mixed Models
# 
# Conditional R2: 0.813
#    Marginal R2: 0.772

# get overall effect of YearCont
dat1 <- data.frame(ggemmeans(m2, terms = "YearCont[1:27]", type = "fixed", back.transform = TRUE))
# get effects of YearCont across ERs and Agencies 
dat2 <- data.frame(ggemmeans(m2, terms = c("YearCont[1:27]", 
                                           "Ecoregion_NARS","Agency"), 
                             type = "fixed", back.transform = TRUE))
# filter dat2, so effects on span the length of which we have data
min_max <- dat_abun %>%
  group_by(Ecoregion_NARS, Agency) %>%
  summarize(min = min(YearCont),
            max = max(YearCont))

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

# rarefied to 300 at site level
dat_alpha_all_Sn <- read.csv("./data/dat_alpha_all_4Jan.csv", row.names=1, header=T)

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
#15.3    0.00426 0.000627 7637  0.00303  0.00549

(exp(0.00426)-1) * 100
# [1] 0.4269087
# So, every year richness increases by 0.43%

((exp(27*0.004269087))-1) *100
# [1] 12.21712
# So, across 27 year richness increases 12.2%

car::Anova(m_Sn, type = "II") 
#                          Chisq Df Pr(>Chisq)    
# Year                   68.0107  1  < 2.2e-16 ***
# Ecoregion_NARS        930.7046  8  < 2.2e-16 ***
# PropID                  7.6081  1  0.0058106 ** 
# AreaSampTot_m2         41.4257  1  1.224e-10 ***
# Agency                715.0728  1  < 2.2e-16 ***
# Gen_ID_Prop           127.3458  1  < 2.2e-16 ***
# Year:Ecoregion_NARS    28.8652  8  0.0003348 ***
# Ecoregion_NARS:Agency  56.1487  8  2.639e-09 ***

performance::r2(m_Sn)
# R2 for Mixed Models
# 
# Conditional R2: 0.682
#    Marginal R2: 0.337

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

### not rarefied - dataset matched to sites-dates to rarefied dataset
dat_alpha_all_S <- read.csv("./data/dat_alpha_all_S4Jan.csv", row.names=1, header=T)

m_S <- glmmTMB::glmmTMB(alpha ~ Year + Ecoregion_NARS + Year:Ecoregion_NARS +
               PropID + AreaSampTot_m2 + Agency + Gen_ID_Prop +
               Ecoregion_NARS:Agency + 
               (1|SiteNumber),
             family = poisson,
             data = dat_alpha_all_S)

sl1 <- emtrends(m_S, specs = ~ Year, var = 'Year',
                pbkrtest.limit = 8240,
                lmerTest.limit = 8240)
summary(sl1)
# Year Year.trend       SE   df lower.CL upper.CL
# 15.3     0.0033 0.000616 7576  0.00209   0.0045


car::Anova(m_S, type = "II")
#Response: alpha
#                           Chisq Df Pr(>Chisq)    
# Year                    49.5947  1  1.890e-12 ***
# Ecoregion_NARS         916.3495  8  < 2.2e-16 ***
# PropID                   0.0652  1     0.7985    
# AreaSampTot_m2          22.9781  1  1.639e-06 ***
# Agency                1318.1461  1  < 2.2e-16 ***
# Gen_ID_Prop             84.3433  1  < 2.2e-16 ***
# Year:Ecoregion_NARS     65.9492  8  3.132e-11 ***
# Ecoregion_NARS:Agency   81.3924  8  2.563e-14 ***

performance::r2(m_S)
# R2 for Mixed Models
# 
# Conditional R2: 0.744
#    Marginal R2: 0.391


# get overall effect of YearCont
dat1 <- data.frame(ggemmeans(m_S, terms = "Year[1:27]", type = "fixed", back.transform = TRUE))
# get effects of YearCont across ERs and Agencies 
dat2 <- data.frame(ggemmeans(m_S, terms = c("Year[1:27]", 
                                             "Ecoregion_NARS","Agency"), 
                             type = "fixed", back.transform = TRUE))
# filter dat2, so effects on span the length of which we have data
min_max <- dat_alpha_all_S %>%  
  group_by(Ecoregion_NARS, Agency) %>%
  summarize(min = min(Year),
            max = max(Year))

dat2 <- dat2 %>% 
  left_join(min_max, by = c("group" = "Ecoregion_NARS",
                            "facet" = "Agency")) %>%
  filter(x >= min,
         x <= max)

m3_fig_b <- ggplot() +
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
  ylab(expression(atop(""*alpha*" diversity","(richness of genera at sites, not rarefied)"))) +
  xlab("Year") +
  theme_bw() +
  theme(text = element_text(size = 18)) 

m3_fig_b


###                 ###
### ALPHA DIVERSITY ###
###                 ###

m1a <- lmer(alpha.bar ~ Year + Gen_ID_Prop +
                (Year + 0 |Ecoregion_NARS) +
                (1|Ecoregion_NARS:Agency), 
              data = dat_div_NRSA_USGS,
              REML = T) 

summary(m1a)
# Fixed effects:
#             Estimate Std. Error t value
# (Intercept) 18.39323    2.24987   8.175
# Year         0.07703    0.03796   2.029
# Gen_ID_Prop  5.35462    3.01198   1.778
# Slope of year is 0.07703, so for every year we add 0.08 genera

car::Anova(m1a, test = "F", type = "II") 
# Analysis of Deviance Table (Type II Wald F tests with Kenward-Roger df)
# 
# Response: alpha.bar
#                  F Df  Df.res  Pr(>F)  
# Year        4.2240  1  23.666 0.05105 .
# Gen_ID_Prop 3.1069  1 262.098 0.07913 .

performance::r2(m1a)
# R2 for Mixed Models
# 
#   Conditional R2: 0.715
#      Marginal R2: 0.033


dat1 <- data.frame(ggpredict(m1a, terms = "Year [1:27]", type = "fixed"))
dat2 <- data.frame(ggpredict(m1a, terms = c("Year [1:27]", "Ecoregion_NARS","Agency"), 
                             type = "random"))

# filter dat2, so effects on span the length of which we have data
min_max <- dat_div_NRSA_USGS %>%
  group_by(Ecoregion_NARS, Agency) %>%
  summarize(min = min(Year),
            max = max(Year))

dat2 <- dat2 %>% 
  left_join(min_max, by = c("group" = "Ecoregion_NARS",
                            "facet" = "Agency")) %>%
  filter(x >= min,
         x <= max) 

m1a_fig <- ggplot() +
  geom_line(data=dat1, 
            aes(x = x + 1992, y = predicted), 
            color = "black", size=2) +
  geom_ribbon(data=dat1, 
              aes(x=x + 1992, ymin = conf.low, ymax = conf.high),
              alpha=0.3) +
  geom_line(data = dat2, 
            aes(x = x + 1992, y = predicted, color = group),
            linetype = "dashed", size = 1) +
  facet_wrap(~facet) +
  scale_color_brewer(palette="Set1", name = "Ecoregion") +
  scale_x_continuous(breaks=c(1993,1998,2003,2008,2013, 2018))+
  scale_y_continuous(breaks=c(15,20,25,30,35,40), limits = c(15,40))+
  ylab(expression(atop(""*alpha*" diversity","(mean richness of genera across sites)"))) +
  xlab("Year") +
  theme_bw() +
  theme(text = element_text(size = 18)) 

m1a_fig

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
#Fixed effects:
#            Estimate Std. Error t value
#(Intercept)  31.6581     9.2402   3.426
#Year          0.1550     0.1886   0.822
#Gen_ID_Prop  45.9899    11.9841   3.838
# slope of year is 0.1550, so for each year we add 0.16 genera


car::Anova(m1g, test = "F", type = "II")
# Analysis of Deviance Table (Type II Wald F tests with Kenward-Roger df)
# 
# Response: sc.gamma
#                   F Df  Df.res    Pr(>F)    
# Year         0.6794  1  13.581 0.4240191    
# Gen_ID_Prop 14.5813  1 180.823 0.0001843 ***

performance::r2(m1g)
# R2 for Mixed Models
# 
#   Conditional R2: 0.821
#      Marginal R2: 0.045

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

# glmmTMB beta model produces conditional R2 > 1.
# mu of 1.7 is too close to zero, estimate of random effect variances may be
# unreliable. Model's distribution-specific variance is negative. Results are not
# reliable. Opting for LMER
m1b <- lmer(prop.beta.sc ~ Year + Gen_ID_Prop +
                (Year + 0 |Ecoregion_NARS) +
                (1|Ecoregion_NARS:Agency), 
              weights = wt.std,
              data = dat_div_NRSA_USGS,
             REML = TRUE)   

car::Anova(m1b, test = "F",type = "II")
# Analysis of Deviance Table (Type II Wald F tests with Kenward-Roger df)
# 
# Response: prop.beta.sc
#                   F Df  Df.res   Pr(>F)    
# Year         0.9717  1  12.118 0.343534    
# Gen_ID_Prop 14.2303  1 183.874 0.000218 ***

performance::r2(m1b)
# R2 for Mixed Models
# 
# Conditional R2: 0.719
#    Marginal R2: 0.029

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

ggsave("./figures_maintext/div_time_20Jan.png", f2_legend, dpi = 200, 
       height = 10.5, width = 16)


#### Supplement Figure 7
A = m1_fig + theme(legend.position = "none")
B = m3_fig_b + theme(legend.position = "none")
C = m1a_fig + theme(legend.position = "none")

legend = cowplot::get_legend(m1_fig + 
                               theme(legend.position = "top",
                                     text = element_text(size = 25)) +
                               guides(colour = guide_legend(nrow = 2)))

f2 <- cowplot::plot_grid(A, B, C, labels = c("A", "B", "C"),
                         label_size = 20,
                         ncol = 1,
                         nrow = 3,
                         align = "v")
f2_legend <- cowplot::plot_grid(legend, f2, ncol = 1, rel_heights = c(1,15))

ggsave("./figures_supplement/div_time_site_dens_alpha_20Jan.png", f2_legend, dpi = 200, 
       height = 15, width = 9)


#### code for plotting predicted to observed for supplement

dat_abun$predicted = (predict(m2, dat_abun, type = "response"))

dat_alpha_all_Sn$predicted_alpha = predict(m_Sn, dat_alpha_all_Sn, type = "response")

dat_div_NRSA_USGS$predicted_gamma = predict(m1g, dat_div_NRSA_USGS, type = "response")

dat_div_NRSA_USGS$predicted_beta = predict(m1b, dat_div_NRSA_USGS, type = "response")


predobs <- data.frame(Predicted = c(dat_abun$predicted,
                                    dat_alpha_all_Sn$predicted_alpha,
                                    dat_div_NRSA_USGS$predicted_gamma,
                                    dat_div_NRSA_USGS$predicted_beta),
                      Observed = c(log(dat_abun$tot_abun),
                                   dat_alpha_all_Sn$alpha,
                                   dat_div_NRSA_USGS$sc.gamma,
                                   dat_div_NRSA_USGS$prop.beta.sc),
                      Endpoint = c(rep("Density", times = nrow(dat_abun)),
                                   rep("Alpha", times = nrow(dat_alpha_all_Sn)),
                                   rep(c("Gamma", "Beta"), each = nrow(dat_div_NRSA_USGS))))

predobs$Endpoint <- factor(predobs$Endpoint, levels = c("Density", "Alpha",
                                                        "Gamma", "Beta"))

theme_nmds <- function(){
  theme_set(theme_classic())+
    theme(
      panel.background=element_rect(color="black",
                                    fill="white"),
      legend.background=element_blank(),
      axis.text = element_text(size=16,color="black"),
      axis.title.y=element_text(vjust=0.42),
      panel.grid.major = element_blank(),
      panel.grid.minor =element_blank(),
      legend.title = element_text(size=18),
      legend.text = element_text(size=16),
      legend.key = element_rect(colour = NA, fill = NA),
      strip.text.x = element_text(size=14),
      strip.background = element_rect(color="black"),
      axis.title=element_text(size=15)
    )
}

ggplot(predobs, aes(x = Predicted, y = Observed))+
  facet_wrap(~Endpoint, scale = "free")+
  geom_point(color = "dark grey")+
  geom_abline(intercept = 0, slope = 1, size = 1)+
  labs(y = "Observed values",
       x = "Predicted values")+
  theme_nmds()+
  theme(strip.background = element_rect(color = "black"),
        axis.text = element_text(size = 12))

ggsave("./figures_supplement/PredVsObsOvrl.jpg", dpi = 300, height = 6.2, width = 5.75)
