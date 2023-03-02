# Samantha Rumschlag
# 31 Jan 2023
# Objective: Diversity analyses including random effect of agency with linear 
# models - without Gen_ID_Prop to evaluate how trends differ without accounting 
# for improvements in taxonomic identification

# load packages
library(tidyverse)
library(lme4)
library(ggeffects)
library(emmeans)

# read in datasets

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

###         ###
### DENSITY ###
###         ###

# create a model
m2 <- lmer(log(tot_abun) ~ YearCont + Ecoregion_NARS + YearCont:Ecoregion_NARS +
             PropID + AreaSampTot_m2 + Agency + #Gen_ID_Prop +
             Ecoregion_NARS:Agency + 
             (1|SiteNumber),
           data = dat_abun,
           REML = TRUE)

sl2 <- emtrends(m2, specs = ~ YearCont, var = 'YearCont', 
                  pbkrtest.limit = 8236,
                  lmerTest.limit = 8236)
summary(sl2)
# YearCont YearCont.trend     SE   df lower.CL upper.CL
#     16.6        0.00832 0.00142 7752  0.00553   0.0111

#YearCont trend is exp(-0.0083) for year 1 year increase
(exp(0.00832)-1) * 100
# [1] 0.8354707
# So, every year abundance increase by 0.84%

(exp(27*0.008354707))-1
# [1] 0.251872
# 25% increases in total density over 27 years

# car::Anova(m2, type = "II", test = "F")  
# Response: log(tot_abun)
#                                  F Df Df.res    Pr(>F)    
# YearCont                   49.8571  1 7609.6 1.799e-12 ***
# Ecoregion_NARS             12.4555  8 4463.0 < 2.2e-16 ***
# PropID                  16057.3648  1 7983.7 < 2.2e-16 ***
# AreaSampTot_m2           1918.3120  1 7443.4 < 2.2e-16 ***
# Agency                     31.8210  1 4572.1 1.792e-08 ***
# YearCont:Ecoregion_NARS     1.8468  8 7644.1   0.06386 .  
# Ecoregion_NARS:Agency      23.9774  8 4312.5 < 2.2e-16 ***

summary(m2)

performance::r2(m2)
# R2 for Mixed Models
# 
# Conditional R2: 0.788
#    Marginal R2: 0.734

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
dat_alpha_all_Sn <- read.csv("./data/dat_alpha_all.csv", row.names=1, header=T)

m_Sn <- glmmTMB::glmmTMB(alpha ~ Year + Ecoregion_NARS + Year:Ecoregion_NARS +
                           PropID + AreaSampTot_m2 + Agency + #Gen_ID_Prop +
                           Ecoregion_NARS:Agency + 
                           (1|SiteNumber),
                         family = poisson,
                         data = dat_alpha_all_Sn)

sl1 <- emtrends(m_Sn, specs = ~ Year, var = 'Year',
                pbkrtest.limit = 8240,
                lmerTest.limit = 8240)
summary(sl1)
#Year Year.trend       SE   df lower.CL upper.CL
#15.3    0.00614 0.000608 7640  0.00495  0.00733

car::Anova(m_Sn, type = "II") 
#Analysis of Deviance Table (Type II Wald chisquare tests)
#
# Response: alpha
#                         Chisq Df Pr(>Chisq)    
# Year                  146.633  1  < 2.2e-16 ***
# Ecoregion_NARS        953.135  8  < 2.2e-16 ***
# PropID                  8.078  1  0.0044805 ** 
# AreaSampTot_m2         35.090  1  3.149e-09 ***
# Agency                848.761  1  < 2.2e-16 ***
# Year:Ecoregion_NARS    31.391  8  0.0001197 ***
# Ecoregion_NARS:Agency  75.494  8  3.927e-13 ***

performance::r2(m_Sn)
# R2 for Mixed Models
# 
# Conditional R2: 0.685
#    Marginal R2: 0.328

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

m1g <- lmer(sc.gamma ~ Year + #Gen_ID_Prop +
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
#(Intercept)  61.3046     5.3698  11.417
#Year          0.5279     0.1718   3.072

# slope of year is 0.5279, so for each year we add 0.53 genera

car::Anova(m1g, test = "F", type = "II")
# Analysis of Deviance Table (Type II Wald F tests with Kenward-Roger df)
# 
# Response: sc.gamma
#          F Df Df.res  Pr(>F)  
# Year 9.545  1  7.697 0.01562 *

performance::r2(m1g)
# R2 for Mixed Models
# 
#   Conditional R2: 0.828
#      Marginal R2: 0.024

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
# reliable.
m1b <- lmer(prop.beta.sc ~ Year + #Gen_ID_Prop +
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
# Year         0.5077  1  8.1368   0.4961    

performance::r2(m1b)
# R2 for Mixed Models
# 
# Conditional R2: 0.718
#    Marginal R2: 0.003

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

ggsave("./figures_maintext/div_time_no_Gen_ID_Prop.png", f2_legend, dpi = 200, 
       height = 10.5, width = 16)
