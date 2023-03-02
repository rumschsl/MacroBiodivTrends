# Objective: Diversity analyses according to landuse and including random effect 
# of agency with linear models

library(tidyverse)
library(lme4)
library(ggeffects)
library(emmeans)
library(glmmTMB)

# read in density data
dat_abun <- read.csv("./data/dat_abun_all.csv", row.names = 1, header =T) 

# read in diversity data
dat_div_all_LU_Ag <- read.csv("./data/dat_div_all_LU_Ag_4Jan.csv", 
                              header = T, 
                              row.names = 1) %>%
  # drop ER_YR_LU_AG combos not estimated within 0.7-0.8 sample coverage
  filter(sc >= 0.7 & sc <= 0.8)

###         ###
### DENSITY ###
###         ###

# this model fits
m2 <- lmer(log(tot_abun) ~ YearCont + landuse  + Ecoregion_NARS +
             PropID + AreaSampTot_m2 + Agency + Gen_ID_Prop +
             Ecoregion_NARS:Agency +
             YearCont:landuse + 
             YearCont:Ecoregion_NARS +
             (1|SiteNumber),
           data = dat_abun,
           REML = TRUE)

# calculate magnitude of trend of YearCont
slm2 <- emtrends(m2, specs =  ~ landuse, var = "YearCont", pbkrtest.limit = 8236,
                lmerTest.limit = 8236)
slm2
# landuse        YearCont.trend      SE   df lower.CL  upper.CL
# Agriculture          -0.00282 0.00261 6950 -0.00794  0.002302
# Forest_Wetland       -0.00419 0.00205 7319 -0.00821 -0.000182
# Other                -0.00661 0.00295 7245 -0.01238 -0.000832
# Urban                -0.00217 0.00322 7737 -0.00848  0.004137

(exp(c(-0.00282, -0.00419, -0.00661, -0.00217)) - 1)*100
# [1] -0.2816028 -0.4181234 -0.6588202 -0.2167647

performance::r2(m2)
# R2 for Mixed Models
# 
# Conditional R2: 0.813
#    Marginal R2: 0.772

summary(m2)

car::Anova(m2, type = "III", test = "F") 
#Analysis of Deviance Table (Type III Wald F tests with Kenward-Roger df)
#
#Response: log(tot_abun)
#                                  F Df Df.res    Pr(>F)    
# (Intercept)              5329.4136  1 7478.0 < 2.2e-16 ***
# YearCont                    0.6556  1 7524.7  0.418135    
# landuse                     0.4076  3 6529.3  0.747516    
# Ecoregion_NARS              2.6775  8 7258.7  0.006163 ** 
# PropID                  18654.0234  1 7925.7 < 2.2e-16 ***
# AreaSampTot_m2           2156.0928  1 7160.8 < 2.2e-16 ***
# Agency                      1.8229  1 4643.4  0.177033    
# Gen_ID_Prop              1298.9723  1 8025.4 < 2.2e-16 ***
# Ecoregion_NARS:Agency      14.5824  8 4078.3 < 2.2e-16 ***
# YearCont:landuse            0.4499  3 7160.1  0.717379    
# YearCont:Ecoregion_NARS     0.9615  8 7491.0  0.464201     

# calculate difference in density by land use with time held at the mean
emm2 <- emmeans(m2, "landuse", 
              pbkrtest.limit = 8236,
              lmerTest.limit = 8236)
emm2
# landuse        emmean      SE   df lower.CL upper.CL
# Agriculture     6.896 0.02066 4264    6.855    6.936
# Forest_Wetland  6.886 0.01606 4330    6.854    6.917
# Other           6.878 0.02253 4156    6.833    6.922
# Urban           6.866 0.02790 3428    6.811    6.921

exp(c(6.896, 6.886, 6.878,6.866))
#[1] 988.3135 978.4797 970.6831 959.1045

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
dat_alpha_all <- read.csv("./data/dat_alpha_all_4Jan.csv", row.names=1, header=T)


m_LU_a <- glmmTMB::glmmTMB(alpha ~ Year + landuse  + Ecoregion_NARS +
             PropID + AreaSampTot_m2 + Agency + Gen_ID_Prop +
             Ecoregion_NARS:Agency +
             Year:landuse + 
             Year:Ecoregion_NARS +
             (1|SiteNumber),
             family = nbinom1,
           data = dat_alpha_all)

# calculate magnitude of trend of YearCont
slm2 <- emtrends(m_LU_a, specs =  ~ landuse, var = "Year", pbkrtest.limit = 8236,
                 lmerTest.limit = 8236)
slm2
#  landuse        Year.trend       SE   df lower.CL upper.CL
# Agriculture       0.00577 0.001274 7630  0.00328  0.00827
# Forest_Wetland    0.00523 0.000938 7630  0.00339  0.00707
# Other             0.00128 0.001350 7630 -0.00137  0.00393
# Urban             0.00706 0.001495 7630  0.00413  0.00999

(exp(c(0.00577, 0.00523, 0.00128, 0.00706)) - 1)*100
# [1] 0.5786679 0.5243700 0.1280820 0.7084981

# LU difference with time held at mean
emm2a <- emmeans(m_LU_a, "landuse") 
emm2a
#landuse        emmean       SE   df lower.CL upper.CL
#Agriculture     3.104 0.010173 7630    3.084    3.124
#Forest_Wetland  3.187 0.007671 7630    3.172    3.202
#Other           3.135 0.010941 7630    3.113    3.156
#Urban           3.001 0.013414 7630    2.975    3.027

(exp(c(3.104, 3.187, 3.135, 3.001)) - 1)
#[1] 21.28692(ag) 23.21567(forest) 21.98864(grass/shrub)  19.10563(Urban)

pairs(emm2a, adjust = "none")
# contrast                     estimate     SE   df t.ratio p.value
# Agriculture - Forest_Wetland  -0.0832 0.0121 7630  -6.899  <.0001
# Agriculture - Other           -0.0313 0.0152 7630  -2.064  0.0391
# Agriculture - Urban            0.1026 0.0162 7630   6.330  <.0001 
# Forest_Wetland - Other         0.0520 0.0134 7630   3.883  0.0001 
# Forest_Wetland - Urban         0.1859 0.0147 7630  12.658  <.0001 
# Other - Urban                  0.1339 0.0175 7630   7.657  <.0001 

performance::r2(m_LU_a)
# R2 for Mixed Models
# 
# Conditional R2: 0.634
#    Marginal R2: 0.358

summary(m_LU_a)

car::Anova(m_LU_a, type = "III") 
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
#Response: alpha
#                           Chisq Df Pr(>Chisq)    
# (Intercept)           2706.5031  1  < 2.2e-16 ***
# Year                     6.4980  1    0.01080 *  
# landuse                 65.6872  3  3.576e-14 ***
# Ecoregion_NARS         162.5553  8  < 2.2e-16 ***
# PropID                   3.7961  1    0.05137 .  
# AreaSampTot_m2          39.2427  1  3.743e-10 ***
# Agency                  19.4696  1  1.022e-05 ***
# Gen_ID_Prop            114.6871  1  < 2.2e-16 ***
# Ecoregion_NARS:Agency   46.2078  8  2.170e-07 ***
# Year:landuse             9.4261  3    0.02413 *  
# Year:Ecoregion_NARS     32.9669  8  6.244e-05 ***  

# Year*LU significant
pairs(emtrends(m_LU_a, "landuse", var = "Year"), adjust = "none")
# contrast                      estimate      SE   df t.ratio p.value
# Agriculture - Forest_Wetland  0.000544 0.00149 7630   0.366  0.7142
# Agriculture - Other           0.004492 0.00188 7630   2.386  0.0171
# Agriculture - Urban          -0.001286 0.00186 7630  -0.690  0.4904
# Forest_Wetland - Other        0.003947 0.00170 7630   2.316  0.0206
# Forest_Wetland - Urban       -0.001830 0.00162 7630  -1.131  0.2579
# Other - Urban                -0.005777 0.00199 7630  -2.905  0.0037

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
# landuse        Year.trend    SE  df lower.CL upper.CL
# Agriculture         0.244 0.150 464  -0.0494   0.5383
# Forest_Wetland      0.362 0.129 596   0.1072   0.6158
# Other               0.532 0.186 412   0.1666   0.8982
# Urban              -0.228 0.162 379  -0.5471   0.0915

# Within ecoregions, urban lost 0.23 genera per year (6.21 genera over 27 years), 
# respectively, while agricultural, forest/wetland, and grassland/shrub streams 
# gained 0.24, 0.36 and 0.53 genera per year (6.48, 9.72, and 14.31 genera over 
# 27 years), respectively. 

performance::r2(m2g_LU) 
# R2 for Mixed Models
# 
# Conditional R2: 0.739
#    Marginal R2: 0.679

summary(m2g_LU)

car::Anova(m2g_LU, test = "F", type = "III")
# Analysis of Deviance Table (Type III Wald F tests with Kenward-Roger df)
# 
# Response: sc.gamma
#                    F Df Df.res    Pr(>F)    
# (Intercept)  78.8069  1  63.31 1.014e-12 ***
# Year          2.6733  1 463.93   0.10272    
# landuse       1.2543  3  52.91   0.29949    
# Agency       96.0296  1   8.31 7.534e-06 *** 
# Gen_ID_Prop   4.3590  1 357.26   0.03752 *   
# Year:landuse  4.7157  3 451.66   0.00298 **  

pairs(emtrends(m2g_LU, "landuse", var = "Year"), adjust = "none")
# contrast                     estimate    SE  df t.ratio p.value
# Agriculture - Forest_Wetland   -0.117 0.175 544  -0.671  0.5026
# Agriculture - Other            -0.288 0.223 428  -1.292  0.1969
# Agriculture - Urban             0.472 0.198 426   2.385  0.0175 #
# Forest_Wetland - Other         -0.171 0.212 440  -0.806  0.4205
# Forest_Wetland - Urban          0.589 0.182 465   3.240  0.0013 #
# Other - Urban                   0.760 0.233 386   3.269  0.0012 #

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
# Conditional R2: 0.429
#    Marginal R2: 0.282

car::Anova(m2b_LU, test = "F", type = "III")
# Analysis of Deviance Table (Type III Wald F tests with Kenward-Roger df)
# 
# Response: prop.beta.sc
#                    F Df Df.res    Pr(>F)    
# (Intercept)  61.7624  1 212.06 1.925e-13 ***
# Year          0.0234  1 464.44    0.8784    
# landuse       0.5424  3  66.07    0.6549    
# Agency       54.7328  1   8.62 5.182e-05 ***
# Gen_ID_Prop   0.0018  1 388.42    0.9665    
# Year:landuse  2.0188  3 454.86    0.1105    

pairs(emtrends(m2b_LU, "landuse", var = "Year"), adjust = "none")
# contrast                      estimate      SE  df t.ratio p.value
# Agriculture - Forest_Wetland -0.000425 0.00146 553  -0.292  0.7706 
# Agriculture - Other           0.001085 0.00187 429   0.581  0.5615
# Agriculture - Urban          -0.003206 0.00165 431  -1.937  0.0533 
# Forest_Wetland - Other        0.001510 0.00177 446   0.853  0.3944
# Forest_Wetland - Urban       -0.002781 0.00152 466  -1.827  0.0683 
# Other - Urban                -0.004291 0.00195 382  -2.201  0.0283 #

slm2b <- emtrends(m2b_LU, specs =  ~ landuse, var = "Year")
slm2b
# landuse        Year.trend      SE  df  lower.CL upper.CL
# Agriculture     -0.000192 0.00126 464 -0.002660  0.00228
# Forest_Wetland   0.000233 0.00108 613 -0.001893  0.00236
# Other           -0.001277 0.00155 421 -0.004330  0.00178
# Urban            0.003014 0.00136 380  0.000341  0.00569

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

ggsave("./figures_maintext/div_time_landuse_26Jan.png", f2_legend, dpi = 200, 
       height = 10.5, width = 16)

#### code for plotting predicted to observed for supplement
dat_abun$predicted2 = predict(m2, dat_abun, type = "response")

# error is coming from 95 sites lacking PropID or TotAreaSamp
dat_alpha_all$predicted_alpha2 = predict(m_LU_a, dat_alpha_all, type = "response")

dat_div_all_LU_Ag$predicted_gamma2 = predict(m2g_LU, dat_div_all_LU_Ag, type = "response")

dat_div_all_LU_Ag$predicted_beta2 = predict(m2b_LU, dat_div_all_LU_Ag, type = "response")

predobs2 <- data.frame(Predicted = c(dat_abun$predicted2,
                                     dat_alpha_all$predicted_alpha2,
                                     dat_div_all_LU_Ag$predicted_gamma2,
                                     dat_div_all_LU_Ag$predicted_beta2),
                       Observed = c(log(dat_abun$tot_abun),
                                    dat_alpha_all$alpha,
                                    dat_div_all_LU_Ag$sc.gamma,
                                    dat_div_all_LU_Ag$prop.beta.sc),
                       Endpoint = c(rep("Density", times = nrow(dat_abun)),
                                    rep("Alpha", times = nrow(dat_alpha_all)),
                                    rep(c("Gamma", "Beta"), each = nrow(dat_div_all_LU_Ag))),
                       Landuse = c(dat_abun$landuse,
                                   dat_alpha_all$landuse,
                                   dat_div_all_LU_Ag$landuse,
                                   dat_div_all_LU_Ag$landuse))

predobs2$Endpoint <- factor(predobs2$Endpoint, levels = c("Density", "Alpha",
                                                          "Gamma", "Beta"))

LU_names <- c(
  Agriculture = "Agriculture",
  Forest_Wetland = "Forest/Wetland",
  Other = "Grassland/Shrub",
  Urban = "Urban"
)

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

ggplot(predobs2, aes(y = Observed, x = Predicted))+
  facet_wrap(Landuse~Endpoint, scale = "free",
             labeller = labeller(Landuse = LU_names))+
  # geom_point(color = "dark grey")+
  geom_point(aes(color = Landuse))+
  geom_abline(intercept = 0, slope = 1, size = 1)+
  labs(x = "Predicted values",
       y = "Observed values")+
  theme_nmds()+
  scale_color_brewer(palette="Dark2", name = "Dominant landuse",
                     limits = c("Forest_Wetland","Agriculture","Urban","Other"),
                     breaks = c("Forest_Wetland","Other","Agriculture","Urban"),
                     labels = c("Forest/wetland","Grassland/shrub","Agriculture","Urban")) +
  theme(strip.background = element_rect(color = "black"),
        strip.text.x = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.position = "none")


ggsave("./figures_supplement/PredVsObsLanduse.jpg", dpi = 300, height = 7.75, width = 6)