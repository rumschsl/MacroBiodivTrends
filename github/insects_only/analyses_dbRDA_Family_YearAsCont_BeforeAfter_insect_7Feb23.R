# Objective: composition analyses through time according to land use for insects

# load packages
library(tidyverse)
library(vegan)

# read in datasets

##Read in the dataset
insectsfam <- read.csv("./insects_only/data/insectfam.csv")

##For the community matrix; remove the first 5 identifying (non-genera) columns
insectscommF <- insectsfam[,-c(1:5)]

##Remove all genera that appear in fewer than 3 ER-Agency-LU-Year combinations
## singleton and doubleton removal
insectscommF2 <- insectscommF[, which((colSums(insectscommF > 0) > 2))]

invertgenmain <- read.csv("./data/invertgenmain.csv", row.names= 1, header = T)


# functions (from Gavin Simpson's ggvegan package) 
## reading in anything
fortify.cca <- function(model, data, axes = 1:6,
                        display = c("sp", "wa", "lc", "bp", "cn"), ...) {
  ## extract scores
  scrs <- scores(model, choices = axes, display = display, ...)
  ## handle case of only 1 set of scores
  if (length(display) == 1L) {
    scrs <- list(scrs)
    nam <- switch(display,
                  sp = "species",
                  species = "species",
                  wa = "sites",
                  sites = "sites",
                  lc = "constraints",
                  bp = "biplot",
                  cn = "centroids",
                  stop("Unknown value for 'display'"))
    names(scrs) <- nam
  }
  miss <- vapply(scrs, function(x ) all(is.na(x)), logical(1L))
  scrs <- scrs[!miss]
  nams <- names(scrs)
  nr <- vapply(scrs, FUN = NROW, FUN.VALUE = integer(1))
  df <- do.call('rbind', scrs)
  rownames(df) <- NULL
  df <- as.data.frame(df)
  df <- cbind(Score = factor(rep(nams, times = nr)),
              Label = unlist(lapply(scrs, rownames), use.names = FALSE),
              df)
  df
}

arrowMul <- function(arrows, data, at = c(0, 0), fill = 0.75) {
  u <- c(range(data[,1], range(data[,2])))
  u <- u - rep(at, each = 2)
  r <- c(range(arrows[, 1], na.rm = TRUE), range(arrows[, 2], na.rm = TRUE))
  rev <- sign(diff(u))[-2]
  if (rev[1] < 0)
    u[1:2] <- u[2:1]
  if (rev[2] < 0)
    u[3:4] <- u[4:3]
  u <- u/r
  u <- u[is.finite(u) & u > 0]
  fill * min(u)
}








##Generate objects that represent landuse, ecoregion, agency, year, and LU-Year
## combinations; to be used in the dbRDA model itself
Ecoregion <- insectsfam$Ecoregion_NARS
NLCD <- insectsfam$landuse
Year <- insectsfam$Year
Agency <- insectsfam$Agency
PRE_POST <- insectsfam$PRE_POST

##Generate continuous variables of time for each landuse category (0 if not in LU)
YearAg <- ifelse(NLCD == "Agriculture",
                 Year,
                 0)
YearFW <- ifelse(NLCD == "Forest_Wetland",
                 Year,
                 0)
YearOt <- ifelse(NLCD == "Other",
                 Year,
                 0)
YearUb <- ifelse(NLCD == "Urban",
                 Year,
                 0)

##Bray-curtis distances to be used as our response
response = vegdist(insectscommF2, dist = "bray")

##the Condition() means that this is a partial dbRDA; such that the influence of
##Agency and Ecoregion and PRE_POST accounted for, before we constrain
##the ordination by the LU-Year combinations

set.seed(1)
##Run this model to test for significance
modA <- dbrda(response ~ NLCD * Year +
                 Condition(Agency + Ecoregion + PRE_POST))

set.seed(1)
##Run perumtational anova on the reduced model (ignoring conditional terms)
##for each term in the model; this is the recommended test in Numerical Ecology
##and the vegan help docs
anova(modA, permutations = how(nperm = 9999), model = "reduced",
      by = "terms")

# Permutation test for dbrda under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
#
# Model: dbrda(formula = response ~ NLCD * Year + Condition(Agency + Ecoregion + PRE_POST))
#            Df SumOfSqs       F Pr(>F)    
# NLCD        3    5.397 27.3214 0.0001 ***
# Year        1    0.115  1.7515 0.0984 .  
# NLCD:Year   3    0.446  2.2568 0.0039 ** 
# Residual  516   33.977                   
                                   
##Run this model to generate biplot for each NLCD-Year combination
##Need to run 2 of these to get the biplot for each of the land use-year combinations
set.seed(1)
modF <- dbrda(response ~ NLCD + Year + YearFW + YearAg + YearUb + YearOt + 
                       Condition(Agency + Ecoregion + PRE_POST))
set.seed(1)
modF2 <- dbrda(response ~ NLCD + Year + YearOt + YearFW + YearAg + YearUb + 
                 Condition(Agency + Ecoregion+ PRE_POST))


##Model summaries are identical, just easier to plot output from modF
summary(modA)$concont$importance[,c(1:5)]
summary(modF)$concont$importance[,c(1:5)]

modA
##Explains 41.89% of variance ## 1-0.5811

summary(modA)$concont$importance[,c(1:5)]
##dbRDA1 explains 84.1% of fitted variance
##dbRDA2 explains 8.1% of fitted variance

summary(modA)$cont$importance[,c(1:5)]
##dbRDA1 explains 12.5% of total variance
##dbRDA2 explains 1.2% of total variance


##Extract biplot vectors for each of the year terms (main plus NLCDs); drop
##first three rows which are the biplots for the main effects of landuse
##Fortify - this is all from Gavin Simpson's autoplot.rda function
obj <- fortify.cca(modF, axes = c(1:2))

##fortify the second figure rda (for the "other" vector that isn't generated)
obj2 <- fortify.cca(modF2, axes = c(1,2))
obj[(nrow(obj) + 1), ] <- obj2[which(obj2$Label == "YearOt"),]

## sort out x, y aesthetics
vars <- colnames(obj)[3:4]

## subset out the layers wanted
obj <- obj[obj[["Score"]] %in% c("sites", "biplot", "centroids"), , drop = FALSE]

## remove biplot arrows for centroids if present
want <- obj[["Score"]] == "biplot"
tmp <- obj[want, ]
obj <- obj[!want, ]
bnam <- tmp[, "Label"]
cnam <- obj[obj[["Score"]] == "centroids", "Label"]
obj <- rbind(obj, tmp[!bnam %in% cnam, , drop = FALSE])

## scale biplot vectors to length 
want <- obj[["Score"]] == "biplot"
## use 'mul' for adjusting vector length of the species!!
mul <- arrowMul(obj[want, vars, drop = FALSE],
                obj[!want, vars, drop = FALSE])

obj[want,vars] <- mul * obj[want, vars]

##Generate the biplot vectors (bpvecs) 
bpvecs <- obj[want,]
bpvecs$Label2 <- c("Year", "Year*Forest/wetland",
                   "Year*Agriculture", "Year*Urban", 
                   "Year*Grassland/shrub")
bpvecs$Label3 <- c("Year", "Year*Forest/wetland",
                   "Year*Agriculture", "Year*Urban", 
                   "Year*Grassland/shrub")


##Generate the centroids (centr) 
wantc <- obj[["Score"]] == "centroids"
centr <- obj[wantc, , drop = FALSE]
centr$Label = substr(centr$Label, 5, nchar(centr$Label))

##Generate the site scores dataset (with identifying info)
sitescores <- obj %>% filter(Score == "sites")
sitescores$Landuse <- NLCD
sitescores$Year <- Year

##Generate the 95% confidence interval ellipses for the centroid of each landuse category
#First, run the ordiellipse function to get the estimates
ord <- ordiellipse(modF, NLCD, kind = "se", conf = 0.95, draw = "none", label = F)

#function to generate ellipses from the ordiellipse function
veganCovEllipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100){
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#Generate ellipse points in a single dataset
df_ell <- data.frame()
for(g in unique(NLCD)){
  df_ell <- rbind(df_ell,
                  data.frame(veganCovEllipse(ord[[g]]$cov,
                                             ord[[g]]$center,
                                             ord[[g]]$scale),
                             Landuse = g))
}

landuse_names <- c(
  `Forest_Wetland` = "Forest/wetland",
  `Agriculture` = "Agriculture",
  `Urban` = "Urban",
  `Other` = "Other"
)

# write a ggplot theme for figures
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

#modf flipped axis 1 (can occur, does not alter the interpretation - 
# this is still the same ordination result as modf), need to flip it here to match
bpvecs$dbRDA1[1:4] <- bpvecs$dbRDA1[1:4]*-1

A <- ggplot()+
  # facet_wrap(~Landuse)+
  #scale_y_continuous(limits = c(-2.5, 2.5))+
  #scale_x_continuous(limits = c(-2.5, 2.5))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  ##Note: r is the radius of the circle and needs to be set to the scaling value
  ##used for the predictor vectors (mul); 
  ##multiple the dbRDA1 and dbRDA2 by this value in the geom_segment code
  ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = mul), color = "black")+
  # geom_point(aes(color = Landuse, alpha = Year),size = 2) +
  geom_segment(data = bpvecs %>% filter(Label == "Year"),
               aes(x = 0, y = 0,
                   xend = dbRDA1, yend = dbRDA2),
               arrow = arrow(length = unit(0.3, "cm")),
               size = 1)+
  geom_segment(data = bpvecs %>% filter(Label == "YearAg"),
               aes(x = 0, y = 0,
                   xend = dbRDA1, yend = dbRDA2),
               color = "#D95F02",
               arrow = arrow(length = unit(0.3, "cm")),
               size = 1)+
  geom_segment(data = bpvecs %>% filter(Label == "YearFW"),
               aes(x = 0, y = 0,
                   xend = dbRDA1, yend = dbRDA2),
               color = "#1C9E77",
               arrow = arrow(length = unit(0.3, "cm")),
               size = 1)+
  geom_segment(data = bpvecs %>% filter(Label == "YearUb"),
               aes(x = 0, y = 0,
                   xend = dbRDA1, yend = dbRDA2),
               color = "#7570B3",
               arrow = arrow(length = unit(0.3, "cm")),
               size = 1)+
  geom_segment(data = bpvecs %>% filter(Label == "YearOt"),
               aes(x = 0, y = 0,
                   xend = dbRDA1, yend = dbRDA2),
               color = "#E72D8A",
               arrow = arrow(length = unit(0.3, "cm")),
               size = 1)+
  ##Code to plot the 95%CIs for the centroids
  geom_path(data = df_ell,
            #need to flip dbrda1 (comes form modF)
            aes(x = dbRDA1*-1, y = dbRDA2, color = Landuse), size=1) +
  geom_point(data = centr,
             #need to flip dbrda1 (comes form modF)
             aes(x = dbRDA1*-1, y = dbRDA2, color = Label), size=1) +
  ##Centroids (putting as text, but could be placed as other geom)
  # ggrepel::geom_text_repel(data = centr,
  #                          aes(x = dbRDA1, y = dbRDA2, label = Label),
  #                          size = 3, max.time = 10, box.padding = 0.1, point.padding = 0.1,
  #                          fontface = "bold",
  #                          #alpha = 0.7,
  #                          bg.color = "white",
  #                          bg.r = .1,
  #                          min.segment.length = 0,
  #                          nudge_y = -0.3)+
  ##Vectors
  ggrepel::geom_text_repel(data = bpvecs,
                           aes(x = dbRDA1, y = dbRDA2, label = Label3),
                           size =4,
                           nudge_y = 0.15,
                           max.time = 10, box.padding = 0.1, point.padding = 0.5,
                           bg.color = "white",
                           bg.r = .1)+
  xlab("dbRDA1 (84.1% of fitted, 12.5% of total variation)") +
  ylab("dbRDA2 (8.1% of fitted, 1.2% of total variation)")+
  scale_color_brewer(palette="Dark2", name = "Dominant landuse",
                     limits = c("Forest_Wetland","Agriculture","Urban","Other"),
                     breaks = c("Forest_Wetland","Other","Agriculture","Urban"),
                     labels = c("Forest/wetland","Grassland/shrub","Agriculture","Urban")) +
  theme_nmds()+
  theme(strip.background = element_rect(fill = "light gray", color = "black"),
        legend.position = "bottom") +
  guides(colour = guide_legend(ncol = 1))

A

##Hellinger transformation of the community data
insectscommF3 <- sqrt(decostand(insectscommF2, method = "hellinger"))

##Use envfit() with ALL families as the predictor (really the response) to
##generate species scores (correlation vectors sensu PRIMER and MJ Anderson)
en1F <- envfit(modA ~ . , data = insectscommF3, permutations = 0)

##Convert this to a dataframe
genscoresF <- as.data.frame(scores(en1F, display = "vectors"))

##add taxa names from the rownmaes and reset rownames as numeric
genscoresF <- cbind(genscoresF, Taxa = rownames(genscoresF))
rownames(genscoresF) <- 1:nrow(genscoresF)

# create a column to color family vectors by order for all but chironomids, 
# which will be colored by family

# match the Taxa (families) in genscoresF to Families in insectsgenmain
# and pull the col_label (order)
genscoresF$col_label <- invertgenmain$col_label[match(genscoresF$Taxa, 
                                                       invertgenmain$Family)]

# above code will create NA's for Chironomids since Chironomids are subfamilies in analyses
# so, for those NA's match and pull the subfamily info from insectsgenmain
# and pull the corresponding col_label (subfamily)
genscoresF$col_label <- ifelse(is.na(genscoresF$col_label), 
                               invertgenmain$col_label[match(genscoresF$Taxa, 
                                                             invertgenmain$Subfamily)],
                                genscoresF$col_label)

# calculate vector lengths for each family
genscoresF$length <- sqrt(genscoresF$dbRDA1^2 + genscoresF$dbRDA2^2)

# object used to color groups in figure
colOrder <- readRDS("./data/colOrder.rds")


B <-genscoresF %>% 
  filter(length > 0.3) %>%
  ggplot()+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  #scale_y_continuous(limits = c(-2.5,2.5))+
  #scale_x_continuous(limits = c(-2.5,2.5))+
  ##Note: r is the radius of the circle and needs to be set to the scaling value
  ##used for the predictor vectors (mul); 
  ##multiple the dbRDA1 and dbRDA2 by this value in the geom_segment code
  ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = mul), color = "black")+
  geom_segment(aes(x = 0, xend = dbRDA1*mul,
                   y = 0, yend = dbRDA2*mul,
                   color = col_label, group = Taxa),
               arrow = arrow(length = unit(0.25, "cm")),
               size = 1) +
  colOrder +
  xlab("dbRDA1 (84.1% of fitted, 12.5% of total variation)") +
  ylab("dbRDA2 (8.1% of fitted, 1.2% of total variation)")+
  theme_nmds()+
  guides(colour = guide_legend(ncol = 2))+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=16))
B

A2 = A + theme(legend.position = "none") 
B2 = B + theme(legend.position = "none") 

f2 <- cowplot::plot_grid(A2, B2, labels = c("A", "B"),
                         label_size = 20,
                         nrow = 1,
                         rel_heights = c(1,1),
                         rel_widths = c(1,1))

legend_A = cowplot::get_legend(A) 
legend_B = cowplot::get_legend(B)

f2_legend = cowplot::plot_grid(legend_A, legend_B, 
                               nrow = 1,
                               rel_heights = c(1,1),
                               rel_widths = c(1,1))

f2_all = cowplot::plot_grid(f2, f2_legend, ncol = 1, rel_heights = c(6,2))

#read out insectsgenamin & colOrder
#saveRDS(colOrder, "./data/colOrder.rds")

ggsave("./insects_only/figures/composition.png", f2_all, dpi = 200, 
       height = 8, width = 12)
