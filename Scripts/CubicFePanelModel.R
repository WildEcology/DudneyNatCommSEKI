## ---------------------------
##
## Script name: Cubic model results
##
## Author: Dr. Joan Dudney
##
## Date Created: 2021-05-02
##
## Copyright (c) Joan Dudney, 2021
## Email: jdudney@berkeley.edu
##
## ---------------------------
##
## Notes:This code runs the FE panel
##       model with a cubic term
##       and estimatese p.p. change 
##       and range shift
##       
##     Includes:
##      1. Cubic FE panel model
##      2. MC simulation
##      3. Estimated range shift
## ---------------------------

library(tidyverse)
require(foreign)
require(nnet)
library(ggeffects)
library(magrittr)
library(margins)
library(ggthemes)
library(modelr)
library(sjmisc)
library(sjlabelled)
library(fixest)
library(plotrix)
library(egg)
library(ggpmisc)
library(mvtnorm)
library(clubSandwich)

select=dplyr::select
rename=dplyr::rename
group_by=dplyr::group_by


theme_set(
  theme_bw(base_size = 11)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5))
)


#DATA CLEANING
pathogen_dat=read_csv("Data/survey1and2_wpbr.csv")

panel00=pathogen_dat%>% #6559 live trees
  filter(status==1)%>%
  group_by(plot, vpd00, density00)%>%
  dplyr::summarize(perinc00=sum(inc_tot, na.rm=T)/sum(status, na.rm=T),
                   dbh00=mean(dbh, na.rm=T))%>%
  ungroup()%>%
  mutate(year=2016)%>%
  set_colnames(c( "plot", "vpd", "density", "perinc",    
                  "dbh", "year"))%>%
  mutate(year=factor(year), plot=factor(plot))

panel95=pathogen_dat%>% #7031 live trees
  filter(era95_status==1)%>%
  group_by(plot, vpd95, density95)%>%
  dplyr::summarize(perinc95=sum(era95_inc_tot, na.rm=T)/sum(era95_status, na.rm=T),
                   dbh95=mean(era95_dbh, na.rm=T))%>%
  ungroup()%>%
  mutate(year=1995)%>%
  set_colnames(c( "plot", "vpd", "density", "perinc",    
                  "dbh", "year"))%>%
  mutate(year=factor(year), plot=factor(plot))

all_pan_dat=panel95%>%
  full_join(panel00)%>%
  mutate(year=factor(year), plot=factor(plot))


#======================================================================================================
##CUBIC FIXED EFFECTS PANEL MODEL
##======================================================================================================

##PANEL MODEL USING FEOLS
fe_mod = feols(perinc~ vpd+I(vpd^2)+I(vpd^3)+dbh+density | plot + year,
               data = na.omit(all_pan_dat))

summary(fe_mod)


V_CR1 <- vcovCR(fe_mod, cluster = c(all_pan_dat$plot), type = "CR1")
V_CR1=as.matrix(V_CR1)

##randomly sample from coefficient estimates and variance-covariance matrix

coef_vector = fe_mod$coefficients

draw = rmvnorm(n = 100000, mean = coef_vector, sigma = V_CR1)


##replacing data00 vpd with vpd from 1995
vpd95=select(panel95, plot, vpd)

dat95=panel00%>%
  ungroup()%>%
  select(-vpd)%>%
  left_join(vpd95)%>%
  select(plot, vpd ,density, perinc , dbh ,year)

##getting elevation data
elevation=pathogen_dat%>%
  select(plot, elevation)%>%
  distinct()%>%
  na.omit()%>%
  mutate(plot=as.factor(plot))

dat=select(panel00,plot)


#======================================================================================================
##CALCULATING P.POINT DIFFERENCE ACROSS ELEVATION TERCILES
##======================================================================================================


##DATA ELEVATION TERCILES
datelev_new=select(panel00,plot)%>%
  left_join(elevation)%>%
  mutate(elev_third=ifelse(elevation<=2655.204, "Low (1387-2655)",
                           ifelse(elevation>2655.204&elevation<=3132.137, "Mid (2656-3132)",
                                  ifelse(elevation>=3132.137, "High (3133-3486)", "NA"))))
df_ranges=data.frame()

for (i in 1:10000){
  
 
  ##now run the monte carlo simulation
  d = draw[i,]
  modified_fe_mod = fe_mod
  modified_fe_mod$coefficients = d
  datelev_new$vals_00 = predict(modified_fe_mod, newdata = panel00)
  datelev_new$vals_95 = predict(modified_fe_mod, newdata = dat95)
  
  dat_range=datelev_new%>%
    mutate(diff00=vals_00-vals_95)%>%
    select(-c(vals_00, vals_95, elevation))%>%
    pivot_longer(-c(elev_third, plot))%>%
    group_by(elev_third, name)%>%
    summarize(mean=mean(value))
  
  tot=data.frame(dat_range)
  df_ranges=rbind(df_ranges, tot)
}


##calculating mean and 95%CI
newranges=df_ranges%>%
  group_by(elev_third, name)%>%
  summarize(allmeans=mean(mean),
            q25=quantile(mean, (.025)),
            q97=quantile(mean, (.975)))

##plot
ggplot(newranges, aes(elev_third, allmeans*100, fill=name))+
  geom_hline(yintercept = 0, linetype="dashed", color="grey")+
  geom_bar(stat="identity",position="dodge", width = 0.5, alpha=.7)+
  geom_errorbar(aes(ymin = q25*100, ymax = q97*100), alpha=.3, 
                width=.2,position = position_dodge(0.5),size=.5)+
  scale_fill_manual(values=c("grey"),
                    labels=c( "Climate change 2016"),
                    name="")+
  xlab("Elevation (m)")+
  guides(fill=F)+
  ylab("Δ pred. prevalence (p.p.)")+
  scale_x_discrete(limits =c("Low (1387-2655)" ="Low (1387-2655)" , 
                             "Mid (2656-3132)" ="Mid (2656-3132)" ,
                             "High (3133-3486)"="High (3133-3486)"))


#======================================================================================================
##CALCULATING P.POINT DIFFERENCE ACROSS ELEVATION TERCILES
##======================================================================================================

dat=select(panel00,plot)
elev_bands=read_csv("Data/elevation_bands.csv")

df_elev=data.frame()

for (i in 1:10000){
  
  ##now run the monte carlo simulation
  d = draw[i,]
  modified_fe_mod = fe_mod
  modified_fe_mod$coefficients = d
  dat$vals_00 = predict(modified_fe_mod, newdata = panel00)
  dat$vals_95 = predict(modified_fe_mod, newdata = dat95)
  
  dat_range=dat%>%
    left_join(elevation)%>%
    pivot_longer(-c(plot, elevation))%>%
    filter(value>0.02)%>%
    group_by(name)%>%
    summarize(minelev=round(min(elevation)), maxelev=round(max(elevation)))%>%
    ungroup()
  
  maxarea_cc=elev_bands%>%
    filter(elev_bands$max_elev>=dat_range$maxelev[2] & elev_bands$max_elev<=dat_range$maxelev[1])%>%
    summarize(km2=sum(area),name="max00")
  
  minarea_cc=elev_bands%>%
    filter(elev_bands$min_elev>=dat_range$minelev[2] & elev_bands$min_elev<=dat_range$minelev[1])%>%
    summarize(km2=sum(area),name="min00")
  
  area_all=rbind(maxarea_cc, minarea_cc)
  
  
  df_elev=rbind(df_elev, area_all)
  
}

head(df_elev)

##estimating the mean and 95%CI of range expansion/contraction
rangedat=df_elev%>%
  group_by(name)%>%
  summarize(means=mean(km2),
            q25max=quantile(km2, .025),
            q97max=quantile(km2, .975))
rangedat



##checking the response forms between cubic and quadratic model
library(effects)
library(gridExtra)
source("Scripts/GLMMS_pathogen.R")

newpan=all_pan_dat%>%
  mutate(plot=factor(plot), year=factor(year))

fe_mod1 = lm(perinc~ vpd+I(vpd^2)+I(vpd^3)+dbh+density + plot+year,
               data = na.omit(newpan))

summary(fe_mod1)
AIC(fe_mod1, fe_mod_3)

plot(effect("vpd", fe_mod1))

feplot2=plot(effect("vpd", fe_mod1),axes=list(
  y=list(lab="P(prevalence)",ticks=list(at=c(-0.8,-.1,-.05, 0,.05,.1, .6))),
  x=list(vpd=list(lab="VPD (hPa)"))),main="FE panel model QUAD")

feplot2

fe_mod_3 = lm(perinc~ vpd+I(vpd^2)+dbh+density +plot+year,
             data = na.omit(newpan))     
summary(fe_mod_3)
plot(effect("vpd", fe_mod_3))  


feplot3=plot(effect("vpd", fe_mod_3),axes=list(
  y=list(lab="P(prevalence)",ticks=list(at=c(-0.8,-.1,-.05, 0,.05,.1, .6))),
  x=list(vpd=list(lab="VPD (hPa)"))),main="FE panel model CUBIC")

feplot3
grid.arrange(q95_plot,plot00,feplot2,plot95_3,plot00_3, feplot3, ncol=3)
           