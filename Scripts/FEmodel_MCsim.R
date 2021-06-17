## ---------------------------
##
## Script name: Fixed effects panel model, MC Simulation, Fig 4
##
## Author: Dr. Joan Dudney
##
## Date Created: 2021-04-27
##
## Copyright (c) Joan Dudney, 2021
## Email: jdudney@berkeley.edu
##
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

##=================================================================================================================
##                      FIXED EFFECTS PANEL MODEL AND MONTE CARLO SIMULATION
##    this code uses the first and second survey data to estimate a fixed effects panel model
##    and conduct a Monte Carlo simulation to account for model uncertainty and future climate change 
##    predictive uncertainty
## ==================================================================================================================


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
  mutate(year=factor(year), plot=factor(plot))%>%
  mutate(perinc2=perinc*100)


##===================================================================================================
##                            THE FIXED EFFECTS PANEL MODEL
##===================================================================================================

fe_mod = feols(perinc~ vpd+I(vpd^2)+dbh+density | plot + year,
               data = na.omit(all_pan_dat))
summary(fe_mod)

##===================================================================================================
##                            
##                                      MONTE CARLO SIMULATION
##      First we draw 10,000 x 4 coeficients from the FE panel model and variance-covariance matrix.
##      Then we randomly sample from 20 CMIP5 experiments of RCP4.5 2056-60
##      And predict change in blister rust prevalence across elevation terciles 10,000x
##
##====================================================================================================

##randomly sample from coefficient estimates and variance-covariance matrix

V_CR1 = vcovCR(fe_mod, cluster = c(all_pan_dat$plot), type = "CR1")
V_CR1=as.matrix(V_CR1)

coef_vector = fe_mod$coefficients

draw = rmvnorm(n = 100000, mean = coef_vector, sigma = V_CR1)


##replacing data00 vpd with vpd from 1995
vpd95=select(panel95, plot, vpd)

dat95=panel00%>%
  ungroup()%>%
  select(-vpd)%>%
  left_join(vpd95)%>%
  select(plot, vpd ,density, perinc , dbh ,year)

##future CMIP5 scenarios
filenames=read_csv("Data/filenames_short.csv")

##estimated % change from mean VPD second survey 
macadat=read_csv("Data/perChangeVPD_fut.csv")

macadiff=macadat%>%
  select(plot, idshort, perchange)%>%
  mutate(vpdchange=perchange+1)%>%
  select(-perchange)%>%
  filter(grepl("rcp45", idshort))

filenames=filter(filenames, grepl("rcp45", idshort))


#======================================================================================================
##CALCULATING P.POINT DIFFERENCE ACROSS ELEVATION TERCILES
##======================================================================================================


##DATA ELEVATION TERCILES
datelev_new=pathogen_dat%>%
  distinct(plot, elevation)%>%
  mutate(elev_third=ifelse(elevation<=2655.204, "Low (1387-2655)",
                           ifelse(elevation>2655.204&elevation<=3132.137, "Mid (2656-3132)",
                                  ifelse(elevation>=3132.137, "High (3133-3486)", "NA"))))

df_ranges=data.frame()

for (i in 1:10000){
  
  ##first, randomly sample different future climate change scenarios
  file=as.character(filenames[sample(nrow(filenames), 1), ])
  
  ##create dataframe
  futvpd=macadiff%>%
    filter(idshort==file)%>%
    mutate(plot=factor(plot))
  
  fut_dat = panel00 %>%##Future vpd
    left_join(futvpd)%>%
    mutate(vpdnew=vpd*vpdchange)%>%
    select(-vpd, -idshort, -vpdchange)%>%
    rename(vpd = vpdnew)
  
  ##now run the monte carlo simulation
  
  d = draw[i,]
  modified_fe_mod = fe_mod
  modified_fe_mod$coefficients = d
  datelev_new$vals_00 = predict(modified_fe_mod, newdata = panel00)
  datelev_new$vals_95 = predict(modified_fe_mod, newdata = dat95)
  datelev_new$vals_fut=predict(modified_fe_mod, newdata=fut_dat)
  
  dat_range=datelev_new%>%
    mutate(diff00=vals_00-vals_95, diff60=vals_fut-vals_95)%>%
    select(-c(vals_00, vals_95, vals_fut, elevation))%>%
    pivot_longer(-c(elev_third, plot))%>%
    group_by(elev_third, name)%>%
    summarize(mean=mean(value))
  
  tot=data.frame(dat_range)
  df_ranges=rbind(df_ranges, tot)
}


#======================================================================================================
## CREATING FIG. 4c, PREDICTED PREVALENCE ACROSS ELEVATION
##======================================================================================================

elevation=pathogen_dat%>%
  distinct(plot, elevation)%>%
  mutate(plot=factor(plot))

datelev=select(panel00,plot)%>%
  left_join(elevation)

df_elev_plots=data.frame()

for (i in 1:10000){
  
  ##first, randomly sample different future climate change scenarios
  file=as.character(filenames[sample(nrow(filenames), 1), ])
  
  ##making dataframe
  futvpd=macadiff%>%
    filter(idshort==file)%>%
    mutate(plot=factor(plot))
  
  fut_dat <- panel00 %>%##Future vpd
    left_join(futvpd)%>%
    mutate(vpdnew=vpd*vpdchange)%>%
    select(-vpd, -idshort, -vpdchange)%>%
    rename(vpd = vpdnew)
  
  ##now run the monte carlo simulation
  
  d = draw[i,]
  modified_fe_mod = fe_mod
  modified_fe_mod$coefficients = d
  datelev$vals_00 = predict(modified_fe_mod, newdata = panel00)
  datelev$vals_95 = predict(modified_fe_mod, newdata = dat95)
  datelev$vals_fut=predict(modified_fe_mod, newdata=fut_dat)
  datelev$iteration=i
  
  df_elev_plots=rbind(df_elev_plots, datelev)
}

head(df_elev_plots)

#======================================================================================================
## ESTIMATING RANGE EXPANSION AND CONTRACTION
##======================================================================================================

elev_bands=read_csv("Data/elevation_bands.csv")
dat=select(panel00,plot)
df_elev=data.frame()

for (i in 1:10000){
  
  ##first, randomly sample different future climate change scenarios
  file=as.character(filenames[sample(nrow(filenames), 1), ])
  
  ##making dataframe
  futvpd=macadiff%>%
    filter(idshort==file)%>%
    mutate(plot=factor(plot))
  
  fut_dat <- panel00 %>%##Future vpd
    left_join(futvpd)%>%
    mutate(vpdnew=vpd*vpdchange)%>%
    select(-vpd, -idshort, -vpdchange)%>%
    rename(vpd = vpdnew)
  
  ##now run the monte carlo simulation
  
  d = draw[i,]
  modified_fe_mod = fe_mod
  modified_fe_mod$coefficients = d
  dat$vals_00 = predict(modified_fe_mod, newdata = panel00)
  dat$vals_95 = predict(modified_fe_mod, newdata = dat95)
  dat$vals_fut=predict(modified_fe_mod, newdata=fut_dat)
  
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
  
  maxarea_fut=elev_bands%>%
    filter(elev_bands$max_elev>=dat_range$maxelev[2] & elev_bands$max_elev<=dat_range$maxelev[3])%>%
    summarize(km2=sum(area), name="futmax")
  
  minarea_fut=elev_bands%>%
    filter(elev_bands$min_elev>=dat_range$minelev[2] & elev_bands$min_elev<=dat_range$minelev[3])%>%
    summarize(km2=sum(area), name="futmin")
  
  area_all=rbind(maxarea_cc, minarea_cc, maxarea_fut, minarea_fut)
  
  
  df_elev=rbind(df_elev, area_all)
  
}

head(df_elev)

##estimating the mean and 95%CI of range expansion/contraction
rangedat=df_elev%>%
  group_by(name)%>%
  summarize(means=mean(km2),
            q25max=quantile(km2, .025),
            q97max=quantile(km2, .975))



#======================================================================================================
## ESTIMATING ELEVATION RANGE FOR MAP FIG 4
##======================================================================================================
dat=select(panel00,plot)

df_range=data.frame()

for (i in 1:10000){
  
  ##first, randomly sample different future climate change scenarios
  file=as.character(filenames[sample(nrow(filenames), 1), ])
  
  ##making dataframe
  futvpd=macadiff%>%
    filter(idshort==file)%>%
    mutate(plot=factor(plot))
  
  fut_dat <- panel00 %>%##Future vpd
    left_join(futvpd)%>%
    mutate(vpdnew=vpd*vpdchange)%>%
    select(-vpd, -idshort, -vpdchange)%>%
    rename(vpd = vpdnew)
  
  ##now run the monte carlo simulation
  
  d = draw[i,]
  modified_fe_mod = fe_mod
  modified_fe_mod$coefficients = d
  dat$vals_00 = predict(modified_fe_mod, newdata = panel00)
  dat$vals_95 = predict(modified_fe_mod, newdata = dat95)
  dat$vals_fut=predict(modified_fe_mod, newdata=fut_dat)
  
  dat_range=dat%>%
    left_join(elevation)%>%
    pivot_longer(-c(plot, elevation))%>%
    filter(value>0.02)%>%
    group_by(name)%>%
    summarize(minelev=round(min(elevation)), maxelev=round(max(elevation)))%>%
    ungroup()
  
  df_range=rbind(df_range, dat_range)
  
}

df_range

##calculating max/min from MC sim
ranges=df_range%>%
  group_by(name)%>%
  summarize(mean_min=mean(minelev),
            mean_max=mean(maxelev))
ranges


##=================================================================================================================
##       ##ESTIMATED ELEVATION SHIFT               
##==================================================================================================================
dat=select(panel00,plot)
df_rangeshift=data.frame()

for (i in 1:10000){
  
  ##first, randomly sample different future climate change scenarios
  file=as.character(filenames[sample(nrow(filenames), 1), ])
  
  ##making dataframe
  futvpd=macadiff%>%
    filter(idshort==file)%>%
    mutate(plot=factor(plot))
  
  fut_dat <- panel00 %>%##Future vpd
    left_join(futvpd)%>%
    mutate(vpdnew=vpd*vpdchange)%>%
    select(-vpd, -idshort, -vpdchange)%>%
    rename(vpd = vpdnew)
  
  ##now run the monte carlo simulation
  
  d = draw[i,]
  modified_fe_mod = fe_mod
  modified_fe_mod$coefficients = d
  dat$vals_00 = predict(modified_fe_mod, newdata = panel00)
  dat$vals_95 = predict(modified_fe_mod, newdata = dat95)
  dat$vals_fut=predict(modified_fe_mod, newdata=fut_dat)
  
  dat_range=dat%>%
    left_join(elevation)%>%
    pivot_longer(-c(plot, elevation))%>%
    filter(value>0.02)%>%
    group_by(name)%>%
    summarize(minelev=min(elevation), maxelev=max(elevation))%>%
    ungroup()%>%
    pivot_wider(names_from=name, values_from=c(minelev, maxelev))%>%
    mutate(diff00max=maxelev_vals_00-maxelev_vals_95,
           diff_futmax=maxelev_vals_fut-maxelev_vals_95,
           diff00min=minelev_vals_00-minelev_vals_95,
           diff_futmin=minelev_vals_fut-minelev_vals_95) %>% 
    select(diff00min,diff00max,diff_futmin,diff_futmax)
  
  tot=data.frame(dat_range)
  df_rangeshift=rbind(df_rangeshift, tot)
}

#write_csv(df_rangeshift, "Data_large/rangeshift_linear.csv")

linearR <- df_rangeshift %>% 
  ungroup() %>% 
  summarize(mdiffmax16=mean(diff00max),
            q25max16=quantile(diff00max, .025),
            q95max16=quantile(diff00max, .975),
            
            diffmin16=mean(diff00min),
            q25min16=quantile(diff00min, .025),
            q95min16=quantile(diff00min, .975),
            
            mdiffmaxfut=mean(diff_futmax),
            q25maxfut=quantile(diff_futmax, .025),
            q95maxfut=quantile(diff_futmax, .975),
            
            diffminfut=mean(diff_futmin),
            q25minfut=quantile(diff_futmin, .025),
            q95minfut=quantile(diff_futmin, .975))

