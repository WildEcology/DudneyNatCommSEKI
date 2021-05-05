## ---------------------------
##
## Script name: FE panel model and figs
##
## Author: Dr. Joan Dudney
##
## Date Created: 2021-05-03
##
## Copyright (c) Joan Dudney, 2021
## Email: jdudney@berkeley.edu
##
## ---------------------------
##
## Notes: 
##          
##Codes FE panel model
##Codes Fig 3 panels e and f
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


select=dplyr::select
rename=dplyr::rename
group_by=dplyr::group_by

theme_set(
  theme_bw(base_size = 15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)))

##=================================================================================================================
##                                      FIXED EFFECTS PANEL MODEL 
##                      this code uses the first and second survey data to estimate a 
##                         fixed effects panel model and Figure 3 panels e and f
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
  mutate(year=factor(year), plot=factor(plot))


##===================================================================================================
##                            THE FIXED EFFECTS PANEL MODEL
##===================================================================================================

fe_mod = feols(perinc~ vpd+I(vpd^2)+dbh+density | plot + year,
               data = na.omit(all_pan_dat))
summary(fe_mod)

##===================================================================================================
##                               FIGURE 3 PANELS E AND F
##===================================================================================================

coefs=data.frame(summary(fe_mod)$coefficients, summary(fe_mod)$se, pvalue=c(.04, .001, .46, .03))
colnames(coefs)=c("estimate", "sterr", "pvalue")

##making figure
fe_coef1=round(coefs, digits=3)
fe_coef1$var=c("VPD",  "VPD^2","DBH", "Density")

##figure
fe_fig=ggplot(fe_coef1, aes(x=var, y=estimate, fill=var))+
  geom_bar(stat="identity", position=position_dodge(), width = 0.5, alpha=.7)+
  scale_fill_manual(values = c("grey", "grey", "black","black")) +
  geom_errorbar(aes(ymin=estimate-sterr, ymax=estimate+sterr), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=ifelse(pvalue<.01, "**",
                             ifelse(pvalue<0.045,"*", ""))), vjust=-3.8, 
            color="black", position=position_dodge(0.9), size=5)+
  #ggtitle(label="Panel model coefficient estimates")+
  guides(fill=FALSE)+
  ylim(-.03,.4)+
  xlab("Variables")+
  ylab("Coefficient estimates")+
  scale_x_discrete(labels = parse(text=c("DBH"="DBH","Density"="Density",
                                         "VPD"="VPD","VPD^2"="VPD^2")))+
  theme(axis.title.x = element_blank())
fe_fig


###change in infections
seconddat=pathogen_dat%>% 
  filter(status==1)%>%
  select(plot,tree_id,inc_tot, status, vpd00,dbh, density00)%>%
  group_by(plot, vpd00, density00)%>%
  dplyr::summarize(perinc00=sum(inc_tot, na.rm=T)/sum(status, na.rm=T),
                   dbh00=mean(dbh, na.rm=T))%>%
  ungroup()


firstdat=pathogen_dat%>% 
  filter(era95_status==1)%>%
  group_by(plot, vpd95, density95)%>%
  dplyr::summarize(perinc95=sum(era95_inc_tot, na.rm=T)/sum(era95_status, na.rm=T),
                   dbh95=mean(era95_dbh, na.rm=T))%>%
  ungroup()


comb_raw=seconddat%>%
  left_join(firstdat)%>%
  mutate(diffvpd=vpd00-vpd95, diffinc=perinc00-perinc95)%>%
  mutate(pos=ifelse(diffinc>0, "pos", "neg"))%>%
  na.omit()%>%
  mutate(quant=ifelse(vpd00<9.6, "0-33rd",
                      ifelse(vpd00>9.61&vpd00<11.7, "34th–67th",
                             ifelse(vpd00>=11.7, "68th–100th", "na"))))

diff_raw=comb_raw %>%
  ggplot(aes(x=diffvpd, y=diffinc, fill=quant, color=quant))+
  scale_color_manual(name= "Terciles", values=c("#477571", "#404080","#d85555" ))+
  scale_fill_manual(name="Terciles",values=c("#477571","#404080", "#d85555"))+
  geom_smooth(method="lm")+
  labs(y=expression(Δprevalence~(t[2]-t[1])),
       x=expression(ΔVPD~(t[2]-t[1]))) +
  ggtitle(label="Change in infections")+
  theme(legend.position = c(0.19, 0.23), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(fill = "transparent",color = "transparent"),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10),plot.title = element_text(hjust = 0.5))

diff_raw
