## ---------------------------
##
## Script name: Mortality and Ribes models and figures
##
## Author: Dr. Joan Dudney
##
## Date Created: 2021-04-30
##
## Copyright (c) Joan Dudney, 2021
## Email: jdudney@berkeley.edu
##
## ---------------------------

library(reshape)
library(ggeffects)
library(tidyverse)
library(lme4)
library(optimx)
library(ggpubr)
library(patchwork)
library(MuMIn)
library(sjPlot)
library(sjmisc)
library(sjlabelled)

select=dplyr::select
rename=dplyr::rename

##=================================================================================================================
##                      
##                                P(MORTALITY) AND RIBES SPP. MODELS
##                  This code estimates infected host mortality and Ribes spp. occurrence
##                  Predictions from these GLMM models are used to create Figure 6a,b
## ==================================================================================================================

##reading in data
pathogen_dat=read_csv("Data/survey1and2_wpbr.csv")

##cleaning data and standardizing variables 
era95=pathogen_dat%>% 
  filter(!is.na(era95_inc_tot))%>%
  filter(!is.na(era95_status))

scaled95=era95%>% 
  dplyr::select(plot, era95_inc_tot, species, domsp, south, slope, vpd00, 
                era95_dbh, era95_ribes,elevation,density95, status)%>%
  mutate_at(scale, .vars = vars(-plot, -domsp,-species,-era95_inc_tot, 
                                -south, -status, -era95_ribes))%>%
  as.data.frame(.)%>%
  mutate(mortality=ifelse(status==1, 0, 1))%>%
  filter(!is.na(status))


##=================================================================================================================
##                                    INFECTED HOST MORTALITY GLMMS
##                                  estimated mortality using two GLMMs, 
##                      one that estimates with VPD and the other with elevation
## ==================================================================================================================

##selecting only infected hosts in the first survey
inf_host=scaled95%>%
  filter(era95_inc_tot==1)


##GLMM using VPD
infmort.vpd=glmer(mortality ~ era95_dbh+vpd00+density95+slope+south+(1|plot)+(1|species),
                  control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                         optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)),
                  family = binomial,data = inf_host)
summary(infmort.vpd)
r.squaredGLMM(infmort.vpd)

##GLMM using elevation
infmort.elev=glmer(mortality ~ era95_dbh+elevation+density95+slope+south+(1|plot)+(1|species),
                   control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                          optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)),
                   family = binomial, data = inf_host)
summary(infmort.elev)
r.squaredGLMM(infmort.elev)


##=================================================================================================================
##                              PREDICTED INFECTED HOST MORTALITY
##                               Using elevation GLMM of mortality
##                to get the predicted probability of mortality across elevation
## ==================================================================================================================

##not scaling elevation to match elevation in m in Fig 5
scaled95_notelev=era95%>% 
  dplyr::select(plot, era95_inc_tot, species, domsp, south, slope, vpd00, 
                era95_dbh, era95_ribes,elevation,density95, status)%>%
  mutate_at(scale, .vars = vars(-plot, -domsp,-species,-era95_inc_tot, 
                                -south, -status, -era95_ribes, -elevation))%>%
  as.data.frame(.)%>%
  mutate(mortality=ifelse(status==1, 0, 1))%>%
  filter(!is.na(status))

##selecting only infected hosts in the first survey
inf_host_elev=scaled95_notelev%>%
  filter(era95_inc_tot==1)


##GLMM with unscaled elevation
infmort.elev1=glmer(mortality ~ era95_dbh+elevation+density95+slope+south+(1|plot)+(1|species),
                   family = binomial, data = inf_host_elev)

##Marginal effect of elevation
predicted_inf=ggpredict(infmort.elev1, terms="elevation[all]")
plot(predicted_inf)



##=================================================================================================================
##                                   RIBES SPP. GLMS
##                        estimated Ribes occurrence using two glms, 
##                one that estimates with VPD and the other with elevation
## ==================================================================================================================


##Ribes data
ribes=pathogen_dat%>%
  select(plot, ribes, domsp,elevation,vpd00, slope, south)%>%
  mutate_at(scale, .vars = vars(-plot,-ribes, -domsp,-south))%>%
  as.data.frame(.)%>%
  distinct()%>%
  mutate(ribes=as.numeric(ribes))%>%
  na.omit()
  

##Ribes glm with elevation
ribes.elev.mod=glm(ribes~elevation+slope+south, family=binomial, data=ribes)
summary(ribes.elev.mod)


##Ribes glm with VPD
ribes.mod=glm(ribes~vpd00+slope+south, family=binomial, data=ribes)
summary(ribes.mod)

##Model summary table
tab_model(ribes.mod,ribes.elev.mod,
          transform = NULL,show.aic = T,  show.re.var = F, show.est=T, 
          show.se=T,auto.label = F, show.ci = F,
          pred.labels = c("Intercept", "VPD", 
                          "Slope","Aspect", "Elevation"),
          dv.labels = c("Ribes spp. with VPD", "Ribes spp. with elevation"),
          #string.pred = "Coeffcient",
          string.est = "Coefficient",
          string.se  = "Std. Error",
          string.p = "P-Value")


##=================================================================================================================
##                              PREDICTED RIBES OCCURRENCE 
##                               Using elevation GLM 
##                to get the predicted probability of Ribes across elevation
## ==================================================================================================================

##Ribes data
ribes_dat_elev=pathogen_dat%>%
  select(plot, ribes, domsp,elevation,slope, south)%>%
  as.data.frame(.)%>%
  distinct()%>%
  mutate(ribes=as.numeric(ribes))%>%
  na.omit()

##Ribes glm with elevation
ribes.elev.mod1=glm(ribes~elevation+slope+south, family=binomial, data=ribes_dat_elev)


ribes_elev=ggpredict(ribes.elev.mod1, terms="elevation[all]")
plot(ribes_elev)


#=================================================================================================================
##                                   FIGURES 6A,B 
##                      here predicted values of Ribes and mortality
##               are combined into one dataset and used to graph Figures 6a,b
## ==================================================================================================================

##Ribes data
ribesdat=ribes_elev%>%
  mutate(name="ribes")%>%
  set_colnames(c("elevation", "predicted", "std.error" ,"conf.low" ,"conf.high", "group" , "name"))

##Combining ribes and mortality data
mortdatpre=predicted_inf%>%
  mutate(name="mort")%>%
  set_colnames(c("elevation", "predicted" , "std.error" ,"conf.low" ,"conf.high" , "group" , "name"  ))%>%
  full_join(ribesdat)%>%
  mutate(newname=ifelse(name=="mort", "Mortality", "  Ribes"))


##FIGURE
newfig=ggplot(mortdatpre, aes(elevation, predicted, fill=newname, color=newname)) +
  facet_wrap(newname~.)+
  geom_smooth(se=F)+
  scale_color_manual(values =  c( "#2B82A1","#23404B"),  name="")+
  scale_fill_manual(values =  c("#2B82A1","#23404B"),  name="")+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, group=newname, fill=newname), 
              alpha = .1, color=c("lightgrey"))+
  xlab("Elevation (m)")+
  guides(fill=F, color=F)+
  ylab("Predicted probability")+
  scale_y_continuous(limits=c(0, 1.03),expand = c(0, 0),
                     labels = scales::number_format(accuracy = 0.1,
                                                    decimal.mark = '.'))+
  scale_x_continuous(name="Elevation (m)", seq(1300, 3500,500),  expand = c(0, 0))+
  scale_y_continuous(name="Predicted probability", seq(0, 1,.4)) +
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank(),
        strip.text.x = element_text(size = 15, colour = "black", angle = 0),
        legend.background = element_blank(),
        legend.position = c(.25,.15),
        strip.background = element_blank(),
        panel.spacing = unit(3, "lines"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  coord_flip()

newfig

#ggsave("newfig5.mortribes.pdf", newfig,width=4, height=8, bg = "transparent")

