## ---------------------------
##
## Script name: GLMM models of blister rust infections
##
## Author: Dr. Joan Dudney
##
## Date Created: 2021-04-27
##
## Copyright (c) Joan Dudney, 2021
## Email: jdudney@berkeley.edu
##
## ---------------------------
##
## Notes: 
##          
##Codes first and second survey 
##pathogen models with VPD and temperature
##Codes Figure 3 panels a-d
##
## ---------------------------

library(tidyverse)
library(lme4)
library(ggpubr)
library(reshape2)
library(ggthemes)
library(optimx)
library(MuMIn)
library(patchwork)
library(effects)
library(gridExtra)

select=dplyr::select
rename=dplyr::rename

##=================================================================================================================
##                      GLMM TESTING FOR CONFOUNDING IN VPD COEFFICIENT ESTIMATES
##    this code uses the first and second survey data to test whether VPD coefficient estimates
##    change through time, indicating that dispersal, local adaptation, or other unmeasured variables 
##    may be biasing estimate of the climate-disease relationship
## ==================================================================================================================

pathogen_dat=read_csv("Data/survey1and2_wpbr.csv")

##scaling across all data
era95=pathogen_dat%>%
  filter(!is.na(era95_inc_tot))%>%
  filter(!is.na(era95_status))%>%
  select(vpd95, era95_inc_tot, tree_id, species, plot, era95_ribes,slope, south, density95, era95_dbh, domsp)%>%
  set_colnames(c("vpd", "inc", "tree_id","species", "plot", "ribes",
                 "slope", "aspect", "density", "dbh", "domsp"))%>%
  mutate(time="first")


era00=pathogen_dat%>%
  filter(status==1)%>%
  filter(era95_inc_tot!=1)%>%
  select(inc_tot, vpd00, tree_id, species, plot, ribes, slope, south, density00, dbh, domsp)%>%
  set_colnames(c("inc", "vpd", "tree_id", "species", "plot", "ribes", 
                 "slope", "aspect", "density","dbh", "domsp"))%>%
  mutate(time="second")


merge.dat=era00%>%
  full_join(era95)%>%
  mutate(ribes=as.numeric(ribes))%>%
  mutate_at(scale, .vars = vars(-plot, -species,-inc, 
                                -tree_id, -ribes,-time, -domsp))%>%
  as.data.frame(.)


##model
models=glmer(inc~vpd*time+I(vpd^2)*time+ribes+aspect+slope+
               density*time+dbh*time+(1|plot)+(1|species),
             control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                    optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)),
             data = merge.dat, family=binomial)

summary(models)

models3=glmer(inc~vpd*time+I(vpd^2)*time+I(vpd^3)*time+ribes+aspect+slope+
               density*time+dbh*time+(1|plot)+(1|species),
             control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                    optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)),
             data = merge.dat, family=binomial)



##===================================================================================================
##                        GLMMs OF FIRST AND SECOND SURVEY BLISTER RUST INFECTIONS
##    this code is for the cross-sectional analyses of first and second survey tree-level infections
## ===================================================================================================

##FIRST SURVEY
dat95=filter(merge.dat, time=="first")

##quadratic model
mod95=glmer(inc~vpd+I(vpd^2)+density+ribes+dbh+slope+aspect+(1|plot)+(1|species),
            control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                   optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)),
            data = dat95, family=binomial)

summary(mod95)
r.squaredGLMM(mod95)

q95_plot=plot(effect("vpd", mod95),axes=list(
  y=list(lab="P(Infection)",ticks=list(at=c(-0.01,.001,.01,.02, .04))),
  x=list(vpd=list(lab="VPD (hPa)"))),main="First survey GLMM QUAD")


##cubic model

mod95_3=glmer(inc~vpd+I(vpd^2)+I(vpd^3)+density+ribes+dbh+slope+aspect+(1|plot)+(1|species),
              control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                     optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)),
              data = dat95, family=binomial)

summary(mod95_3)

plot95_3=plot(effect("vpd", mod95_3),axes=list(
  y=list(lab="P(Infection)",ticks=list(at=c(-0.01,.001,.01,.02, .04))),
  x=list(vpd=list(lab="VPD (hPa)"))),main="First survey GLMM CUBIC")

##SECOND SURVEY
dat00=filter(merge.dat, time=="second")


##quadratic model
mod00=glmer(inc~vpd+I(vpd^2)+density+ribes+dbh+slope+aspect+
              (1|plot)+(1|species),
            control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                   optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)),
            data = dat00, family=binomial)

summary(mod00)
r.squaredGLMM(mod00)

plot00=plot(effect("vpd", mod00),axes=list(
  y=list(lab="P(Infection)",ticks=list(at=c(-0.01,.001,.01,.02, .03))),
  x=list(vpd=list(lab="VPD (hPa)"))),main="Second survey GLMM QUAD")

##cubic model
mod00_3=glmer(inc~vpd+I(vpd^2)+I(vpd^3)+density+ribes+dbh+slope+aspect+
                (1|plot)+(1|species),
              control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                     optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)),
              data = dat00, family=binomial)


summary(mod00_3)

plot00_3=plot(effect("vpd", mod00_3),axes=list(
  y=list(lab="P(Infection)",ticks=list(at=c(-0.01,.001,.01,.02, .03))),
  x=list(vpd=list(lab="VPD (hPa)"))),main="Second survey GLMM CUBIC")
plot00_3



grid.arrange(q95_plot,plot95_3,plot00,plot00_3, ncol=2)


##=============================================================================================
##                        FIGURES FOR FIRST AND SECOND SURVEY MODELS
##                         this code creates Fig 4 panels b and d
## ============================================================================================

theme_set(
  theme_bw(base_size = 15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text=element_text(family ="Helvetica"))
)


#SECOND SURVEY
blist00=data.frame(summary(mod00)$coefficients)
blist00=round(blist00, digits=3)
names(blist00)[1:4]=c("estimate", "sterr", "zscore", "pvalue")
blist00$var=c("intercept",  "VPD","VPD^2","Density", "Ribes", "DBH",
              "Slope","Aspect")
blist00=blist00%>%
  filter(var!="intercept")

rust00=ggplot(blist00, aes(x=var, y=estimate, fill=var))+
  geom_bar(stat="identity", position=position_dodge(), width = 0.8, alpha=.8)+
  scale_fill_manual(values = c( "grey","grey","grey", "grey", "grey", "#d85555","#d85555")) +
  geom_errorbar(aes(ymin=estimate-sterr, ymax=estimate+sterr), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=ifelse(pvalue<.01, "**",
                             ifelse(pvalue<.05,"*", ""))), vjust=-2, 
            color="black", position=position_dodge(0.9), size=5)+
  guides(fill=FALSE)+
  ylim(-2,3)+
  xlab("Variables")+
  ylab("Coefficient estimates")+
  scale_x_discrete(labels = parse(text=c("Aspect"="Aspect","DBH"="DBH","Density"="Density",
                                         "Ribes"="Ribes","Slope"="Slope",
                                         "VPD"="VPD","VPD^2"="VPD^2")))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x = element_blank())
        #axis.text=element_text(size=11),
        #axis.title=element_text(size=12))
rust00

##FIRST SURVEY
blist95=data.frame(summary(mod95)$coefficients)
blist95=round(blist95, digits=3)
names(blist95)[1:4]=c("estimate", "sterr", "zscore", "pvalue")
blist95$var=c("intercept",  "VPD","VPD^2","Density", "Ribes", "DBH",
              "Slope","Aspect")
blist95=blist95%>%
  filter(var!="intercept")

rust95=ggplot(blist95, aes(x=var, y=estimate, fill=var))+
  geom_bar(stat="identity", position=position_dodge(), width = 0.8,alpha=.9)+
  scale_fill_manual(values = c( "grey","grey","grey", "grey", "grey","#404080","#404080")) +
  geom_errorbar(aes(ymin=estimate-sterr, ymax=estimate+sterr), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=ifelse(pvalue<.01, "**",
                             ifelse(pvalue<.05,"*", ""))), vjust=-2, 
            color="black", position=position_dodge(0.9), size=5)+
  guides(fill=FALSE)+
  ylim(-3,6)+
  ylab("Coefficient estimates")+
  scale_x_discrete(labels = parse(text=c("Aspect"="Aspect","DBH"="DBH","Density"="Density",
                                         "Ribes"="Ribes","Slope"="Slope",
                                         "VPD"="VPD","VPD^2"="VPD^2")))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x = element_blank())

rust95

rust95/
  rust00



##=============================================================================================
##                        FIGURES OF RAW BLISTER RUST AND TREE DATA
##                                creates Fig 3 panels a and c
## ============================================================================================


##FIRST SURVEY FIGURE

#Creating a second axis
ylim.prim2 = c(0, 1)   
ylim.sec2 = c(0, 2000) 

ylim.prim2 = c(0, 2000)   
ylim.sec2 = c(0, 1) 

b2 = diff(ylim.prim2)/diff(ylim.sec2)
a2 = b2*(ylim.prim2[1] - ylim.sec2[1])

era95_new=era95 %>%
  ggplot(aes(x=vpd))+
  geom_histogram(aes(y = stat((count / sum(count))*8)),
                 color="#477571",fill="#477571", alpha=0.5,  binwidth = .2)+
  stat_smooth(aes(x=vpd, y=inc), method = "glm", formula = y ~ x + I(x^2),
              method.args=list(family = "binomial"), fill= "#404080",color="#404080")+
  scale_y_continuous("Prop. infected stems",
                     sec.axis = sec_axis(~. *7031/8, name = "White pine stems (#)",
                                         breaks = seq(0,600,100)), n.breaks = 3)+
  ylab("Prop. infected & stem count")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(name="VPD (hPa)", seq(6,18,2))+
  ggtitle(label="First survey infections")+
  theme(axis.line.y.right = element_line(color = "#477571"), 
        axis.ticks.y.right = element_line(color = "#477571"),
        axis.text.y.right = element_text(color = "#477571"), 
        axis.title.y.right = element_text(color = "#477571"),
        axis.line.y.left = element_line(color = "#404080"), 
        axis.ticks.y.left = element_line(color = "#404080"),
        axis.text.y.left = element_text(color = "#404080"), 
        axis.title.y.left = element_text(color = "#404080"),
        plot.title = element_text(hjust = 0.5))

era95_new

era95_new_cube=era95 %>%
  ggplot(aes(x=vpd))+
  geom_histogram(aes(y = stat((count / sum(count))*8)),
                 color="#477571",fill="#477571", alpha=0.5,  binwidth = .2)+
  stat_smooth(aes(x=vpd, y=inc), method = "glm", formula = y ~ x + I(x^2)+I(x^3),
              method.args=list(family = "binomial"), fill= "#404080",color="#404080")+
  scale_y_continuous("Prop. infected stems",
                     sec.axis = sec_axis(~. *7031/8, name = "White pine stems (#)",
                                         breaks = seq(0,600,100)), n.breaks = 3)+
  ylab("Prop. infected & stem count")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(name="VPD (hPa)", seq(6,18,2))+
  ggtitle(label="First survey infections")+
  theme(axis.line.y.right = element_line(color = "#477571"), 
        axis.ticks.y.right = element_line(color = "#477571"),
        axis.text.y.right = element_text(color = "#477571"), 
        axis.title.y.right = element_text(color = "#477571"),
        axis.line.y.left = element_line(color = "#404080"), 
        axis.ticks.y.left = element_line(color = "#404080"),
        axis.text.y.left = element_text(color = "#404080"), 
        axis.title.y.left = element_text(color = "#404080"),
        plot.title = element_text(hjust = 0.5))

era95_new_cube

##===========================
##SECOND SURVEY FIGURE
##===========================

era00_new=era00 %>%
  ggplot(aes(x=vpd))+
  geom_histogram(aes(y = stat((count / sum(count))*5)),
                 color="#477571",fill="#477571", alpha=0.5,  binwidth = .2)+
  stat_smooth(aes(x=vpd, y=inc), method = "glm", formula = y ~ x + I(x^2),
              method.args=list(family = "binomial"), fill= "#d85555",color="#d85555")+
  scale_y_continuous("Prop. infected stems",
                     sec.axis = sec_axis(~. *5458/5, name = "White pine stems (#)",
                                         breaks = seq(0,600,100)), n.breaks = 6) +
  ylab("Prop. infected & stem count")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(name="VPD (hPa)", seq(6,18,2))+
  ggtitle(label="Second survey infections")+
  theme(axis.line.y.right = element_line(color = "#477571"), 
        axis.ticks.y.right = element_line(color = "#477571"),
        axis.text.y.right = element_text(color = "#477571"), 
        axis.title.y.right = element_text(color = "#477571"),
        axis.line.y.left = element_line(color = "#d85555"), 
        axis.ticks.y.left = element_line(color = "#d85555"),
        axis.text.y.left = element_text(color = "#d85555"), 
        axis.title.y.left = element_text(color = "#d85555"),
        plot.title = element_text(hjust = 0.5))

era00_new

era00_new_cube=era00 %>%
  ggplot(aes(x=vpd))+
  geom_histogram(aes(y = stat((count / sum(count))*5)),
                 color="#477571",fill="#477571", alpha=0.5,  binwidth = .2)+
  stat_smooth(aes(x=vpd, y=inc), method = "glm", formula = y ~ x + I(x^2)+I(x^3),
              method.args=list(family = "binomial"), fill= "#d85555",color="#d85555")+
  scale_y_continuous("Prop. infected stems",
                     sec.axis = sec_axis(~. *5458/5, name = "White pine stems (#)",
                                         breaks = seq(0,600,100)), n.breaks = 6) +
  ylab("Prop. infected & stem count")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(name="VPD (hPa)", seq(6,18,2))+
  ggtitle(label="Second survey infections")+
  theme(axis.line.y.right = element_line(color = "#477571"), 
        axis.ticks.y.right = element_line(color = "#477571"),
        axis.text.y.right = element_text(color = "#477571"), 
        axis.title.y.right = element_text(color = "#477571"),
        axis.line.y.left = element_line(color = "#d85555"), 
        axis.ticks.y.left = element_line(color = "#d85555"),
        axis.text.y.left = element_text(color = "#d85555"), 
        axis.title.y.left = element_text(color = "#d85555"),
        plot.title = element_text(hjust = 0.5))



era95_new+era00_new


##for diff_raw/fe_fig calling script FE_panelmodel.R
source("Scripts/FE_panelmodel.R")

allfigs=(era95_new/rust95)|(era00_new/rust00)|(diff_raw/fe_fig)

allfigs+plot_annotation(tag_levels="a") & theme(plot.tag.position = c(.05, 1),
                                                plot.tag = element_text(face = 'bold', size=15, family ="Helvetica", 
                                                                        hjust = 2, vjust = -.1),text=element_text(family ="Helvetica"))



##=================================================================================================================
##                                    GLMM MODELS USING TEMPERATURE   
##    
##    This code produces the model summaries in supplementary Table 5. It codes the quadratic model 
##    of first and second survey infections using maximum temperatures instead of VPD.
##
## ==================================================================================================================


##reading in temperature data
tempdata=read_csv("Data/temp.data.csv")

#first and second survey
temp95=filter(tempdata, time=="first")
temp00=filter(tempdata, time=="second")

##combinging temp data and estimating model with temperature

##FIRST SURVEY MODELS
comb_dat95=dat95%>%
  left_join(temp95)

mod95temps=glmer(inc~tmax+I(tmax^2)+density+ribes+dbh+slope+aspect+(1|plot)+(1|species),
                 control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                        optCtrl = list(method = "hjkb", starttests = FALSE, kkt = FALSE)),
                 data = comb_dat95, family=binomial)

summary(mod95temps)
r.squaredGLMM(mod95temps)

##SECOND SURVEY MODELS

comb_dat00=dat00%>%
  left_join(temp00)

mod00temps=glmer(inc~tmax+I(tmax^2)+density+ribes+dbh+slope+aspect+
                   (1|plot)+(1|species),
                 control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                        optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)),
                 data = comb_dat00, family=binomial)

summary(mod00temps)
r.squaredGLMM(mod00temps)



##=================================================================================================================
##                      
##                              CODE FOR SUPPLEMENTARY FIGURE 1
##    
## ==================================================================================================================

## set the theme for the figures
theme_set(
  theme_bw(base_size = 14)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text=element_text(family ="Helvetica"))
)



## =================================
##FIRST SURVEY FIGURE
## =================================

ylim.prim <- c(0, 900)  
ylim.sec <- c(0, 1) 

b = diff(ylim.prim)/diff(ylim.sec)

a = b*(ylim.prim[1] - ylim.sec[1])


era95_tempsfig=comb_dat95 %>%
  ggplot( aes(x=tmax))+
  geom_histogram(aes(x=vpd), color="#477571",fill="#477571", alpha=0.8,  binwidth = .2)+
  stat_smooth(aes(x=vpd, y=a+inc*b), method = "glm", formula = y ~ x + I(x^2),
              method.args=list(family = quasi(variance = "mu", link = "log")), size = 1,fill= "#404080",color="#404080")+
  scale_y_continuous("White pine stems (#)", 
                     sec.axis = sec_axis(~ (. - a)/b, name = "Proportion infected"), 
                     n.breaks = 8)+
  theme(axis.line.y.left = element_line(color = "#477571"), 
        axis.ticks.y.left = element_line(color = "#477571"),
        axis.text.y.left = element_text(color = "#477571"), 
        axis.title.y.left = element_text(color = "#477571"),
        plot.title = element_text(hjust = 0.5)
  )+
  xlab("Standardized maximum temperature (°C)")+
  ggtitle(label="First survey temperature")

era95_tempsfig


## =================================
##SECOND SURVEY FIGURE
## =================================

ylim.prim2 <- c(0, 2000)   
ylim.sec2 <- c(0, 1) 

b2 = diff(ylim.prim2)/diff(ylim.sec2)
a2 = b2*(ylim.prim2[1] - ylim.sec2[1])

era00_tempsfig= comb_dat00 %>%
  ggplot( aes(x=tmax))+
  geom_histogram( color="#477571",fill="#477571", alpha=0.8,  binwidth = .2)+
  stat_smooth(aes(y=a2+inc*b2), method = "glm", formula = y ~ x + I(x^2),
              method.args=list(family = quasi(variance = "mu", link = "log")), fill= "#d85555",color="#d85555")+
  scale_y_continuous("White pine stems (#)", 
                     sec.axis = sec_axis(~ (. - a2)/b2, name = "Proportion infected"), 
                     n.breaks = 8)+
  theme(axis.line.y.left = element_line(color = "#477571"), 
        axis.ticks.y.left = element_line(color = "#477571"),
        axis.text.y.left = element_text(color = "#477571"), 
        axis.title.y.left = element_text(color = "#477571"),
        plot.title = element_text(hjust = 0.5))+
  xlab("Standardized maximum temperature (°C)")+
  ggtitle(label="Second survey temperature")

era95_tempsfig/era00_tempsfig + 
  plot_annotation(tag_levels = 'a') & theme(plot.tag.position = c(0, .90),
    plot.tag = element_text(face = 'bold', size=12, family ="Helvetica", 
          hjust = -1, vjust = -1),text=element_text(family ="Helvetica"))

