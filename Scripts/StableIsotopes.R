
## ---------------------------
##
## Script name: Needle stable isotope study
##
## Author: Dr. Joan Dudney
##
## Date Created: 2021-05-01
##
## Copyright (c) Joan Dudney, 2021
## Email: jdudney@berkeley.edu
##
## ---------------------------
##
## Notes: This code produces the following:
##   
##    1. Runs the sugar pine drought study analysis
##    2. Creates Figure 5
##    3. Creates Figure 6 panels c,d
## ---------------------------


library(optimx)
library(tidyverse)
library(lme4)
library(ggtext) 
library(ggpubr)
library(patchwork)

##=================================================================================================================
##                              SUGAR PINE STABLE ISOTOPE STUDY ANALYSIS
##                This code runs six repeated measures analysis of variance models
##                including c13, c15, %carbon, %nitrogen, needle length, and #fascicles
##==================================================================================================================

isodat=read_csv("Data/isotopes_sugarpine.csv")

##Checking whether data are normally distributed
ggdensity(isodat$delta_c) ##c13
ggqqplot(isodat$delta_c)
ggdensity(isodat$n) ##n15
ggqqplot(isodat$n)
ggdensity(isodat$pern) #percent nitrogen
ggqqplot(isodat$pern)
ggdensity(isodat$perc) #percent carbon
ggqqplot(isodat$perc)
ggdensity(isodat$mean_num) #mean # fascicles
ggqqplot(isodat$mean_num)
ggdensity(isodat$mean_length) #mean needle length
ggqqplot(isodat$mean_length)

## conducting tests
analysis_dat=isodat%>%
  mutate(year=factor(year), treeid=factor(treeid), pair=factor(pair), 
         wpbr=factor(wpbr), drought=factor(drought))

##Stable carbon isotope model
##Anova
aovc <- aov(delta_c ~ wpbr+Error(pair/year), data = analysis_dat)
summary(aovc)


##check with glm which should produce the same results
glmmod=glm(delta_c~wpbr+year+pair, data=analysis_dat)


##Stable nitrogen isotope model
aovn=aov(n ~ wpbr+ Error(pair/year), data = analysis_dat)
summary(aovn)

##Percent needle carbon model
aovperc=aov(perc ~ wpbr+Error(pair/year), data = analysis_dat)
summary(aovperc)

summary(glm(perc~wpbr+year+pair, data=analysis_dat))

##Percent nitrogen model
aovpern=aov(pern ~ wpbr+Error(pair/year), data = analysis_dat)
summary(aovpern)
glance(aovpern) 

summary(glm(pern~wpbr+year+pair, data=analysis_dat))

##Mean needle length model
aovl=aov(mean_length ~ wpbr+Error(pair/year), data = analysis_dat)
summary(aovl)

summary(glm(mean_length~wpbr+year+pair, data=analysis_dat))

##Mean fascicle number
aovnum=aov(mean_num ~ wpbr+Error(pair/year), data = analysis_dat)
summary(aovnum)

##=================================================================================================================
##                              MORTALITY MODEL OF SUGAR PINE
##==================================================================================================================

##Mortality
mortanaly=analysis_dat%>% 
  distinct(treeid, wpbr, dbh, pair, mortality)%>%
  mutate_at(scale, .vars = vars(dbh))%>%
  as.data.frame(.)

summary(glmer(mortality~wpbr+(1|pair), data=mortanaly, 
              control=glmerControl(optimizer="bobyqa"),family = binomial))


##=================================================================================================================
##                               PRODUCES FIGURE 5
##                This code produces all panels for Figure 5 
##                        of sugar pine drought study
##==================================================================================================================

theme_set(
  theme_bw(base_size = 15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text=element_text(family ="Helvetica"))
)


isofign=ggplot(isodat, aes(x=Rust, y=n, fill=Rust, color=Rust))+
  geom_boxplot(width=0.5,position=position_dodge(1), alpha=.5)+
  geom_jitter(position=position_jitter(0.2),alpha=.1,aes(color=Rust))+
  scale_color_manual(values=c("#477571","#EAB948"))+
  scale_fill_manual(values=c("#477571","#EAB948"))+
  labs(y="Needle *&delta;*<sup>15</sup>N (&permil;)",
       x="Blister rust")+
  guides(fill=F, color=F)+
  annotate("text", y=1,x=1.5, label="F = 4.85, P = 0.03", 
           alpha=1,color="black", size=2.5)+
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown()
  )

isofign


pernfig=ggplot(isodat, aes(x=Rust, y=pern, fill=Rust, color=Rust))+
  geom_boxplot(width=0.5,position=position_dodge(1), alpha=.5)+
  geom_jitter(position=position_jitter(0.2),alpha=.1,aes(color=Rust))+
  scale_color_manual(values=c("#477571","#EAB948"))+
  scale_fill_manual(values=c("#477571","#EAB948"))+
  ylab("% Nitrogen")+
  xlab("Blister rust")+
  guides(fill=F)+
  guides(color=F)+
  annotate("text", y=1.7,x=1.5, label="F = 4.71, P = 0.03", 
           alpha=1,color="black", size=2.5)


pernfig

needfig=ggplot(isodat, aes(x=Rust, y=mean_length, fill=Rust, color=Rust))+
  geom_boxplot(width=0.5,position=position_dodge(1), alpha=.5)+
  geom_jitter(position=position_jitter(0.2),alpha=.1,aes(color=Rust))+
  scale_color_manual(values=c("#477571","#EAB948"))+
  scale_fill_manual(values=c("#477571","#EAB948"))+
  ylab("Needle length (cm)")+
  xlab("Blister rust")+
  guides(fill=F, color=F)+
  annotate("text", y=12,x=1.5, label="F = 12.74, P < 0.001", 
           alpha=1,color="black", size=2.5)


needfig

numfig=ggplot(isodat, aes(x=Rust, y=mean_num, fill=Rust, color=Rust))+
  geom_boxplot(width=0.5,position=position_dodge(1), alpha=.5)+
  scale_color_manual(values=c("#477571","#EAB948"))+
  scale_fill_manual(values=c("#477571","#EAB948"))+
  geom_jitter(position=position_jitter(0.2),alpha=.1,aes(color=Rust))+
  ylab("Fascicle #")+
  xlab("Blister rust")+
  guides(fill=F, color=F)+
  annotate("text", y=57,x=1.5, label="F = 26.05, P < 0.001", 
           alpha=1,color="black", size=2.5)

numfig


isofigyearsc=ggplot(isodat, aes(x=year, y=delta_c, color=Rust, fill=Rust))+
  geom_smooth(alpha=.3)+
  geom_jitter(width = .2, alpha=.2)+
  scale_color_manual(values=c("#477571","#EAB948"),
                     labels = c("No", "Yes"),
                     limits = c("No", "Yes"),
                     name="Blister rust") +
  scale_fill_manual(values=c("#477571","#EAB948"),
                    labels = c("No", "Yes"),
                    limits = c("No", "Yes"),
                    name="Blister rust")+
  labs(y="Needle *&delta;*<sup>13</sup>C (&permil;)",
       x="Year")+
  #guides(fill=F)+
  scale_x_continuous(name="Year", seq(2012,2017,1))+
  annotate("text", y=-25,x=2016, label="F = 11.28, P = 0.001",
      alpha=1,color="black", size=3)+
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
    legend.position = c(.85,.15),
    legend.background = element_blank()
  )

isofigyearsc



(isofigyearsc|isofign+pernfig+numfig+needfig)+
  plot_annotation(tag_levels="a") & theme(plot.tag.position = c(.05, 1),
          plot.tag = element_text(face = 'bold', size=12, family ="Helvetica", 
                 hjust = -1, vjust = -.25),text=element_text(family ="Helvetica"))

(isofigyearsc|isofign+pernfig+numfig+needfig)+
  plot_annotation(tag_levels="a") & theme(plot.tag = element_text(face = 'bold', size=12, family ="Helvetica"), 
                text=element_text(family ="Helvetica", size = 13))



##reading in needle length data
needles=read_csv("Data/needle_lengths.csv")

##WHITEBARK figure
whitebarkfig=ggplot(filter(needles, name=="whitebark"), aes(x=factor(year), y=length, fill=drought))+
  #stat_compare_means()+
  geom_jitter(position=position_jitter(0.2),alpha=.5,aes(color=drought))+
  geom_boxplot(alpha=.7, width=.4, lwd=.4)+
  scale_color_manual(values=c("#EAB948","#477571"))+
  scale_fill_manual(values=c( "#EAB948","#477571"))+
  ylab("Length (cm)")+
  xlab("Drought")+
  guides(color=F)+
  theme_bw(base_size = 20)+
  ggtitle("Whitebark")+
  scale_x_discrete(name="", labels=c("2011",  "","2013", "", "2015", "", "2017"))+
  theme(legend.position="bottom",legend.title = element_blank(),
        legend.margin=margin(-10, 0, 0, 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=20),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

whitebarkfig


##SUGAR PINE figuree
sugarpinefig=ggplot(filter(needles, name=="sugar pine"), aes(x=factor(year), y=length, fill=drought))+
  geom_jitter(position=position_jitter(0.2),alpha=.5,aes(color=drought))+
  geom_boxplot(alpha=.7, width=.4, lwd=.4)+
  scale_color_manual(values=c("#EAB948","#477571"))+
  scale_fill_manual(values=c("#EAB948","#477571"))+
  ylab("Length (cm)")+
  xlab("Drought")+
  guides(color=F)+
  theme_bw(base_size = 25)+
  ggtitle("Sugar pine")+
  scale_x_discrete(name="", labels=c("","2013", "", "2015", "", "2017"))+
  theme(legend.position = "bottom",legend.title = element_blank(),
        legend.margin=margin(-10, 0, 0, 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=25),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

sugarpinefig








