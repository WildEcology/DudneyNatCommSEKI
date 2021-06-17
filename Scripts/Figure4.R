## ---------------------------
##
## Script name: Fig 4 code
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
## this code produces Figure 4  
## panels a-d; Fig 6 panel e
## ---------------------------


library(ggpubr)
library(tidyverse)
library(patchwork)
library(ggpmisc)

select=dplyr::select
rename=dplyr::rename
group_by=dplyr::group_by


theme_set(
  theme_bw(base_size = 20)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5))
)

##reading in data
pathogen_dat=read_csv("Data/survey1and2_wpbr.csv")

#====================================================================================================
##FIG 4 VPD PANEL A
#====================================================================================================

macadat=read_csv("Data/perChangeVPD_fut.csv")

macadiff=macadat%>%
  select(plot, idshort, perchange)%>%
  mutate(vpdchange=perchange+1)%>%
  select(-perchange)%>%
  filter(grepl("rcp45", idshort))%>%
  group_by(plot)%>%
  summarize(mean.pc=mean(vpdchange))


vpd_dat=pathogen_dat%>%
  select(plot,vpd95, vpd00)%>%
  distinct()

all_vpd=vpd_dat%>%
  left_join(macadiff)%>%
  mutate(vpdfut=vpd00*mean.pc)%>%
  select(-mean.pc)

vpd_dat=all_vpd%>%
  pivot_longer(-plot)%>%
  mutate(name=factor(name, levels=c("vpd95", "vpd00", "vpdfut"), ordered=T))

vpdfig=ggplot(vpd_dat, aes(x=value, fill=name, color=name))+
  geom_density(alpha=.3)+
  scale_color_manual(values=c("white", "white", "white"))+
  scale_fill_manual(values=c("#2d2960","#d85555","#ff9538"),
                    breaks=c( "vpd95", "vpd00", "vpdfut"),
                    labels=c("Counterfactual","Climate change 2016", "RCP4.5 2056-60"),
                    name="")+
  guides(color=F)+
  ylab("Density (VPD)")+
  ggtitle("VPD distributions")+
  #ylim(0, .22)+
  scale_x_continuous(name="VPD (hPa)", breaks = scales::pretty_breaks(n = 10))+
  theme(legend.position = c(.73,.91), legend.text=element_text(size=12),
        legend.background=element_blank())

vpdfig


#====================================================================================================
## PREVALENCE FIG 4 PANEL B
#====================================================================================================

inc00=pathogen_dat%>% 
  filter(status==1)%>%
  group_by(plot, elevation)%>%
  summarize(perinc00=sum(inc_tot, na.rm=T)/sum(status, na.rm=T))

inc95=pathogen_dat%>% 
  filter(era95_status==1)%>%
  group_by(plot, elevation)%>%
  summarize(perinc95=sum(era95_inc_tot, na.rm=T)/sum(era95_status, na.rm=T))

combdat=inc95%>%
  full_join(inc00)%>%
  pivot_longer(-c(plot, elevation))

combdat$name = factor(combdat$name, levels = c("perinc95", "perinc00"), ordered = TRUE)


justinc=combdat %>%
  ggplot(aes(x=elevation, y=(value), color=name, fill=name))+
  geom_hline(yintercept = 0, linetype="dashed", color="grey")+
  geom_smooth(method = "glm", formula = y ~ x + I(x^2),
              method.args=list(family = "quasibinomial"))+
  ylab("Prop. infected")+
  scale_fill_manual(values=c("#2d2960", "#d85555"),
                    labels=c("First survey", "Second survey"),
                    name="")+
  scale_color_manual(values=c("#2d2960", "#d85555"),
                     labels=c("First survey", "Second survey"),
                     name="")+
  stat_peaks(alpha=.4)+
  stat_peaks(aes(y=value), alpha=.4,
             span=1)+
  ggtitle("Observed prevalence")+
  xlab("Elevation (m)")+
  theme(legend.position = c(.75,.92), legend.text=element_text(size=12),
        legend.background=element_blank())

justinc

#====================================================================================================
##FIG 4 PREDICTED PREVALENCE PANEL C
#====================================================================================================

mcdat=read_csv("Data_large/mcsim_elevation.csv")

##random sample figure
data_samp=data.frame(mcdat[sample(nrow(mcdat), 1000), ])

plotsamp_dat=data_samp%>%
  select(-iteration)%>%
  pivot_longer(-c(elevation, plot))

sample_elev_plot=ggplot(plotsamp_dat, aes(x=elevation, y=(value), color=name, fill=name))+
  geom_smooth(method="loess",span=1) +
  stat_peaks(alpha=.2)+
  scale_fill_manual(values=c("#2d2960","#d85555","#ff9538"),
                    breaks=c( "vals_95" ,"vals_00" , "vals_fut"),
                    labels=c("Counterfactual","Climate change 2016", "RCP4.5 2056-60"),
                    name="")+
  scale_color_manual(values=c("#2d2960","#d85555","#ff9538"),
                     breaks=c( "vals_95" ,"vals_00" , "vals_fut"),
                     labels=c("Counterfactual","Climate change 2016", "RCP4.5 2056-60"),
                     name="")+
  guides(color = guide_legend(override.aes = 
                                list(fill=c("#2d2960","#d85555","#ff9538"), 
                                     color=c("#2d2960","#d85555","#ff9538"))))+
  ylab("P(prop. infected)")+
  ggtitle("Predicted prevalence")+
  geom_hline(yintercept = 0, linetype="dashed", color="grey")+
  xlab("Elevation (m)")+
  guides(fill=F)+
  theme(legend.position = c(.73,.91), legend.text=element_text(size=12),
        legend.background=element_blank())

sample_elev_plot


#====================================================================================================
##FIG 4 PANEL D P.P. DIFFERENCES
#====================================================================================================


##MC simulated data

df_ranges=read_csv("Data_large/p.point_elev_MCdat.csv")

head(df_ranges)

newranges=df_ranges%>%
  group_by(elev_third, name)%>%
  summarize(allmeans=mean(mean),
            q25=quantile(mean, (.025)),
            q95=quantile(mean, (.975)))

newranges$name = factor(newranges$name, levels = c("diff00", "diff60"), ordered = TRUE)


elevationfig=ggplot(newranges, aes(elev_third, y=allmeans*100, fill=name))+
  geom_hline(yintercept = 0, linetype="dashed", color="grey")+
  geom_bar(stat="identity", position="dodge", width = 0.5, alpha=.7)+
  geom_errorbar(aes(ymin = q25*100, ymax = q95*100), alpha=.3, 
                width=.2,position = position_dodge(0.5),size=.5)+
  scale_fill_manual(values=c("#d85555","#ff9538"),
                    limits=c( "diff00", "diff60"),
                    labels=c( "Climate change 2016", "RCP4.5 2056-60"),
                    name="")+
  ggtitle("Predicted differences")+
  xlab("Elevation terciles (m)")+
  ylim(-38,15)+
  ylab("Δ P(prop. infected) (p.p.)")+
  scale_x_discrete(limits =c("Low (1387-2655)" ="Low (1387-2655)" , 
                             "Mid (2656-3132)" ="Mid (2656-3132)" ,
                             "High (3133-3486)"="High (3133-3486)"),
                   labels=c("Low", "Mid", "High"))+
  theme(legend.position = c(.70,.2), legend.text=element_text(size=12),
        legend.background=element_blank(),
        axis.text.x = element_text(size=17))

elevationfig


elevationfig=ggplot(df_ranges, aes(elev_third, y=mean*100, fill=factor(name)))+
  geom_hline(yintercept = 0, linetype="dashed", color="grey")+
  geom_point(position=position_jitterdodge(),alpha=0.01, aes(color=factor(name))) +
  geom_boxplot(alpha=.8, color="#2d2960", outlier.shape = NA)+
  scale_fill_manual(values=c("#d85555","#ff9538"),
                    limits=c( "diff00", "diff60"),
                    labels=c( "Climate change 2016", "RCP4.5 2056-60"),
                    name="")+
  scale_color_manual(values=c("#d85555","#ff9538"),
                    limits=c( "diff00", "diff60"),
                    labels=c( "Climate change 2016", "RCP4.5 2056-60"),
                    name="")+
  ggtitle("Predicted differences")+
  xlab("Elevation terciles (m)")+
  ylim(-38,15)+
  ylab("Δ P(prop. infected) (p.p.)")+
  scale_x_discrete(limits =c("Low (1387-2655)" ="Low (1387-2655)" , 
                             "Mid (2656-3132)" ="Mid (2656-3132)" ,
                             "High (3133-3486)"="High (3133-3486)"),
                   labels=c("Low", "Mid", "High"))+
  theme(legend.position = c(.70,.2), legend.text=element_text(size=12),
        legend.background=element_blank(),
        axis.text.x = element_text(size=17))

elevationfig


##COMBINING FIGURES

figs4=(vpdfig/justinc)|(sample_elev_plot/elevationfig)

figs4+
    plot_annotation(tag_levels = 'a') & theme(plot.tag.position = c(0, .90),
            plot.tag = element_text(face = 'bold', size=23, family ="Helvetica", 
                hjust = 0, vjust = -1),text=element_text(family ="Helvetica"))

#export 10x12

#====================================================================================================
##FIG 6 PANEL E
#====================================================================================================

justincpivot=combdat %>%
  ggplot(aes(x=elevation, y=value, color=name, fill=name))+
  geom_smooth(method = "glm", formula = y ~ x + I(x^2),
              method.args=list(family = "quasibinomial"))+
  ylab("Prop.(infected)")+
  scale_fill_manual(values=c("#2d2960", "#d85555"),
                    labels=c("First survey", "Second survey"),
                    name="")+
  scale_color_manual(values=c("#2d2960", "#d85555"),
                     labels=c("First survey", "Second survey"),
                     name="")+
  ggtitle("")+
  xlab("Elevation (m)")+
  theme_bw(base_size = 17)+
  scale_x_continuous(name="Elevation (m)", seq(1300, 3500,500),  expand = c(0, 0))+
  theme(legend.position = c(.55,.95), legend.text=element_text(size=12),
        legend.background=element_blank(),
        #axis.text.x =element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  coord_flip()

justincpivot

#ggsave("incidence_fig6_nopoint.pdf", justincpivot,width=3, height=8, bg = "transparent")


