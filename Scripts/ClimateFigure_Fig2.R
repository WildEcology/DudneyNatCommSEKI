## ---------------------------
##
## Script name: Fig 2 Climate Data
##
## Author: Dr. Joan Dudney
##
## Date Created: 2021-04-27
##
## Copyright (c) Joan Dudney, 2021
## Email: jdudney@berkeley.edu
##
## ---------------------------


library(ggpubr)
library(tidyverse)
library(patchwork)

select=dplyr::select
rename=dplyr::rename

theme_set(
  theme_bw(base_size = 17)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

##=================================================================================================================
##                            FIG 2 CLIMATE FIGURE
##    This code produces Fig 2 panels b-d, including VPD figure, annual precip figure
##    and temperature figure. Data was downloaded from https://prism.oregonstate.edu/ for
##    each plot between 1975-2016. The map of Fig 2 was produced in QGIS and data from that
##    map is available upon request.
##    
## ==================================================================================================================

##Reading in data

climdat=read_csv("Data/climatedata_2021.csv")

##Calculating the mean values across the time period
maxmin_temps=data.frame(mean=mean(climdat$meantmax), mean2=mean(climdat$meantmin))
meanppt=data.frame(mean=mean(climdat$ppt)) 
meanvpd=data.frame(mean=mean(climdat$meanvpd)) 


## Temperature figure
tempsdat=climdat%>%
  select(year, meantmin, meantmax)%>%
  pivot_longer(-year)

minmaxT_fig=ggplot(tempsdat, aes(x=year, y=value, color=name))+
  scale_color_manual(values = c("#d85555","#404080"),
                     labels=c("Maximum","Minimum"),
                     name="")+
  geom_point(alpha=.6)+
  geom_line(alpha=.5)+
  geom_hline(data=maxmin_temps, aes(yintercept = mean), 
             colour = "#d85555", linetype=2)+
  geom_hline(data=maxmin_temps, aes(yintercept = mean2), 
             colour = "#404080", linetype=2)+
  theme(strip.background = element_blank(),
        strip.placement = "outside")+
  ylab("Temperature (°C)")+
  scale_x_continuous(name="Year", breaks = scales::pretty_breaks(n = 6))+
  theme(legend.title = element_blank(),
        legend.position = c(.17,.5),
        legend.background=element_blank())

minmaxT_fig

##VPD figure

vpd_fig=ggplot(climdat, aes(x=year, y=meanvpd, alpha=.5))+
  scale_color_manual(values = c("black","black"))+
  geom_line(color="black")+
  geom_point()+
  geom_hline(data=meanvpd, aes(yintercept = mean), 
             colour = 'black', linetype=2)+
  theme(strip.background = element_blank(),
        strip.placement = "outside")+
  ylab("Max. VPD (hPa)")+
  geom_vline(xintercept = 1995, colour = "#495DB2", linetype=2)+
  geom_vline(xintercept = 1999, colour = "#495DB2", linetype=2)+
  geom_vline(xintercept = 2013, colour = "#495DB2", linetype=2)+
  geom_vline(xintercept = 2017, colour = "#495DB2", linetype=2)+
  guides(color=F)+
  guides(alpha=F)+
  scale_x_continuous(name="Year", breaks = scales::pretty_breaks(n = 6))+
  geom_rect(aes(xmin=1994, xmax=2000, ymin=6.5, ymax=7.5),
            fill="white")+
  geom_rect(aes(xmin=2010, xmax=2018, ymin=6.5, ymax=7.5),
            fill="white")+
  annotate("text", x = 1997, y = 7.0, size=c(5),label = "First survey", family ="Helvetica")+
  annotate("text", x = 2014, y = 7.0, size=c(5),label = "Second survey",family ="Helvetica")

vpd_fig


##Precipitation figure
precip_fig=ggplot(climdat, aes(x=year, y=ppt, fill=ppt,alpha=.8))+
  geom_rect(aes(xmin=2011.5, xmax=2015.5, ymin=0, ymax=max(ppt)),
            linetype=0, alpha=.01,fill="#CCCCCC")+
  geom_rect(aes(xmin=2006.5, xmax=2007.5, ymin=0, ymax=max(ppt)),
            linetype=0, alpha=.01,fill="#CCCCCC")+
  geom_bar(stat="identity", width=.5)+
  geom_hline(data=meanppt, aes(yintercept = mean), 
             colour = 'black', linetype=2)+
  theme(strip.background = element_blank(),
        strip.placement = "outside")+
  ylab("Precip. (mm)")+
  guides(fill=F)+
  guides(alpha=F)+
  scale_x_continuous(name="Year", breaks = scales::pretty_breaks(n = 6))

precip_fig

##combining figures using patchwork
vpd_fig/precip_fig/minmaxT_fig+ 
  plot_annotation(tag_levels="a") & theme(plot.tag.position = c(.05, 1),
   plot.tag = element_text(face = 'bold', size=15, family ="Helvetica", 
        hjust = -1, vjust = 0),text=element_text(family ="Helvetica"))


