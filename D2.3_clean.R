#D2.3 Open source software for decoupling low and fast dynamics of data time series.
#Development of time series analysis methodology in an open source environment to 
#decompose time series of hyperspectral data or vegetation indices at daily, 
#seasonal and, if possible, interannual time scale.

#Steps:
#1)Quality check of the data 
#2)Gap filling of the night-time data
#3)Aggregattion of data with a custom period length
#4)Application of SSA - definining spectral bands, window length, number of components 


#Required packages:
if(!require(spectral.methods)) install.packages("spectral.methods")
library(spectral.methods)
if(!require(Rssa)) install.packages("Rssa")
library("Rssa")
if(!require(data.table)) install.packages("data.table")
library(data.table)
if(!require(zoo)) install.packages("zoo")
library(zoo)
if(!require(dplyr))install.packages("dplyr")
library(dplyr)
if(!require(FieldSpectroscopyCC)) install_github("tommasojulitta/FieldSpectroscopyCC") 
if(!require(FieldSpectroscopyDP)) install_github("tommasojulitta/FieldSpectroscopyDP") 
library(FieldSpectroscopyCC)
library(FieldSpectroscopyDP)
if(!require(timezone))install.packages("timezone")
library(timezone)
if(!require(GeoLight)) install.packages("GeoLight")
library(GeoLight)
if(!require(oce)) install.packages("oce")
library(oce)
if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)
if(!require(viridis)) install.packages("virids")
library(viridis)
if(!require(patchwork)) install.packages("patchwork")
library(patchwork)

#---user input:------------------------------------------------------
#output directory
out_directory<-choose.dir(caption = "select output folder")

#select time series for the decomposition
ts<-read.csv(choose.files(default = paste(out_directory, "\\*.*", sep = ""),caption = "select time series for the decomposition"),
             header=TRUE,sep=";",dec = ".",stringsAsFactors=FALSE,na.strings = "#N/D",fill=TRUE)

#---------------------------FILTER DATA------------------------------
#1. To filter out radiometric data acquired under unstable illumination conditions
#we use E_stability index computed as the percentage change between the first (E???1) 
#and the second (E???2) irradiance measurements within once cycle 

#---user input:------------------------------------------------------
#define the threshold for the percetange of irradiance change
stab_index<-1
#--------------------------------------------------------------------


#2. We also filter out saturated radiance data of both downward and upward channels of 
#full range and fluo range spectrometers
#sat value L [boolean]	Saturation value of downward channel. Fluo range
#sat value E [boolean]	Saturation value of upward channel 1. Fluo range
#sat value E2 [boolean]	Saturation value of upward channel 2. Fluo rang
#sat value L full [boolean]	Saturation value of downward channel. Full range
#sat value E full [boolean]	Saturation value of upward channel 1. Full range
#sat value E2 full [boolean]	Saturation value of upward channel 2. Full range


#4. In case of time series of fluorescence retrieved with SFM, is needed
#to filter out fluoresnece values when optimization of spectral fiting
#did not converge to a reasonable result (in this case A_sfm_conv==1, B_sfm_conv==1)

#Create an index of the good quality data
quality_ok<-which(ts$E_stability....<=stab_index & ts$E_stability.full.... <=stab_index
                  & ts$sat.value.L==0 & ts$sat.value.E==0 & ts$sat.value.E2==0 & ts$sat.value.L.full==0 
                  & ts$sat.value.E.full==0 & ts$sat.value.E2.full==0 & ts$A_sfm_conv==0 & ts$B_sfm_conv==0 &
                    ts$SIF_A_sfm..mW.m.2nm.1sr.1.>=0 & ts$SIF_B_sfm..mW.m.2nm.1sr.1.>=0)

ts_filt<-ts[quality_ok,]

#--------------------------------------------------------------------

#-------------------------AGGREGATE DATA-----------------------------
#---user input:------------------------------------------------------
#Define aggregation interval.Choose between "10 min", "20 min", "30 min"...
time_interval<-30
#--------------------------------------------------------------------
interval<-paste(as.character(time_interval), "min", sep=" ")
interval_backwards<-paste(as.character(-time_interval), "min", sep=" ")
interval_cut<-as.POSIXlt(as.character(cut(as.POSIXlt(ts_filt$datetime..UTC.,format="%d/%m/%Y %H:%M:%S", tz="GMT"),interval)),format="%Y-%m-%d %H:%M:%S", tz="GMT")
ts_filt_aggr<-aggregate(ts_filt,by=list(as.character(interval_cut)),FUN=mean,na.action=na.pass,na.rm=TRUE)
names(ts_filt_aggr)[1]<-"datetime"

#----------------------GAP-FILL NIGHT-TIME DATA----------------------
first_day_seq<-rev(seq(from=as.POSIXct(ts_filt_aggr$datetime[1], tz="GMT"), 
                       to=as.POSIXct(as.Date(ts_filt_aggr$datetime[1]),tz="GMT"), by=interval_backwards))

timestamp_seq<-seq(from=as.POSIXct(ts_filt_aggr$datetime[1], tz="GMT"),
                   to=as.POSIXct(format(as.POSIXct(ts_filt_aggr$datetime[length(ts_filt_aggr$datetime)], tz="GMT"), format = "%Y-%m-%d"),tz="GMT")+24*3600,
                   by=interval)

timestamp_seq_all<-c(as.character(first_day_seq[1:length(first_day_seq)-1]),as.character(timestamp_seq))
#---------------------------------------------------------------------
no_data<-vector()
no_data[which(timestamp_seq_all %in% ts_filt_aggr$datetime)] <-0
no_data[which(!timestamp_seq_all %in% ts_filt_aggr$datetime)]<-1

df_no_data<-data.frame(timestamp_seq_all, no_data); names(df_no_data)<-c("datetime", "no_data")
ts_filt_aggr_ext<-merge(ts_filt_aggr, df_no_data, by="datetime",all=TRUE)
#---------------------------------------------------------------------
#Choose variables of interest : SZA; PAR; datetime; doy.fraction, SIFA, SIFB, PRI
ts_chosen_var<-subset(ts_filt_aggr_ext, select=c("no_data","datetime","Lat","Lon","doy.dayfract","SZA","SAA","PAR..W.m.2.",
                                                 "SIF_A_sfm..mW.m.2nm.1sr.1.", 
                                                 "SIF_B_sfm..mW.m.2nm.1sr.1.", "PRI", "NDVI"))


ts_chosen_var$doy.dayfract<-DateToDOY(as.POSIXct(ts_chosen_var$datetime, tz="GMT")+2*3600)
ts_chosen_var$SZA<-zenith(solar(as.POSIXct(ts_chosen_var$datetime, tz="GMT")), lon=ts_chosen_var$Lon[25], lat=ts_chosen_var$Lat[25])
ts_chosen_var$SAA<-sunAngle(as.POSIXct(ts_chosen_var$datetime, tz="GMT"), lat=ts_chosen_var$Lat[25], lon=ts_chosen_var$Lon[25],
                            useRefraction=FALSE)$azimuth

ts_chosen_var$SIF_A_sfm..mW.m.2nm.1sr.1.[ts_chosen_var$SZA>=90]<-0
ts_chosen_var$SIF_B_sfm..mW.m.2nm.1sr.1.[ts_chosen_var$SZA>=90]<-0
ts_chosen_var$PRI[ts_chosen_var$SZA>=90]<-0
ts_chosen_var$PAR..W.m.2.[ts_chosen_var$SZA>=90]<-0
ts_chosen_var$Lat<-ts_chosen_var$Lat[25]       
ts_chosen_var$Lon<-ts_chosen_var$Lon[25]
ts_chosen_var$night_time[ts_chosen_var$SZA>=90]<-1
ts_chosen_var$night_time[ts_chosen_var$SZA<90]<-0

#-------------------Interpolate-------------------------------------------
ts_chosen_var$SIF_A_sfm..mW.m.2nm.1sr.1.<-na.approx(ts_chosen_var$SIF_A_sfm..mW.m.2nm.1sr.1.)
ts_chosen_var$SIF_B_sfm..mW.m.2nm.1sr.1.<-na.approx(ts_chosen_var$SIF_B_sfm..mW.m.2nm.1sr.1.)
ts_chosen_var$PRI<-na.approx(ts_chosen_var$PRI)
ts_chosen_var$PAR..W.m.2.<-na.approx(ts_chosen_var$PAR..W.m.2.)

#-------------------ADD NOISE TO NIGHT-TIME DATA-------------------------------------------
#---user input:----------------------------------------------------------------------------
#Choose variables to decompose:
variables_vector<-c(scan("", what = "char"))
variables_ssa<-list()
#------------------------------------------------------------------------------------------
for (v in 1:length(variables_vector)){
  cur<-ts_chosen_var[,c(variables_vector[v])]
  if (mean(cur[ts_chosen_var$night_time==0])>0){
    q<-0.95
  }else{
    q<-0.05
  }
  quant<-as.numeric(quantile(cur, probs=c(q)))
  noise<-rnorm(length(cur[ts_chosen_var$night_time==1]),0,abs(quant*0.01))
  cur[ts_chosen_var$night_time==1]<-noise
  nextcol<-data.frame(cur)
  colnames(nextcol)<-c(paste(variables_vector[v],"noise",sep=""))
  variables_ssa[[v]]<-nextcol
}

df_ssa = do.call(cbind,variables_ssa)

#----------------------DECOMPOSE TS WITH SSA--------------------------
#1. Set parameter for the decomposition
n_points_per_day<-24/(time_interval/60)
n_points<-nrow(ts_chosen_var)
n_days<-n_points/n_points_per_day

#bands
bands_list<-list(longterm=c(n_points/3,Inf), 
                 week_month=c(3*n_points_per_day,n_points/3), 
                 diurnal=c(n_points_per_day/3,3*n_points_per_day), 
                 subdiurnal=c(0,n_points_per_day/3))

#window length
wl_max<-n_points/2
wl_vector<-c(wl_max,n_points/3,3*n_points_per_day,n_points_per_day/3)

#number of components
n.comp_vector<-c(30,30,30,n_points_per_day/3)

#2. Run SSA
var_decom_list<-list()

for (i in 1:ncol(df_ssa)){
  cur<-df_ssa %>% select(i)
  cur_decomp<-filterTSeriesSSA(series = cur[,1],
                               borders.wl = bands_list,
                               M = wl_vector,
                               n.comp =n.comp_vector,
                               repeat.extr=c(2,1,1,1),
                               center.series=TRUE,
                               harmonics = c(0,0,0,0),
                               SSA.methods = "auto",
                               plot.spectra = FALSE,
                               open.plot =TRUE,
                               second.axis=TRUE)
  var_decom_list[[i]]<-cur_decomp$dec.series
  names(var_decom_list)[[i]]<-names(cur)
}

df_reconstr_ssa<-data.frame(t(do.call(rbind,var_decom_list)))
names(df_reconstr_ssa)<-c(paste(rep(names(var_decom_list), each = length(names(bands_list))), names(bands_list), sep = "_"))
ts_chosen_var<-data.frame(ts_chosen_var, df_reconstr_ssa)

#----------------------PLOT AND SAVE RESULTS--------------------------
#Tilter ts to exclude night-time data and interpolted data, and data acuired under sza > 70

ts_plots<-ts_chosen_var[ts_chosen_var$night_time==0 & ts_chosen_var$no_data==0 & ts_chosen_var$SZA<=70,]

#Plot F760 vs PRI original 
f760_pri_orig<-ggplot(data= ts_plots,aes(SIF_A_sfm..mW.m.2nm.1sr.1.,PRI,
                               colour=NDVI))+geom_point()+
  scale_color_viridis()+
  xlab("F 760")+ylab("PRI")+
  #ylim(-1.5,2)+xlim(-0.02,0.08)+
  #labs(subtitle = "Time scale: 30 min - 7 hh")+
  theme(plot.title = element_text(size=19,face="bold"))+
  theme(legend.title = element_text(size=19))+
  theme(legend.text = element_text(size=17))+
  theme(plot.subtitle = element_text(size=17))+
  theme(axis.text=element_text(size=17),axis.title=element_text(size=19))

#Plot F687 vs PRI original 
f687_pri_orig<-ggplot(data= ts_plots,aes(SIF_B_sfm..mW.m.2nm.1sr.1.,PRI,
                                         colour=NDVI))+geom_point()+
  scale_color_viridis()+
  xlab("F 687")+ylab("PRI")+
  #ylim(-1.5,2)+xlim(-0.02,0.08)+
  #labs(subtitle = "Time scale: 30 min - 7 hh")+
  theme(plot.title = element_text(size=19,face="bold"))+
  theme(legend.title = element_text(size=19))+
  theme(legend.text = element_text(size=17))+
  theme(plot.subtitle = element_text(size=17))+
  theme(axis.text=element_text(size=17),axis.title=element_text(size=19))

#Plot F760 vs F687 original 
f760_f687_orig<-ggplot(data= ts_plots,aes(SIF_A_sfm..mW.m.2nm.1sr.1.,SIF_B_sfm..mW.m.2nm.1sr.1.,
                                         colour=NDVI))+geom_point()+
  scale_color_viridis()+
  xlab("F 760")+ylab("F 687")+
  #ylim(-1.5,2)+xlim(-0.02,0.08)+
  #labs(subtitle = "Time scale: 30 min - 7 hh")+
  theme(plot.title = element_text(size=19,face="bold"))+
  theme(legend.title = element_text(size=19))+
  theme(legend.text = element_text(size=17))+
  theme(plot.subtitle = element_text(size=17))+
  theme(axis.text=element_text(size=17),axis.title=element_text(size=19))

#Plot F760 vs PRI sub-diurnal ssa 
f760_pri_sub_diur_ssa<-ggplot(data= ts_plots,aes(SIF_A_sfm..mW.m.2nm.1sr.1.noise_subdiurnal,
                                          PRInoise_subdiurnal,
                                          colour=NDVI))+geom_point()+
  scale_color_viridis()+
  xlab("F 760 sub-diurnal ssa")+ylab("PRI sub-diurnal ssa")+
  #ylim(-1.5,2)+xlim(-0.02,0.08)+
  #labs(subtitle = "Time scale: 30 min - 7 hh")+
  theme(plot.title = element_text(size=19,face="bold"))+
  theme(legend.title = element_text(size=19))+
  theme(legend.text = element_text(size=17))+
  theme(plot.subtitle = element_text(size=17))+
  theme(axis.text=element_text(size=17),axis.title=element_text(size=19))


#Plot F760 vs PRI sub-diurnal ssa 
f687_pri_sub_diur_ssa<-ggplot(data= ts_plots,aes(SIF_B_sfm..mW.m.2nm.1sr.1.noise_subdiurnal,
                                                  PRInoise_subdiurnal,
                                                  colour=NDVI))+geom_point()+
  scale_color_viridis()+
  xlab("F 687 sub-diurnal ssa")+ylab("PRI sub-diurnal ssa")+
  #ylim(-1.5,2)+xlim(-0.02,0.08)+
  #labs(subtitle = "Time scale: 30 min - 7 hh")+
  theme(plot.title = element_text(size=19,face="bold"))+
  theme(legend.title = element_text(size=19))+
  theme(legend.text = element_text(size=17))+
  theme(plot.subtitle = element_text(size=17))+
  theme(axis.text=element_text(size=17),axis.title=element_text(size=19))


plot_all<-f760_pri_orig+f687_pri_orig+f760_pri_sub_diur_ssa+f687_pri_sub_diur_ssa+plot_layout(guides ="collect")
ggsave(plot_all,path =out_directory, file=paste0("plot_all",".png"), 
       width = 30, height = 25, units = "cm")


