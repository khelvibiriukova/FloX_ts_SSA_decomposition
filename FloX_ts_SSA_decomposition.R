#Open source code for decoupling low and fast dynamics in time series of hyperspectral data.
#Decomposition of the time series of solar-induced chlorophyll fluorescence 
#retrieved in O2A and O2B absorption bands (SIF A/SIF B), and Photochemical 
#reflectance index (PRI) with Singular Spectrum Analysis (SSA).

#The script works only with the output of the FloX processing chain (https://www.jb-hyperspectral.com/products/flox/) on Windows OS

#Steps:
#1) Quality check of the spectral data and filtering
#2) Aggregation of data by a 30 min period
#3) Gap filling of the night-time data
#4) Adding noise to the night-time data
#4) Application of SSA on SIF A, SIF B and PRI - decomposition on long-term (seasonal), diurnal and sub-diurnal time scales
#5) Visualization of the results


#Required packages:
if(!require(spectral.methods)) install.packages("spectral.methods")
library(spectral.methods)
if(!require(data.table)) install.packages("data.table")
library(data.table)
if(!require(zoo)) install.packages("zoo")
library(zoo)
if(!require(dplyr))install.packages("dplyr")
library(dplyr)
if (!require(devtools)) install.packages("devtools")
library(devtools)
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
#we use E_stability index computed as the percentage change between the first (E1) 
#and the second (E2) irradiance measurements within once cycle 

#---user input:------------------------------------------------------
#define the threshold for the percentage of irradiance change
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


#Create an index of the good quality data
quality_ok<-which(ts$E_stability....<=stab_index & ts$E_stability.full.... <=stab_index
                  & ts$sat.value.L==0 & ts$sat.value.E==0 & ts$sat.value.E2==0 & ts$sat.value.L.full==0
                  & ts$sat.value.E.full==0 & ts$sat.value.E2.full==0 & 
                    ts$SIF_A_sfm..mW.m.2nm.1sr.1.>=0 & ts$SIF_B_sfm..mW.m.2nm.1sr.1.>=0)


ts_filt<-ts[quality_ok,]

#--------------------------------------------------------------------

#-------------------------AGGREGATE DATA-----------------------------
#---user input:------------------------------------------------------
#1. Define aggregation interval.Choose between "10 min", "20 min", "30 min"...
time_interval<-30
#--------------------------------------------------------------------
interval<-paste(as.character(time_interval), "min", sep=" ")
interval_backwards<-paste(as.character(-time_interval), "min", sep=" ")
interval_cut<-as.POSIXlt(as.character(cut(as.POSIXlt(ts_filt$datetime..UTC.,tryFormats=c("%d/%m/%Y %H:%M:%S","%d/%m/%Y %H:%M",
                                                                                         "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M"),tz="GMT"),interval)),format="%Y-%m-%d %H:%M:%S", tz="GMT")
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
ts_chosen_var<-subset(ts_filt_aggr_ext, select=c("no_data","datetime","Lat","Lon","doy.dayfract","SZA","PAR..W.m.2.",
                                                 "SIF_A_sfm..mW.m.2nm.1sr.1.",
                                                 "SIF_B_sfm..mW.m.2nm.1sr.1.", "PRI", "NDVI"))

x<-as.POSIXct(ts_chosen_var$datetime, tz="GMT")
doy<-difftime(x,as.POSIXct(paste(format(x, "%Y"),"-01-01 00:00",sep=""), tz="GMT"),units='days')
ts_chosen_var$doy.dayfract<-as.numeric(doy)+1
ts_chosen_var$SZA<-zenith(solar(as.POSIXct(ts_chosen_var$datetime, tz="GMT")), lon=unique(ts_chosen_var$Lon)[!is.na(unique(ts_chosen_var$Lon))],
                          lat=unique(ts_chosen_var$Lat)[!is.na(unique(ts_chosen_var$Lat))])

ts_chosen_var$SIF_A_sfm..mW.m.2nm.1sr.1.[ts_chosen_var$SZA>=90]<-0
ts_chosen_var$SIF_B_sfm..mW.m.2nm.1sr.1.[ts_chosen_var$SZA>=90]<-0
ts_chosen_var$PRI[ts_chosen_var$SZA>=90]<-0
ts_chosen_var$PAR..W.m.2.[ts_chosen_var$SZA>=90]<-0
ts_chosen_var$Lat<-unique(round(ts_chosen_var$Lat, digits = 2))[!is.na(unique(round(ts_chosen_var$Lat, digits = 2)))] 
ts_chosen_var$Lon<-unique(round(ts_chosen_var$Lon, digits = 2))[!is.na(unique(round(ts_chosen_var$Lon, digits = 2)))]
ts_chosen_var$night_time[ts_chosen_var$SZA>=90]<-1
ts_chosen_var$night_time[ts_chosen_var$SZA<90]<-0

#-------------------Interpolate-------------------------------------------
ts_chosen_var$SIF_A_sfm..mW.m.2nm.1sr.1.<-na.approx(ts_chosen_var$SIF_A_sfm..mW.m.2nm.1sr.1.)
ts_chosen_var$SIF_B_sfm..mW.m.2nm.1sr.1.<-na.approx(ts_chosen_var$SIF_B_sfm..mW.m.2nm.1sr.1.)
ts_chosen_var$PRI<-na.approx(ts_chosen_var$PRI)
ts_chosen_var$PAR..W.m.2.<-na.approx(ts_chosen_var$PAR..W.m.2.)

#-------------------ADD NOISE TO NIGHT-TIME DATA-------------------------------------------
#---user input:----------------------------------------------------------------------------
#1. Choose variables to decompose:
variables_vector<-c("SIF_A_sfm..mW.m.2nm.1sr.1.", "SIF_B_sfm..mW.m.2nm.1sr.1.","PRI")

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
bands_list<-list(longterm=c(7*n_points_per_day,Inf), 
                 diurnal=c(n_points_per_day/4,7*n_points_per_day), 
                 subdiurnal=c(0,n_points_per_day/4))

#window length
wl_max<-round(n_points/3)
wl_vector<-c(wl_max,7*n_points_per_day,n_points_per_day)

#number of components
n.comp_vector<-c(30,30,30)

#2. Run SSA
var_decom_list<-list()

for (i in 1:ncol(df_ssa)){
  cur<-df_ssa %>% select(i)
  cur_decomp<-filterTSeriesSSA(series = cur[,1],
                               borders.wl = bands_list,
                               M = wl_vector,
                               n.comp =n.comp_vector,
                               repeat.extr=c(2,1,1),
                               center.series=TRUE,
                               harmonics = c(0,0,0),
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
write.csv(ts_chosen_var, file.path(out_directory, "ssa_decomposed_time_series.csv"))

#----------------------PLOT AND SAVE RESULTS--------------------------
#1. Plot decomposition results
png(file.path(out_directory, "SSA_decomposition_output.png"), units="px", width=2000, height=2000, res=300)
par(oma = c(4,2,2,2), mfrow = c(3, 1))

plot(ts_chosen_var$doy.dayfract, ts_chosen_var$SIF_A_sfm..mW.m.2nm.1sr.1., type="l", xlab=c("DOY"), ylab=c("SIF A [mW/m2/nm/sr]"), cex.axis=1.5, cex.lab=1.5)
lines(ts_chosen_var$doy.dayfract, ts_chosen_var$SIF_A_sfm..mW.m.2nm.1sr.1.noise_longterm, col="red")
lines(ts_chosen_var$doy.dayfract, ts_chosen_var$SIF_A_sfm..mW.m.2nm.1sr.1.noise_diurnal, col="blue")
lines(ts_chosen_var$doy.dayfract, ts_chosen_var$SIF_A_sfm..mW.m.2nm.1sr.1.noise_subdiurnal, col="darkviolet")

plot(ts_chosen_var$doy.dayfract, ts_chosen_var$SIF_B_sfm..mW.m.2nm.1sr.1., type="l", xlab=c("DOY"), ylab=c("SIF B [mW/m2/nm/sr]"),cex.axis=1.5, cex.lab=1.5)
lines(ts_chosen_var$doy.dayfract, ts_chosen_var$SIF_B_sfm..mW.m.2nm.1sr.1.noise_longterm, col="red")
lines(ts_chosen_var$doy.dayfract, ts_chosen_var$SIF_B_sfm..mW.m.2nm.1sr.1.noise_diurnal, col="blue")
lines(ts_chosen_var$doy.dayfract, ts_chosen_var$SIF_B_sfm..mW.m.2nm.1sr.1.noise_subdiurnal, col="darkviolet")

plot(ts_chosen_var$doy.dayfract, ts_chosen_var$PRI, type="l", xlab=c("DOY"), ylab=c("PRI"),cex.axis=1.5, cex.lab=1.5)
lines(ts_chosen_var$doy.dayfract, ts_chosen_var$PRInoise_longterm, col="red")
lines(ts_chosen_var$doy.dayfract, ts_chosen_var$PRInoise_diurnal, col="blue")
lines(ts_chosen_var$doy.dayfract, ts_chosen_var$PRInoise_subdiurnal, col="darkviolet")

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend=c("Original", "Longterm", "Diurnal", "Sub-diurnal"),
       col=c("black", "red", "blue", "darkviolet"), lty=1, cex=1.5, xpd = TRUE, horiz = TRUE)
while (!is.null(dev.list()))  dev.off()
#------------------------------------------------------------------------------------------
#2. Scatter plots between decomposed variables

#Filter time series to exclude night-time data and interpolated data, and data acquired under SZA > 70

ts_plots<-ts_chosen_var[ts_chosen_var$night_time==0 & ts_chosen_var$no_data==0 & ts_chosen_var$SZA<=70,]

#Plot SIF A vs PRI original 
sifa_pri_orig<-ggplot(data= ts_plots,aes(SIF_A_sfm..mW.m.2nm.1sr.1.,PRI,
                               colour=NDVI))+geom_point()+
  scale_color_viridis()+
  xlab("SIF A [mW/m2/nm/sr]")+ylab("PRI")+
  labs(subtitle = "SZA<=70")+
  theme(plot.title = element_text(size=17,face="bold"))+
  theme(legend.title = element_text(size=17))+
  theme(legend.text = element_text(size=15))+
  theme(plot.subtitle = element_text(size=15))+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=17))

#Plot SIF B vs PRI original 
sifb_pri_orig<-ggplot(data= ts_plots,aes(SIF_B_sfm..mW.m.2nm.1sr.1.,PRI,
                                         colour=NDVI))+geom_point()+
  scale_color_viridis()+
  xlab("SIF B [mW/m2/nm/sr]")+ylab("PRI")+
  labs(subtitle = "SZA<=70")+
  theme(plot.title = element_text(size=17,face="bold"))+
  theme(legend.title = element_text(size=17))+
  theme(legend.text = element_text(size=15))+
  theme(plot.subtitle = element_text(size=15))+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=17))

#Plot SIF A vs SIF B original 
sifa_sifb_orig<-ggplot(data= ts_plots,aes(SIF_A_sfm..mW.m.2nm.1sr.1.,SIF_B_sfm..mW.m.2nm.1sr.1.,
                                         colour=NDVI))+geom_point()+
  scale_color_viridis()+
  xlab("SIF A [mW/m2/nm/sr]")+ylab("SIF B [mW/m2/nm/sr]")+
  labs(subtitle = "SZA<=70")+
  theme(plot.title = element_text(size=17,face="bold"))+
  theme(legend.title = element_text(size=17))+
  theme(legend.text = element_text(size=15))+
  theme(plot.subtitle = element_text(size=15))+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=17))

#Plot SIF A vs PRI sub-diurnal SSA 
sifa_pri_sub_diur_ssa<-ggplot(data= ts_plots,aes(SIF_A_sfm..mW.m.2nm.1sr.1.noise_subdiurnal,
                                          PRInoise_subdiurnal,
                                          colour=NDVI))+geom_point()+
  scale_color_viridis()+
  xlab("SIF A [mW/m2/nm/sr] Sub-Diurnal (SSA)")+ylab("PRI Sub-Diurnal (SSA)")+
  labs(subtitle = "SZA<=70")+
  theme(plot.title = element_text(size=17,face="bold"))+
  theme(legend.title = element_text(size=17))+
  theme(legend.text = element_text(size=15))+
  theme(plot.subtitle = element_text(size=15))+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=17))


#Plot SIF B vs PRI sub-diurnal ssa 
sifb_pri_sub_diur_ssa<-ggplot(data= ts_plots,aes(SIF_B_sfm..mW.m.2nm.1sr.1.noise_subdiurnal,
                                                  PRInoise_subdiurnal,
                                                  colour=NDVI))+geom_point()+
  scale_color_viridis()+
  xlab("SIF B [mW/m2/nm/sr] Sub-Diurnal (SSA)")+ylab("PRI Sub-Diurnal (SSA)")+
  labs(subtitle = "SZA<=70")+
  theme(plot.title = element_text(size=17,face="bold"))+
  theme(legend.title = element_text(size=17))+
  theme(legend.text = element_text(size=15))+
  theme(plot.subtitle = element_text(size=15))+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=17))

#Plot SIF A vs SIF B sub-diurnal ssa 
sifa_sifb_sub_diur_ssa<-ggplot(data= ts_plots,aes(SIF_A_sfm..mW.m.2nm.1sr.1.noise_subdiurnal,
                                                  SIF_B_sfm..mW.m.2nm.1sr.1.noise_subdiurnal,
                                                 colour=NDVI))+geom_point()+
  scale_color_viridis()+
  xlab("SIF A [mW/m2/nm/sr] Sub-Diurnal (SSA)")+ylab("SIF B [mW/m2/nm/sr] Sub-Diurnal (SSA)")+
  labs(subtitle = "SZA<=70")+
  theme(plot.title = element_text(size=17,face="bold"))+
  theme(legend.title = element_text(size=17))+
  theme(legend.text = element_text(size=15))+
  theme(plot.subtitle = element_text(size=15))+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=17))

plot_all<-sifa_sifb_orig+sifa_pri_orig+sifb_pri_orig+
  sifa_sifb_sub_diur_ssa+sifa_pri_sub_diur_ssa+sifb_pri_sub_diur_ssa+
  plot_layout(ncol=3, guides ="collect")

ggsave(plot_all,path =out_directory, file=paste0("SSA_decomposition_scatterplots",".png"), 
       width = 42, height = 25, units = "cm")

