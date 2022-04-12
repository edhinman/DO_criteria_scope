library(wqTools)
library(tidyr)
library(dplyr)
library(leaflet)
### read ATTAINS and find river/stream AUs that are 3A
# au_uses = read.csv("UTAHDWQ_uses_9706.csv")
# aus_3A = subset(au_uses, au_uses$USE_NAME=="Aquatic Wildlife (Cold Water)"&au_uses$WATER_TYPE_NAME=="RIVER")

# Read in sites and assign beneficial uses, then filter to data at 3A sites
# https://www.waterqualitydata.us/#statecode=US%3A49&siteType=Stream&startDateLo=10-01-2014&startDateHi=09-30-2020&sampleMedia=Water&mimeType=csv&sorted=no&providers=NWIS&providers=STEWARDS&providers=STORET
sites = read.csv("station.csv")
sites_aus = assignUses(sites, lat = "LatitudeMeasure",long= "LongitudeMeasure", flatten=TRUE)
sites_3A = subset(sites_aus, sites_aus$BeneficialUse=="3A")

#https://www.waterqualitydata.us/#statecode=US%3A49&siteType=Stream&sampleMedia=Water&characteristicName=Dissolved%20oxygen%20(DO)&characteristicName=Temperature%2C%20water&characteristicName=Specific%20conductance&startDateLo=10-01-2014&startDateHi=09-30-2020&minactivities=10&mimeType=csv&dataProfile=narrowResult&providers=NWIS&providers=STEWARDS&providers=STORET
dat = read.csv("narrowresult.csv")
dat_3A = subset(dat, dat$MonitoringLocationIdentifier%in%sites_3A$MonitoringLocationIdentifier)

dat_3A_wide = dat_3A%>%
  group_by(ActivityStartDate,MonitoringLocationIdentifier,CharacteristicName,ResultMeasure.MeasureUnitCode)%>%
  summarise(IR_Value = mean(ResultMeasureValue))%>%
  pivot_wider(id_cols = c(ActivityStartDate,MonitoringLocationIdentifier),names_from = CharacteristicName, values_from = IR_Value)

dat_comp = subset(dat_3A_wide,!is.na(dat_3A_wide$`Dissolved oxygen (DO)`)&!is.na(dat_3A_wide$`Specific conductance`)&!is.na(dat_3A_wide$`Temperature, water`))

dat1 = merge(dat_comp, sites, all.x = TRUE)

lat_long = unique(dat1[,c("LatitudeMeasure","LongitudeMeasure","MonitoringLocationIdentifier")])
write.csv(lat_long,"elevation_sites.csv", row.names = FALSE)
# # Get elevation for each site - doesn't work if not admin
# library(elevatr, sp)
# lat_long1 = sp::SpatialPointsDataFrame(coords=lat_long[,c("LatitudeMeasure","LongitudeMeasure")],lat_long)
# lat_long2 = get_elev_point(lat_long[,c("LatitudeMeasure","LongitudeMeasure")],prj="EPSG:4326",src="epqs")

# Used USGS site to assign elevations because GIS wasn't working either. https://apps.nationalmap.gov/elevation/##bottom4
site_elev = read.csv("UT_Statewide_30m_DEM/elevation_sites_m.csv")
site_elev = merge(site_elev, sites_3A, all.x = TRUE)

# equation for calculating barometric pressure in Pascals
baro_p = function(x,temp=20){
  tempK = temp+273.15
  p0 = 1 # atm at sea level
  g = 9.80665 # m/s^2
  M = 0.0289644 # kg/mol
  R = 8.3144598 # J/(mol*K)
  ph = p0*exp((-g*M*x)/(R*tempK))
  return(ph)
}

site_elev$Barometric_Pressure_atm = baro_p(site_elev$Elev_m)

dat2 = merge(dat1, site_elev, all.x = TRUE)
dat2$Temperature_K = dat2$`Temperature, water`+273.15

# Use Benson and Krause equation to determine DO concentration at 100% saturation
DO_baseline <- function(tempK){
  out = exp(-139.34411+(1.575701*10^5/tempK)-(6.642308*10^7/(tempK^2))+(1.243800*10^10/(tempK^3))-(8.621949*10^11/(tempK^4)))
  return(out)
}

Salinity <- function(SC){
  out = 5.572*10^(-4)*SC+2.02*10^(-9)*SC^2
  return(out)
}

Fs <- function(S,TK){
  out = exp(-S*(0.017674-(10.754/TK)+(2140.7/(TK^2))))
  return(out)
}

Fp <- function(P,t,TK){
  theta = 0.000975-1.426*10^(-5)*t+6.436*10^(-8)*t^2
  u = exp(11.8571-(3840.70/TK)-(216961/(TK^2)))
  out = ((P-u)*(1-theta*P))/((1-u)*(1-theta))
  return(out)
}

dat2$DO_baseline = DO_baseline(dat2$Temperature_K)
dat2$Salinity = Salinity(dat2$`Specific conductance`)
dat2$Fs = Fs(dat2$Salinity,dat2$Temperature_K)
dat2$Fp = Fp(dat2$Barometric_Pressure_atm,dat2$`Temperature, water`,dat2$Temperature_K)
dat2$DO_100sat = dat2$DO_baseline*dat2$Fs*dat2$Fp
dat2$DO_90sat = dat2$DO_100sat*0.9
dat2$Crit_7day = 9.5
dat2$Crit_min = 8
dat2$Crit_7day_sat = ifelse(dat2$DO_100sat<1.1*9.5,dat2$DO_90sat,9.5)
dat2$Crit_min_sat = ifelse(dat2$DO_100sat<1.1*8,dat2$DO_90sat,8)

dat3 = dat2[,c("MonitoringLocationIdentifier","ActivityStartDate","Dissolved oxygen (DO)","Crit_7day","Crit_min","Crit_7day_sat","Crit_min_sat")]
dat3_long = dat3%>%pivot_longer(cols=c("Crit_7day","Crit_min","Crit_7day_sat","Crit_min_sat"),names_to = "Criterion_Name",values_to = "Criterion")
dat3_long$Exceeds = ifelse(dat3_long$`Dissolved oxygen (DO)`<dat3_long$Criterion,1,0)

dat4 = dat3_long%>%group_by(MonitoringLocationIdentifier,Criterion_Name)%>%summarise(Ncount = length(Exceeds), Num_Exceeds = sum(Exceeds), Perc_Exc = sum(Exceeds)/length(Exceeds)*100)
dat4 = subset(dat4, dat4$Ncount>9)
dat4$Impaired = ifelse(dat4$Perc_Exc>10,1,0)
dat4$Type = ifelse(grepl("sat",dat4$Criterion_Name), "Saturation-based","Original concentration")

dat5 = dat4%>%group_by(MonitoringLocationIdentifier,Type)%>%summarise(Impaired = ifelse(any(Impaired==1),"Yes","No"))%>%pivot_wider(id_cols = MonitoringLocationIdentifier,names_from = Type, values_from = Impaired)
dat5 = merge(dat5, unique(dat4[,c("MonitoringLocationIdentifier","Ncount")]), all.x=TRUE)

dat5 = merge(dat5, site_elev, all.x=TRUE)
dat5$color = "purple"
dat5$color[dat5$`Original concentration`=="Yes"&dat5$`Saturation-based`=="Yes"] = "red"
dat5$color[dat5$`Original concentration`=="No"&dat5$`Saturation-based`=="No"] = "blue"
dat5$note = "Impaired with original concentration criteria, not impaired using saturation-based criteria"
dat5$note[dat5$`Original concentration`=="Yes"&dat5$`Saturation-based`=="Yes"] = "Impaired using both types of criteria"
dat5$note[dat5$`Original concentration`=="No"&dat5$`Saturation-based`=="No"] = "Not impaired"


m = wqTools::baseMap()%>%
  addCircleMarkers(data=dat5, lat=dat5$LatitudeMeasure, lng=dat5$LongitudeMeasure, color=dat5$color,
                   popup = paste0(
                     "Monitoring Location ID: ",dat5$MonitoringLocationIdentifier,
                     "<br> Monitoring Location Name: ",dat5$MonitoringLocationName,
                     "<br> Elevation: ",dat5$Elev_m,
                     "<br> Ncount: ",dat5$Ncount,
                     "<br> Description: ",dat5$note
                   ))
