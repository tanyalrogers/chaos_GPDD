# Processes, filters, and cleans GPDD data
# Tanya Rogers

library("rgpdd")
library("dplyr")
library("tidyr")
library("purrr")

#Filtering used by prior studies

#Sibly et al 2007 (634)
#graded 2-5, at least 10 obs, some regression constraints

#Knape and Valpine 2012 MainID (627)
#"removing harvest and nonindex based data, data sampled at non-annual intervals 
#and time series taking less than 15 unique values" #no length constraint, median 23.

#Clark Luis 2019 (640)
#graded 3 or higher, at least 30 obs, taxonomic constraints, constraints on repeating zeros
#at least 5 unique values

#load life history info (includes just focal series)
gpdd_lifehistory=read.csv("data/gpdd_lifehistory.csv")

#function to get longest run of zeros
zero_run <- function(x){
  k<-rle(x$PopulationUntransformed)
  if(is.element(0, k$values)){
    m<-k$lengths[k$values==0]
    return(max(m))
  } else {
    return(0)
  }
}

#calculate filtering metrics
gpdd_temp <- gpdd_main %>%
  inner_join(gpdd_data,by="MainID") %>% 
  group_by(MainID) %>%
  filter(SeriesStep!=-9999) %>% 
  summarize(
    datasetlength=max((max(SeriesStep)-min(SeriesStep)+1),length(PopulationUntransformed)),
    ndatapoints = length(which(!is.na(PopulationUntransformed))),
    uniquevals = n_distinct(PopulationUntransformed, na.rm = T),
    zerorun = zero_run(.data),
    zeros = length(which(PopulationUntransformed==0)),
    propzeros = zeros/ndatapoints,
    missingvals = datasetlength-ndatapoints,
    propmissing = 1-(ndatapoints/datasetlength)) %>% 
  ungroup() 

#join filtering metrics and taxon info to main table
#remove some columns
gpdd_join <- gpdd_main %>%
  left_join(gpdd_temp,by="MainID") %>% 
  inner_join(gpdd_taxon,by="TaxonID") %>%
  select(-SiblyFittedTheta, -SiblyThetaCILower, -SiblyThetaCIUpper,
         -SiblyExtremeNEffect, -SiblyReturnRate, 
         -WoldaCode, -Authority, - Notes.y)

#filter time series
gpdd_filter <- gpdd_join %>%
  filter(Reliability>=2) %>%  
  filter(ndatapoints>=30) %>% 
  filter(uniquevals>=5) %>% 
  filter(SamplingProtocol!="Harvest") %>% 
  filter(propzeros<0.6) %>% 
  filter(propmissing<=0.22) %>% #missing data filter
  filter(MainID!=2774) %>% #exclude shorter of duplicate series
  filter(MainID!=1870) %>% #exclude shorter of duplicate series
  filter(MainID!=9308) %>% #exclude lower quality of duplicate series
  filter(MainID!=9685) %>% #exclude disease series
  filter(MainID!=9686) %>% #exclude disease series
  droplevels()

#filter data, join timeperiod info, nest by MainID
gpdd_d_nest <- gpdd_data %>% 
  filter(MainID %in% unique(gpdd_filter$MainID)) %>% 
  left_join(select(gpdd_timeperiod, TimePeriodID, TimePeriod, TimePeriodGroup),by="TimePeriodID")  %>% 
  droplevels() %>% 
  group_by(MainID) %>%
  nest()

#fix error in ts 9833
gpdd_d_nest[gpdd_d_nest$MainID==9833,]$data[[1]]$SeriesStep[94]<-93

#calculate unique time period IDs
gpdd_d_nest$unTimePeriodID=map_dbl(gpdd_d_nest$data, ~ n_distinct(.x$TimePeriod)) 

#get correct sampling intervals
sampling_interval=function(unTimePeriodID, data) {
  if(unTimePeriodID==1) {return("annual")}
  if(unTimePeriodID==2) {return("seasonal")}
  if(unTimePeriodID==4) {return("monthly")}
  if(unTimePeriodID==6) {return("bimonthly")}
  if(unTimePeriodID==7) {return("4-week")}
  if(unTimePeriodID==9) {
    if(any(data$TimePeriodGroup=="month", na.rm=T)) {return("monthly")}
    else {return("4-week")}
  }
  if(unTimePeriodID==12) {
    if(any(data$TimePeriodGroup=="week", na.rm=T)) {return("weekly")}
    else {return("monthly")}
  }
  if(unTimePeriodID>200) {return("daily")}
}
gpdd_d_nest$SamplingInterval=map2_chr(gpdd_d_nest$unTimePeriodID, gpdd_d_nest$data, sampling_interval) 

#function to compute timescale ratios
timescale_ratio=function(SamplingInterval, num, den_mo) {
  ifelse(SamplingInterval=="annual", num/(den_mo/12),
         ifelse(SamplingInterval=="seasonal", num/(den_mo/6),
                ifelse(SamplingInterval=="bimonthly", num/(den_mo/2),
                       ifelse(SamplingInterval=="monthly", num/(den_mo),
                              ifelse(SamplingInterval=="4-week", num/(den_mo*1.07),
                                     ifelse(SamplingInterval=="weekly", num/(den_mo*4.286),NA))))))
  
}

#identify taxon classes, all others are categorized as zooplankton
focal_taxa=c("Aves","Osteichthyes", "Mammalia", "Bacillariophyceae", "Dinophyceae", "Insecta")

#PRIMARY NESTED TABLE
gpdd_d <- gpdd_filter %>%
  #join nested data, location info, life history info to filtered main table
  left_join(gpdd_d_nest,by="MainID") %>% 
  left_join(select(gpdd_lifehistory, MainID, Biome, Mass_g:TrL),by="MainID") %>% 
  left_join(select(gpdd_location, LocationID, ExactName, Country, Continent, LongDD, LatDD),by="LocationID") %>% 
  #relabel taxon classes
  #compute timescale metrics
  mutate(TaxonomicClass2=ifelse(TaxonomicClass %in% focal_taxa, as.character(TaxonomicClass), "Zooplankton"),
         TaxonomicClass3=recode(TaxonomicClass2,Aves="Birds",Osteichthyes="Bony fishes",Mammalia="Mammals",Bacillariophyceae="Phytoplankton",Dinophyceae="Phytoplankton",Insecta="Insects"),
         timescale_MinAge=timescale_ratio(SamplingInterval, datasetlength, MinAge_mo), #generations sampled
         timescale_Lifespan=timescale_ratio(SamplingInterval, datasetlength, Lifespan_mo), #lifespans sampled
         timestep_MinAge=timescale_ratio(SamplingInterval, 1, MinAge_mo), #generations per timestep
         timestep_Lifespan=timescale_ratio(SamplingInterval, 1, Lifespan_mo)) %>% #lifespans per timestep
  #subset/reorder columns
  select(MainID:LocationID, TaxonName, CommonName, Reliability, SamplingInterval, 
         datasetlength:propmissing, TaxonomicPhylum:TaxonomicGenus, TaxonomicClass2, TaxonomicClass3, 
         TaxonomicLevel, Mass_g:TrL, timescale_MinAge:timestep_Lifespan, ExactName:LatDD, Biome, SamplingUnits, SourceTransform, Notes=Notes.x,
         data) %>% 
  #remove daily dataset
  filter(SamplingInterval!="daily") %>% 
  droplevels()

#fill in missing timepoints and drop columns
filltimepoints=function(data) {
  data2=filter(data, SeriesStep!=-9999)
  allsteps=data.frame(SeriesStep=min(data2$SeriesStep, na.rm=T):max(data2$SeriesStep, na.rm=T))
  dfill=full_join(data, allsteps, by = "SeriesStep") %>% arrange(SeriesStep)
  dfill=filter(dfill, SeriesStep!=-9999) %>% 
    select(SeriesStep, SampleYear, Population, PopulationUntransformed)
  return(as.data.frame(dfill))
}
gpdd_d$data_fill=map(gpdd_d$data, filltimepoints) 

#rescale abundance
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
rescale=function(data, diag=F) {
  mind=min(data$PopulationUntransformed, na.rm=T)
  maxd=max(data$PopulationUntransformed, na.rm=T)
  if(mind>0) {
    data$PopRescale=data$PopulationUntransformed/sd(data$PopulationUntransformed, na.rm=T)
    case=1 #no zeros, nothing added
  } else {
    if(all(is.wholenumber(data$PopulationUntransformed), na.rm=T)) {
      data$PopRescale=(data$PopulationUntransformed+1)/sd(data$PopulationUntransformed, na.rm=T)
      case=2 #zeros, integers, plus 1
    } else {
      minnz=min(data$PopulationUntransformed[data$PopulationUntransformed>0], na.rm=T)
      data$PopRescale=(data$PopulationUntransformed+minnz)/sd(data$PopulationUntransformed, na.rm=T)
      case=3 #zeros, non-integers, plus min non-zero value
    }
  }
  data$PopRescale_fd=c(NA, diff(data$PopRescale))
  data$PopRescale_log=log(data$PopRescale)
  data$PopRescale_gr=c(NA, diff(data$PopRescale_log))
  if(diag) return(case) 
  else return(as.data.frame(data))
}
gpdd_d$data_rescale=map(gpdd_d$data_fill, rescale)
gpdd_d$data_rescale_case=map_dbl(gpdd_d$data_fill, rescale, diag=T)
#table(gpdd_d$data_rescalediag)

#quantify monotonic trend
monotonic_eval=function(data) {
  (cor(data$SeriesStep, data$PopRescale, use="p", method="spearman"))^2
}
gpdd_d$monotonicR2=map_dbl(gpdd_d$data_rescale, monotonic_eval)

#unnest table
gpdd_dun = unnest(select(gpdd_d, MainID, CommonName, data_rescale))

#export files
write.csv(gpdd_dun, "./data/gpdd_timeseries.csv", row.names = F)
write.csv(select(gpdd_d, MainID:Notes, data_rescale_case, monotonicR2), "./data/gpdd_ts_metadata.csv", row.names = F)

# #save nested table
# save(gpdd_d, file = "./data/gpdd_d.Rdata")

# #plot a time series
# plotMainID=function(ID) {
#   testplot=filter(gpdd_d, MainID==ID)
#   plot(testplot$data_rescale[[1]]$SeriesStep, testplot$data_rescale[[1]]$PopRescale, 
#        ylab="PopRescale", xlab="SeriesStep", main=paste(ID, testplot$CommonName))
#   lines(testplot$data_rescale[[1]]$SeriesStep, testplot$data_rescale[[1]]$PopRescale)
# }
# plotMainID(9953)
# plotMainID(56)
