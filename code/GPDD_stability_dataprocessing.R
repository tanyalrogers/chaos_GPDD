# Data processing for GPDD stability analysis

library("rgpdd")
library("ggplot2")
library("dplyr")
library("tidyr")
library("rEDM")
library("purrr")

# Prior studies ####

#Sibly et al 2007 (634)
#graded 2-5, at least 10 obs, some regression constraints

#Knape and Valpine 2012 MainID (627)
#"removing harvest and nonindex based data, data sampled at non-annual intervals 
# and time series taking less than 15 unique values" #no length constraint, median 23.
#    kvMainID=read.csv("Data/KnapeValpineMainID.csv")

#Clark Luis 2019 (640)
#graded 3 or higher, at least 30 obs, taxonomic constraints, constraints on repeating zeros
#at least 5 unique values
#file that was posted
#    gpdd_st=read.csv("Data/Raw TimeSeries Data.csv")

# Processing ####
#list of IDs, life history, and output of gpdd ts used (from Clark Luis 2019)
gpdd_lifehistory=read.csv("Data/lifehistory_output.csv")

#function for longest run of zeros
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
  summarize(#sumpop = sum(PopulationUntransformed, na.rm=T),
    #minpop = min(PopulationUntransformed, na.rm=T), 
    datasetlength=(max(SeriesStep)-min(SeriesStep)+1),
    ndatapoints = length(which(!is.na(PopulationUntransformed))),
    uniquevals = n_distinct(PopulationUntransformed, na.rm = T),
    zerorun = zero_run(.data),
    zeros = length(which(PopulationUntransformed==0)),
    propzeros = zeros/ndatapoints,
    propnonmissing = ndatapoints/datasetlength) %>% 
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
#549
gpdd_filter <- gpdd_join %>%
  filter(Reliability>=2) %>%  #alternatively 3 (517)
  filter(ndatapoints>=30) %>% 
  filter(uniquevals>=5) %>% #alternatively 15
  filter(SamplingProtocol!="Harvest") %>% 
  filter(propzeros<0.6) %>% 
  droplevels()
#minpop and sumpop exclusions drop out after this filtering

#filter data, join timeperiod info, nest by MainID
gpdd_d_nest <- gpdd_data %>% 
  filter(MainID %in% unique(gpdd_filter$MainID)) %>% 
  left_join(select(gpdd_timeperiod, TimePeriodID, TimePeriod, TimePeriodGroup),by="TimePeriodID")  %>% 
  droplevels() %>% 
  group_by(MainID) %>%
  nest()

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

#timescale of observations
gpdd_d$timescale_MinAge=ifelse(gpdd_d$SamplingInterval=="annual", gpdd_d$datasetlength/(gpdd_d$MinAge/12),
                               ifelse(gpdd_d$SamplingInterval=="seasonal", gpdd_d$datasetlength/(gpdd_d$MinAge/6),
                                      ifelse(gpdd_d$SamplingInterval=="bimonthly", gpdd_d$datasetlength/(gpdd_d$MinAge/2),
                                             ifelse(gpdd_d$SamplingInterval=="monthly", gpdd_d$datasetlength/(gpdd_d$MinAge),
                                                    ifelse(gpdd_d$SamplingInterval=="4-week", gpdd_d$datasetlength/(gpdd_d$MinAge*1.07),
                                                           ifelse(gpdd_d$SamplingInterval=="weekly", gpdd_d$datasetlength/(gpdd_d$MinAge*4.286),NA))))))
gpdd_d$timescale_Lifespan=ifelse(gpdd_d$SamplingInterval=="annual", gpdd_d$datasetlength/(gpdd_d$Lifespan/12),
                                 ifelse(gpdd_d$SamplingInterval=="seasonal", gpdd_d$datasetlength/(gpdd_d$Lifespan/6),
                                        ifelse(gpdd_d$SamplingInterval=="bimonthly", gpdd_d$datasetlength/(gpdd_d$Lifespan/2),
                                               ifelse(gpdd_d$SamplingInterval=="monthly", gpdd_d$datasetlength/(gpdd_d$Lifespan),
                                                      ifelse(gpdd_d$SamplingInterval=="4-week", gpdd_d$datasetlength/(gpdd_d$Lifespan*1.07),
                                                             ifelse(gpdd_d$SamplingInterval=="weekly", gpdd_d$datasetlength/(gpdd_d$Lifespan*4.286),NA))))))

#join nested data, location info, life history info to filtered main table
#remove daily dataset
#compute timescale metrics
#subset/reorder columns
#PRIMARY NESTED TABLE
gpdd_d <- gpdd_filter %>%
  left_join(gpdd_d_nest,by="MainID") %>% 
  left_join(select(gpdd_lifehistory, MainID, Mass:TrL),by="MainID") %>% 
  left_join(select(gpdd_location, LocationID, ExactName, Country, Continent, LongDD, LatDD, SpatialAccuracy, LocationExtent),by="LocationID") %>% 
  mutate(timescale_MinAge=ifelse(SamplingInterval=="annual", datasetlength/(MinAge/12),
                                 ifelse(SamplingInterval=="seasonal", datasetlength/(MinAge/6),
                                        ifelse(SamplingInterval=="bimonthly", datasetlength/(MinAge/2),
                                               ifelse(SamplingInterval=="monthly", datasetlength/(MinAge),
                                                      ifelse(SamplingInterval=="4-week", datasetlength/(MinAge*1.07),
                                                             ifelse(SamplingInterval=="weekly", datasetlength/(MinAge*4.286),NA)))))),
         timescale_Lifespan=ifelse(SamplingInterval=="annual", datasetlength/(Lifespan/12),
                                   ifelse(SamplingInterval=="seasonal", datasetlength/(Lifespan/6),
                                          ifelse(SamplingInterval=="bimonthly", datasetlength/(Lifespan/2),
                                                 ifelse(SamplingInterval=="monthly", datasetlength/(Lifespan),
                                                        ifelse(SamplingInterval=="4-week", datasetlength/(Lifespan*1.07),
                                                               ifelse(SamplingInterval=="weekly", datasetlength/(Lifespan*4.286),NA)))))),
         timestep_MinAge=ifelse(SamplingInterval=="annual", 1/(MinAge/12),
                                 ifelse(SamplingInterval=="seasonal", 1/(MinAge/6),
                                        ifelse(SamplingInterval=="bimonthly", 1/(MinAge/2),
                                               ifelse(SamplingInterval=="monthly", 1/(MinAge),
                                                      ifelse(SamplingInterval=="4-week", 1/(MinAge*1.07),
                                                             ifelse(SamplingInterval=="weekly", 1/(MinAge*4.286),NA)))))),
         timestep_Lifespan=ifelse(SamplingInterval=="annual", 1/(Lifespan/12),
                                   ifelse(SamplingInterval=="seasonal", 1/(Lifespan/6),
                                          ifelse(SamplingInterval=="bimonthly", 1/(Lifespan/2),
                                                 ifelse(SamplingInterval=="monthly", 1/(Lifespan),
                                                        ifelse(SamplingInterval=="4-week", 1/(Lifespan*1.07),
                                                               ifelse(SamplingInterval=="weekly", 1/(Lifespan*4.286),NA))))))
  ) %>% 
  select(MainID:LocationID, TaxonName, CommonName, Reliability, SamplingInterval, 
         datasetlength:propnonmissing, TaxonomicPhylum:TaxonomicGenus, 
         TaxonomicLevel, Mass:TrL, timescale_MinAge:timestep_Lifespan, ExactName:LocationExtent, SamplingUnits, SourceTransform, Notes=Notes.x,
         data) %>% 
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

#old rescaling function
# rescale=function(data) {
#   mind=min(data$PopulationUntransformed, na.rm=T)
#   maxd=max(data$PopulationUntransformed, na.rm=T)
#   data$PopRescale=((data$PopulationUntransformed-mind)/(maxd-mind))+1
#   data$PopRescale_fd=c(NA, diff(data$PopRescale))
#   data$PopRescale_log=log(data$PopRescale)
#   data$PopRescale_gr=c(NA, diff(data$PopRescale_log))
#   return(as.data.frame(data))
# }

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

#unnest table
gpdd_dun = unnest(select(gpdd_d, MainID, CommonName, data_rescale))

#files for export
write.csv(gpdd_dun, "./Data/gpdd_timeseries.csv", row.names = F)
write.csv(select(gpdd_d, MainID:Notes, data_rescale_case), "./Data/gpdd_ts_metadata.csv", row.names = F)

save(gpdd_d, file = "./Data/gpdd_d_v2.Rdata")

#function for body masses
cvol=function(mm, ratio=1/3.5){mm/10*pi*(mm/10/2*ratio)^2}
