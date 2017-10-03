# summarizing characteristics of EPA NLA 2007 lakes 
# and Northern Highland Lakes District lakes
# 11-17-16
rm(list=ls())

# load data from Renee Brooks
# Brooks et al. 2014. Stable isotope estimates of evaporation:inflow and water 
#residence time for lakes across the United States as a tool for national lake 
#water quality assessments.Limnology and Oceanography, 59: 2150-2165.
# Brooks et al. data available at: 
brooks=read.csv("Brooksetal.csv",header=TRUE,stringsAsFactors=FALSE)
brooks=brooks[brooks$VISIT_NO==1,]

# load data from Read et al.
# Read et al. 2014. The importance of lake-specific characteristics for water quality
# across the continental United States. Ecological Applications, 25: 943-955.
# **NOTE** Residence Time measures in this dataset are derived from McDonald et al.
# McDonald et al. 2013. Inorganic carbon loading as a primary driver of dissolved carbon
# dioxide concentrations in the lakes and reservoirs of the contiguous United States.
# Global Biogeochemical Cycles, 27: 285-295.
# Read et al. data available at: https://portal.lternet.edu/nis/mapbrowse?scope=knb-lter-ntl&identifier=10000&revision=1
read=read.csv("Read_et_al_2014_NLA_data.csv",header=TRUE,stringsAsFactors=FALSE)
read=read[!duplicated(read$SITE_ID),]

# merge two NLA datasets
brooks=brooks[brooks$SITE_ID%in%read$SITE_ID,]

brooks=brooks[order(brooks$SITE_ID),]
read=read[order(read$SITE_ID),]

dim(brooks)
dim(read)
sum(brooks$SITE_ID==read$SITE_ID)

NLA=cbind(brooks,read[,c(30,36,161,181)])

# calculate mean depth according to McDonald et al.
# volume = cone with lake area and max observed depth
# zbar = volume/area
NLA$zbar=(NLA$AREA_HA*1e4*NLA$DEPTHMAX/3)/(NLA$AREA_HA*1e4)


# load data from Hanson et al.
# Hanson et al. 2007. Small lakes dominate a random sample of regional lake characteristics.
# Freshwater Biology, 52: 814-822.
# Hanson et al. data available at: https://lter.limnology.wisc.edu/dataset/fluxes-project-north-temperate-lakes-lter-random-lake-survey-2004
hansonLakes=read.csv("random_lake_survey_lakes.csv",header=TRUE,stringsAsFactors=FALSE)
hansonParams=read.csv("random_lake_survey_measured_parameters.csv",header=TRUE,stringsAsFactors=FALSE)
NHLD=merge(hansonParams,hansonLakes,)

# calculate catchment-to-lake area ratio (drainage ratio)
NHLD$DR=NHLD$watershed_area/NHLD$area

# calculate quantiles
quantile(NLA$AREA_HA/100,probs=c(0,0.025,0.25,0.5,0.75,0.975,1))  #km2
quantile(NLA$BasinAreaToLakeArea_HA,probs=c(0,0.025,0.25,0.5,0.75,0.975,1))
quantile(NLA$zbar,probs=c(0,0.025,0.25,0.5,0.75,0.975,1),na.rm=TRUE)  #m
quantile(NLA$RT*365,probs=c(0,0.025,0.25,0.5,0.75,0.975,1)) #days
quantile(NLA$ResTime,probs=c(0,0.025,0.25,0.5,0.75,0.975,1),na.rm=TRUE) #days
quantile(NLA$E_I,probs=c(0,0.025,0.25,0.5,0.75,0.975,1))

quantile(NHLD$area/1e6,probs=c(0,0.025,0.25,0.5,0.75,0.975,1),na.rm=TRUE) #km2
quantile(NHLD$watershed_area/NHLD$area,probs=c(0,0.025,0.25,0.5,0.75,0.975,1))

