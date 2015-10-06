library("ggplot2", lib.loc="/opt/quarantine/R-extras/3.1.1/build")
library("reshape2", lib.loc="/opt/quarantine/R-extras/3.1.1/build")

## read in the data
setwd("/projects/edickie/code/spins-1yr-analysis/bin")

subjects <- read.table("list.csv", header=TRUE, quote="\"")

## freesurfer derived measures
LandRvolumes <- read.csv("/projects/edickie/analysis/freesurfer_QA/SPINS/LandRvolumes.csv")
CorticalMeasuresENIGMA_ThickAvg <- read.csv("/projects/edickie/analysis/freesurfer_QA/SPINS/CorticalMeasuresENIGMA_ThickAvg.csv")
CorticalMeasuresENIGMA_SurfAvg <- read.csv("/projects/edickie/analysis/freesurfer_QA/SPINS/CorticalMeasuresENIGMA_SurfAvg.csv")

names(subjects)<-c("subid","code")
subjects$Dx[as.character(subjects$code)=="p"]<-"Patient"
subjects$Dx[subjects$code=="c"]<-"Control"
subjects$Dx[subjects$code=="h"]<-"Control"
subjects$Dx[subjects$code=="x"]<-NA
subjects$Dx <- factor(subjects$Dx)

subjects$SubjID <- paste0(subjects$subid,'_01')
subjects$site <- NA
subjects$site[grep('CMH',as.character(subjects$SubjID))]<-"CMH"
subjects$site[grep('ZHH',as.character(subjects$SubjID))]<-"ZHH"
subjects$site[grep('MRC',as.character(subjects$SubjID))]<-"MRC"
subjects$site <- as.factor(subjects$site)

subjects <- subjects[,c("SubjID","Dx","site")]

volumes <- merge(subjects,LandRvolumes,by="SubjID")
volumes <- subset(volumes, !is.na(volumes$Dx))
volume.ROIs <-names(LandRvolumes)[2:(length(LandRvolumes)-2)]
### run a linear model in each site to residulized for ICV and Site
#take the orig matrix and set the cols from colnames to NA before we start
# get the start
df <-volumes
colnames <- volume.ROIs

residualize_ICV <- function(df,colnames) {
  ## define a new dataframe
  results_df<-df
  results_df[ ,c(colnames)] <- NA
  ### run many models  
  ICV_lm <- lapply(colnames, function(x) {
    form <-as.formula(paste0(x, '~ 1 + ICV'))
    lm(form, data=df)
  })
  ## output the residuals into a dataframe
  for (i in 1:length(colnames)){
    results_df[ ,c(colnames[i])] <- ICV_lm[[i]]$residuals
  } 
  ## return the results
  return(results_df)
}

## run residualize_ICV separately for each site them recombine 
res_CMH <- residualize_ICV(subset(volumes, site=="CMH"),volume.ROIs)
res_ZHH <- residualize_ICV(subset(volumes, site=="ZHH"),volume.ROIs)
res_MRC <- residualize_ICV(subset(volumes, site=="MRC"),volume.ROIs)
volumes_resICV <- rbind(res_CMH,res_MRC,res_ZHH)

ggplot(volumes_resICV, aes(x=site, y=LLatVent, shape=Dx, fill=Dx)) +
  geom_boxplot(outlier.shape = NA) + geom_point(alpha = 0.5, position=position_jitterdodge(jitter.width=0.01))

ggplot(volumes, aes(x=site, y=ICV, shape=Dx, fill=Dx)) +
  geom_boxplot(outlier.shape = NA) + geom_point(alpha = 0.5, position=position_jitterdodge(jitter.width=0.01))



volumesICV.melted <- melt(volumes_resICV[ ,c(names(subjects),volume.ROIs)],
                       id.vars=names(subjects),
                       measure.vars=volume.ROIs,
                       variable.name="ROI",
                       value.name="Volume")

volumesICV.melted$Hemisphere <- NA
volumesICV.melted$Hemisphere[substring(volumesICV.melted$ROI,1,1)=="L"] <- "Left"
volumesICV.melted$Hemisphere[substring(volumesICV.melted$ROI,1,1)=="R"] <- "Right"

volumesICV.melted$Region <- NA
volumesICV.melted$Region[grepl('LatVent',volumesICV.melted$ROI)]<-'Lateral Ventricle'
volumesICV.melted$Region[grepl('thal',volumesICV.melted$ROI)]<-'Thalamus'
volumesICV.melted$Region[grepl('caud',volumesICV.melted$ROI)]<-'Caudate'
volumesICV.melted$Region[grepl('put',volumesICV.melted$ROI)]<-'Putamen'
volumesICV.melted$Region[grepl('pal',volumesICV.melted$ROI)]<-'Palidum'
volumesICV.melted$Region[grepl('hippo',volumesICV.melted$ROI)]<-'Hippocampus'
volumesICV.melted$Region[grepl('amyg',volumesICV.melted$ROI)]<-'Amygdala'
volumesICV.melted$Region[grepl('accumb',volumesICV.melted$ROI)]<-'Accumbens'


ggplot(subset(volumesICV.melted,Region != 'Lateral Ventricle'), aes(x=Dx, y=Volume, fill=Dx)) +
  geom_boxplot(outlier.shape = NA) + geom_point(alpha = 0.5, position=position_jitterdodge(jitter.width=0.01),aes(colour=site)) +
  facet_grid(Region ~ Hemisphere)

ggplot(subset(volumesICV.melted,Region == 'Lateral Ventricle'), aes(x=Dx, y=Volume)) +
  geom_boxplot(outlier.shape = NA) + geom_point(aes(colour=site)) +
  facet_grid(Region ~ Hemisphere)

ggplot(subset(volumesICV.melted,Region == 'Amygdala' | Region == 'Hippocampus'), aes(x=Dx, y=Volume)) +
  geom_boxplot(outlier.shape = NA) + geom_point(aes(colour=site)) +
  facet_grid(Region ~ Hemisphere)

ggplot(subset(volumesICV.melted,Region == 'Caudate' | Region == 'Putamen' | Region == 'Palidum' | Region =='Accumbens'), aes(x=Dx, y=Volume)) +
  geom_boxplot(outlier.shape = NA) + geom_point(aes(colour=site)) +
  facet_grid(Region ~ Hemisphere)


