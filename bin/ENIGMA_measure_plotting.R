library("ggplot2", lib.loc="/opt/quarantine/R-extras/3.1.1/build")

## read in the data
setwd("/projects/edickie/code/spins-1yr-analysis")
subjects <- read.table("list.csv", quote="\"")

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


volumes <- merge(subjects,LandRvolumes,by="SubjID")
volumes <- subset(volumes, !is.na(volumes$Dx))
ggplot(volumes, aes(x=site, y=LLatVent, shape=Dx, fill=Dx)) +
  geom_boxplot(outlier.shape = NA) + geom_point(alpha = 0.5, position=position_jitterdodge(jitter.width=0.01))
