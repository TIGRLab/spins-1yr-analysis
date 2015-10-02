
## read in the data
subjects <- read.table("/projects/jdv/data/spins/1yr/list.csv", quote="\"")
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
subjects$site[grep('CMH',as.character(subjects$SubjID))]<-"CMH"
subjects$site[grep('ZHH',as.character(subjects$SubjID))]<-"ZHH"
subjects$site[grep('MRC',as.character(subjects$SubjID))]<-"MRC"
subjects$site <- as.factor(subjects$site)


volumes <- merge(subjects,LandRvolumes,by="SubjID")
