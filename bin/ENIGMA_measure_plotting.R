library("ggplot2", lib.loc="/opt/quarantine/R-extras/3.1.1/build")
library("reshape2", lib.loc="/opt/quarantine/R-extras/3.1.1/build")
library("effsize", lib.loc="~/R/x86_64-unknown-linux-gnu-library/3.1")

## read in the data
setwd("/projects/edickie/code/spins-1yr-analysis/freesurfer")

subjects <- read.table("list.csv", header=TRUE, quote="\"")

## freesurfer derived measures
LandRvolumes <- read.csv("LandRvolumes.csv")
CorticalMeasuresENIGMA_ThickAvg <- read.csv("CorticalMeasuresENIGMA_ThickAvg.csv")
CorticalMeasuresENIGMA_SurfAvg <- read.csv("CorticalMeasuresENIGMA_SurfAvg.csv")

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

# t-tests http://stats.stackexchange.com/questions/128894/p-value-correction-for-multiple-t-tests
## note this is expecting your data to already be melted so that splitvar is the old column names...
ttestpvals <- function(df, splitvar, ttestformula) {  
  myformula = formula(ttestformula)
  p.ttest = do.call("rbind", 
                    lapply(split(df, df[ ,c(splitvar)]), 
                           function(x) t.test(myformula,x)$p.value))
  p.ttest = cbind(p.ttest, 
                  p.adjust(p.ttest, method="fdr"), 
                  p.adjust(p.ttest, method="bonferroni"))
  colnames(p.ttest) = c("raw", "fdr", "bonferroni")
  p.ttest = as.data.frame(p.ttest)
  p.ttest$fdrsig <- ''
  p.ttest$fdrsig[p.ttest$fdr < 0.05] <-'*'
  return(p.ttest)
}

# independant 
run_many_lm <- function(df,splitvar1,splitvar2,lmformula) {
  p.value = c()
  beta = c()
  n = c()
  svar1 = c()
  svar2 = c()
  myformula <- formula(lmformula)
  index = 1
  for (site_x in levels(df[ ,c(splitvar1)])) {
    site_df = subset(df, df[ ,c(splitvar1)] == site_x)
    for (ROI_x in levels(df[ ,c(splitvar2)])) {
      roi_df = subset(site_df, site_df[,c(splitvar2)] == ROI_x)
      mod = lm(myformula, roi_df)
      p.value[index] = anova(mod)[1,5]
      beta[index] = mod$coef[2]
      n[index] = nrow(roi_df)  # missing values counted, get this from model? 
      svar1[index] = site_x
      svar2[index] = ROI_x
      index = index + 1
    }
  }
  df.ind = data.frame(x1 = as.factor(svar1), 
                      x2 = as.factor(svar2), 
                      beta = beta, 
                      p.value = p.value, 
                      n = n)
  names(df.ind) <- c(names(df[ ,c(splitvar1,splitvar2)]),'beta','p.value','n')
  df.ind$p0.5 = ifelse(df.ind$p.value < 0.05, T, F)
  return(df.ind)
}


## run residualize_ICV separately for each site them recombine 
vols_CMH <- subset(volumes, site=="CMH")
vols_ZHH <- subset(volumes, site=="ZHH")
vols_MRC <- subset(volumes, site=="MRC")
res_CMH <- residualize_ICV(vols_CMH,volume.ROIs)
res_ZHH <- residualize_ICV(vols_ZHH,volume.ROIs)
res_MRC <- residualize_ICV(vols_MRC,volume.ROIs)

## rbind the pooled results
volumes_resICV <- rbind(res_CMH,res_MRC,res_ZHH)




# ## do effect size calculations
# yvars =  
# coD <- cohen.d(res_CMH$LLatVent,res_CMH$Dx,hedges.correction=TRUE)
# coD$estimate
# coD$conf.int

## box plot thing of the ICV
ggplot(volumes, aes(x=site, y=ICV, shape=Dx, fill=Dx)) +
  geom_boxplot(outlier.shape = NA) + geom_point(alpha = 0.5, position=position_jitterdodge(jitter.width=0.01))


### now do the melting for both
volumesICV.melted <- melt(volumes_resICV[ ,c(names(subjects),volume.ROIs)],
                          id.vars=names(subjects),
                          measure.vars=volume.ROIs,
                          variable.name="ROI",
                          value.name="Volume_ICVresid")

volumes.melted <- melt(volumes[ ,c(names(subjects),volume.ROIs)],
                          id.vars=names(subjects),
                          measure.vars=volume.ROIs,
                          variable.name="ROI",
                          value.name="Volume")

### adding pretty labels for graphs
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


### run the defined function of the data
results_volumes_resICV <- run_many_lm(volumesICV.melted,'site','ROI','Volume ~ Dx')
ggplot(results_volumes_resICV, aes(x = ROI, 
                   y = beta, 
                   size = -log(p.value, base =10), 
                   shape = p0.5, 
                   color = as.factor(site))) + 
  geom_point() + 
  geom_hline(aes(yintercept=0)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("beta coefficients and p-values for independent tests across site/ROI for Volumes (ICV residualized)")


results_volumes_raw <- run_many_lm(volumes.melted,'site','ROI','Volume ~ Dx')
ggplot(results_volumes_raw, aes(x = ROI, 
                                   y = beta, 
                                   size = -log(p.value, base =10), 
                                   shape = p0.5, 
                                   color = as.factor(site))) + 
  geom_point() + 
  geom_hline(aes(yintercept=0)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("beta coefficients and p-values for independent tests across site/ROI for Volumes (raw)")

################################## now for thinkness
thickness <- merge(subjects,CorticalMeasuresENIGMA_ThickAvg,by="SubjID")
thickness <- subset(thickness, !is.na(thickness$Dx))
thickness.ROIs <-names(CorticalMeasuresENIGMA_ThickAvg)[2:(length(CorticalMeasuresENIGMA_ThickAvg)-2)]
thickness.ROIs.nopooled <- thickness.ROIs[1:68]
thickness.ROIs.pooled <- thickness.ROIs[69:70]


## run residualize_ICV separately for each site them recombine 
thickness_CMH <- subset(thickness, site=="CMH")
thickness_ZHH <- subset(thickness, site=="ZHH")
thickness_MRC <- subset(thickness, site=="MRC")
res_CMH <- residualize_ICV(thickness_CMH,thickness.ROIs)
res_ZHH <- residualize_ICV(thickness_ZHH,thickness.ROIs)
res_MRC <- residualize_ICV(thickness_MRC,thickness.ROIs)

## rbind the pooled results
thickness_resICV <- rbind(res_CMH,res_MRC,res_ZHH)

### now do the melting for both
thicknessICV.melted <- melt(thickness_resICV[ ,c(names(subjects),thickness.ROIs)],
                          id.vars=names(subjects),
                          measure.vars=thickness.ROIs,
                          variable.name="ROI",
                          value.name="Thickness_ICVresid")

thickness.melted <- melt(thickness[ ,c(names(subjects),thickness.ROIs)],
                       id.vars=names(subjects),
                       measure.vars=thickness.ROIs,
                       variable.name="ROI",
                       value.name="Thickness")

### run the defined function of the data
results_thickness_resICV <- run_many_lm(thicknessICV.melted,'site','ROI','Thickness_ICVresid ~ Dx')
results_thickness_resICV_nopooled <-subset(results_thickness_resICV, 
                                           ROI != "LThickness" & ROI != "RThickness" & ROI != "LSurfArea" & ROI != "RSurfArea")
ggplot(results_thickness_resICV_nopooled
                , aes(x = ROI, 
                      y = beta, 
                      size = -log(p.value, base =10), 
                      shape = p0.5, 
                      color = as.factor(site))) + 
  geom_point() + 
  geom_hline(aes(yintercept=0)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("beta coefficients and p-values for independent tests across site/ROI for Thickness (resid)")


results_thickness_raw <- run_many_lm(thickness.melted,'site','ROI','Thickness ~ Dx')
results_thickness_raw_nopooled <-subset(results_thickness_raw, 
                                           ROI != "LThickness" & ROI != "RThickness" & ROI != "LSurfArea" & ROI != "RSurfArea")
ggplot(results_thickness_raw_nopooled, aes(x = ROI, 
                                y = beta, 
                                size = -log(p.value, base =10), 
                                shape = p0.5, 
                                color = as.factor(site))) + 
  geom_point() + 
  geom_hline(aes(yintercept=0)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("beta coefficients and p-values for independent tests across site/ROI for Thickness (raw)")


################################## now for Surface Area
surfarea <- merge(subjects,CorticalMeasuresENIGMA_SurfAvg,by="SubjID")
surfarea <- subset(surfarea, !is.na(surfarea$Dx))
surfarea.ROIs <-names(CorticalMeasuresENIGMA_SurfAvg)[2:(length(CorticalMeasuresENIGMA_SurfAvg)-2)]
surfarea.ROIs.nopooled <- surfarea.ROIs[1:68]
surface.ROIs.pooled <- surfarea.ROIs[71:72]
## run residualize_ICV separately for each site them recombine 
surfarea_CMH <- subset(surfarea, site=="CMH")
surfarea_ZHH <- subset(surfarea, site=="ZHH")
surfarea_MRC <- subset(surfarea, site=="MRC")
res_CMH <- residualize_ICV(surfarea_CMH,surfarea.ROIs)
res_ZHH <- residualize_ICV(surfarea_ZHH,surfarea.ROIs)
res_MRC <- residualize_ICV(surfarea_MRC,surfarea.ROIs)

## rbind the pooled results
surfarea_resICV <- rbind(res_CMH,res_MRC,res_ZHH)

### now do the melting for both
surfareaICV.melted <- melt(surfarea_resICV[ ,c(names(subjects),surfarea.ROIs)],
                            id.vars=names(subjects),
                            measure.vars=surfarea.ROIs,
                            variable.name="ROI",
                            value.name="surfarea_ICVresid")

surfarea.melted <- melt(surfarea[ ,c(names(subjects),surfarea.ROIs)],
                         id.vars=names(subjects),
                         measure.vars=surfarea.ROIs,
                         variable.name="ROI",
                         value.name="surfarea")

### run the defined function of the data
results_surfarea_resICV <- run_many_lm(surfareaICV.melted,'site','ROI','surfarea_ICVresid ~ Dx')
results_surfarea_resICV_nopooled <-subset(results_surfarea_resICV, 
                                       ROI != "LThickness" & ROI != "RThickness" & ROI != "LSurfArea" & ROI != "RSurfArea")
ggplot(results_surfarea_resICV, aes(x = ROI, 
                                     y = beta, 
                                     size = -log(p.value, base =10), 
                                     shape = p0.5, 
                                     color = as.factor(site))) + 
  geom_point() + 
  geom_hline(aes(yintercept=0)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("beta coefficients and p-values for independent tests across site/ROI for Surface Area (resid)")


results_surfarea_raw <- run_many_lm(surfarea.melted,'site','ROI','surfarea ~ Dx')
results_surfarea_raw_nopooled <-subset(results_surfarea_raw, 
                                       ROI != "LThickness" & ROI != "RThickness" & ROI != "LSurfArea" & ROI != "RSurfArea")
ggplot(results_surfarea_raw_nopooled, aes(x = ROI, y = beta, size = -log(p.value, base =10), shape = p0.5, color = as.factor(site))) + 
  geom_point() + 
  geom_hline(aes(yintercept=0)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("beta coefficients and p-values for independent tests across site/ROI for Surface Area (raw)")
