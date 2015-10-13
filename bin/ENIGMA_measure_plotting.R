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
subjects <- subset(subjects, !is.na(subjects$Dx))


### run a linear model in each site to residulized for ICV and Site
#take the orig matrix and set the cols from colnames to NA before we start
# get the start

### my old residualize - works on data in wide format - not going to use it here
# residualize_ICV <- function(df,colnames) {
#   ## define a new dataframe
#   results_df<-df
#   results_df[ ,c(colnames)] <- NA
#   ### run many models  
#   ICV_lm <- lapply(colnames, function(x) {
#     form <-as.formula(paste0(x, '~ 1 + ICV'))
#     lm(form, data=df)
#   })
#   ## output the residuals into a dataframe
#   for (i in 1:length(colnames)){
#     results_df[ ,c(colnames[i])] <- ICV_lm[[i]]$residuals
#   } 
#   ## return the results
#   return(results_df)
# }

# t-tests http://stats.stackexchange.com/questions/128894/p-value-correction-for-multiple-t-tests
## note this is expecting your data to already be melted so that splitvar is the old column names...
# df=volumessurf.melted
# splitvar = 'ROI'
# ttestformula = 'Volume ~ Dx'
# ## do effect size calculations
# yvars =  
# coD <- cohen.d(res_CMH$LLatVent,res_CMH$Dx,hedges.correction=TRUE)
# coD$estimate
# coD$conf.int
ttestpvals <- function(df, splitvar, ttestformula) {  
  myformula = formula(ttestformula)
  coD = do.call("rbind", 
                    lapply(split(df, df[ ,c(splitvar)]), 
                           function(x) cohen.d(myformula,x,hedges.correction=TRUE)$estimate))
  coDSE = do.call("rbind", 
                lapply(split(df, df[ ,c(splitvar)]), 
                       function(x) cohen.d(myformula,x,hedges.correction=TRUE)$var))
  t.ttest = do.call("rbind", 
                    lapply(split(df, df[ ,c(splitvar)]), 
                           function(x) t.test(myformula,x)$statistic))
  df.ttest = do.call("rbind",
                    lapply(split(df, df[ ,c(splitvar)]), 
                           function(x) t.test(myformula,x)$parameter))
  p.ttest = do.call("rbind",
                    lapply(split(df, df[ ,c(splitvar)]), 
                           function(x) t.test(myformula,x)$p.value))
  n.ttest = do.call("rbind", lapply(split(df, df[ ,c(splitvar)]), nrow))
  p.ttest = cbind(coD, coDSE,t.ttest, df.ttest, p.ttest, 
                  p.adjust(p.ttest, method="fdr"), 
                  p.adjust(p.ttest, method="bonferroni"),
                  n.ttest)
  colnames(p.ttest) = c('cohen.D','cohenD.se','tstat','df',"raw", "fdr", "bonferroni","n")
  p.ttest = as.data.frame(p.ttest)
  p.ttest$p0.05 = ifelse(p.ttest$fdr < 0.05, T, F)
  return(p.ttest)
}

# independant 
run_many_lm <- function(df,splitvar1,splitvar2,lmformula) {
  p.value = c()
  beta = c()
  coD = c()
  n = c()
  svar1 = c()
  svar2 = c()
  myformula <- formula(lmformula)
  index = 1
  for (site_x in levels(df[ ,c(splitvar1)])) {
    site_df = subset(df, df[ ,c(splitvar1)] == site_x)
    for (ROI_x in levels(df[ ,c(splitvar2)])) {
      roi_df = subset(site_df, site_df[,c(splitvar2)] == ROI_x)
      coD[index] = cohen.d(myformula,roi_df,hedges.correction=TRUE)$estimate ## now estimates cohen's D for ttest model
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
                      cohen.D = coD,
                      beta = beta, 
                      p.value = p.value, 
                      n = n)
  names(df.ind) <- c(names(df[ ,c(splitvar1,splitvar2)]),'cohen.D','beta','p.value','n')
  df.ind$p0.05 = ifelse(df.ind$p.value < 0.05, T, F)
  return(df.ind)
}

# independant 
# df=volumes_surf.melted
# splitvar = "ROI"
# lmformula = 'Volume ~ site'

run_many_sitelm_byROI <- function(df,splitvar,lmformula) {
  p.value = c()
  Fvalue = c()
  n = c()
  svar = c()
  myformula <- formula(lmformula)
  index = 1
  for (ROI_x in levels(df[ ,c(splitvar)])) {
    roi_df = subset(df, df[,c(splitvar)] == ROI_x)
    mod = lm(myformula, roi_df)
    p.value[index] = anova(mod)[1,5]
    Fvalue[index] = anova(mod)[1,4]
    n[index] = nrow(roi_df)  # missing values counted, get this from model? 
    svar[index] = ROI_x
    index = index + 1
  }
  df.ind = data.frame(x1 = as.factor(svar), 
                      Fvalue = Fvalue, 
                      p.value = p.value, 
                      n = n)
  names(df.ind) <- c(names(df[ ,c(splitvar)]),'Fstat','p.value','n')
  df.ind$p0.05 = ifelse(df.ind$p.value < 0.05, T, F)
  return(df.ind)
}


# residuals of FA ~ site

residualize.nested <- function(df,splitvar1,splitvar2,lmformula,residvarname) {
  myformula = formula(lmformula)
  mycolnames <- c(names(df),residvarname)
  results <- df
  results$residuals <- NA
  results <- NULL
  for (site_x in levels(df[ ,c(splitvar1)])) {
    site_df = subset(df, df[ ,c(splitvar1)] == site_x)
    for (ROI_x in levels(df[ ,c(splitvar2)])) {
      roi_df = subset(site_df, site_df[,c(splitvar2)] == ROI_x)      
      roi_df$residuals =  resid(lm(myformula, roi_df)) 
      names(roi_df) <- mycolnames
      results <- rbind(results,roi_df) 
    }
  }
  names(results) <- mycolnames
  return(results)
}

### function that runs run_many_lm, then plots and saves
# dataframe.melted  = volumessurf.melted
# yvarname = 'Volume'
run_nestedDxlm_and_plot <- function(dataframe.melted,yvarname, saveplot=TRUE, pdfwidth=9,pdfheight=7) {
  results_df <- run_many_lm(dataframe.melted,'site','ROI',paste(yvarname,'~ Dx'))
  thisplot <- ggplot(results_df, aes(x = ROI, 
                         y = cohen.D, 
                         size = -log(p.value, base =10), 
                         shape = p0.05, 
                         color = as.factor(site))) + 
           geom_point() + 
           geom_hline(aes(yintercept=0)) + 
           theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
           ggtitle(paste("beta coefficients and p-values for independent tests across site/ROI for",yvarname))
  plot(thisplot)
  if (saveplot==TRUE) {ggsave(paste(yvarname,'nestedDxlmresults.pdf',sep='_'),thisplot,width=pdfwidth,height=pdfheight)}
  return(results_df)
}

### function to pool L and R roi's 
# this need the data in wide format - do this before melting
# df = volumesplussurf
# idvars = names(subjects)
# roivars = setdiff(names(df),names(subjects))

sumLeftandRight <- function(df,idvars,roivars) {
  LeftNames <-grep('L',roivars,value=T)
  results = df
  pooledvars = c()
  for (lname in LeftNames) {
    if (substr(lname,1,1) == 'L') {
      pooledname = substr(lname,2,nchar(lname))
      rname = paste0('R',pooledname)
      results$pooledvar = results[,c(lname)]+resutls[,c(rname)]
      names(results)[names(results)=='pooledvar'] <- pooledname
      pooledvars = c(pooledvars, pooledname) 
    }
  }
  results <- results[ ,c(idvars, pooledvars)]
  return(results)
}



################################## now for thinkness
thickness <- merge(subjects,CorticalMeasuresENIGMA_ThickAvg,by="SubjID")
thickness <- subset(thickness, !is.na(thickness$Dx))
thickness.ROIs <-names(CorticalMeasuresENIGMA_ThickAvg)[2:(length(CorticalMeasuresENIGMA_ThickAvg)-2)]
thickness.ROIs.notpooled <- thickness.ROIs[1:68]
thickness.ROIs.pooled <- thickness.ROIs[69:70]

### now do the melting
thickness.melted <- melt(thickness[ ,c(names(subjects),thickness.ROIs,'ICV')],
                       id.vars=c(names(subjects),thickness.ROIs.pooled,'ICV'),
                       measure.vars=thickness.ROIs.notpooled,
                       variable.name="ROI",
                       value.name="Thickness")

##resid for ICV
thickness.melted <- residualize.nested(thickness.melted,"site","ROI",'Thickness ~ 1 + ICV','Thickness_ICVresid') 

##resid for total thickness
thickness.melted$TotalThickness <- thickness.melted$LThickness + thickness.melted$RThickness
thickness.melted <- residualize.nested(thickness.melted,"site","ROI",'Thickness ~ 1 + TotalThickness','Thickness_Totalresid') 

## run linear models on everybody - with wite, within ROI
results_thickness_raw <- run_nestedDxlm_and_plot(thickness.melted,"Thickness", save=F, pdfwidth=20,pdfheight=8)
results_thickness_residicv <- run_nestedDxlm_and_plot(thickness.melted,"Thickness_ICVresid", save=F, pdfwidth=20,pdfheight=8)
results_thickness_residtotal <- run_nestedDxlm_and_plot(thickness.melted,"Thickness_Totalresid", save=F, pdfwidth=20,pdfheight=8)

# run ttests on everybody  - pooling across sites 
ttestres_thickness_raw = ttestpvals(thickness.melted,'ROI','Thickness ~ Dx')
ttestres_thickness_residicv = ttestpvals(thickness.melted,'ROI','Thickness_ICVresid ~ Dx')                                           
ttestres_thickness_residtotal = ttestpvals(thickness.melted,'ROI','Thickness_Totalresid ~ Dx') 

# write.csv(ttestres_thickness_raw,'ttestres_thickness_raw.csv')
# write.csv(ttestres_thickness_residicv,'ttestres_thickness_residicv.csv')
# write.csv(ttestres_thickness_residtotal,'ttestres_thickness_residtotal.csv')


################################## now for Surface Area
surfarea <- merge(subjects,CorticalMeasuresENIGMA_SurfAvg,by="SubjID")
surfarea <- subset(surfarea, !is.na(surfarea$Dx))
surfarea.ROIs <-names(CorticalMeasuresENIGMA_SurfAvg)[2:(length(CorticalMeasuresENIGMA_SurfAvg)-2)]
surfarea.ROIs.notpooled <- surfarea.ROIs[1:68]
surfarea.ROIs.pooled <- surfarea.ROIs[71:72]

### now do the melting
surfarea.melted <- melt(surfarea[ ,c(names(subjects),surfarea.ROIs,'ICV')],
                         id.vars=c(names(subjects),surfarea.ROIs.pooled,'ICV'),
                         measure.vars=surfarea.ROIs.notpooled,
                         variable.name="ROI",
                         value.name="SurfArea")

##resid for ICV
surfarea.melted <- residualize.nested(surfarea.melted,"site","ROI",'SurfArea ~ 1 + ICV','SurfArea_ICVresid') 

##resid for total surfarea
surfarea.melted$TotalSurfArea <- surfarea.melted$LSurfArea + surfarea.melted$RSurfArea
surfarea.melted <- residualize.nested(surfarea.melted,"site","ROI",'SurfArea ~ 1 + TotalSurfArea','SurfArea_Totalresid') 

## run linear models on everybody - with wite, within ROI
results_surfarea_raw <- run_nestedDxlm_and_plot(surfarea.melted,"SurfArea", pdfwidth=20,pdfheight=8)
results_surfarea_residicv <- run_nestedDxlm_and_plot(surfarea.melted,"SurfArea_ICVresid", pdfwidth=20,pdfheight=8)
results_surfarea_residtotal <- run_nestedDxlm_and_plot(surfarea.melted,"SurfArea_Totalresid", pdfwidth=20,pdfheight=8)

# run ttests on everybody  - pooling across sites 
ttestres_surfarea_raw = ttestpvals(surfarea.melted,'ROI','SurfArea ~ Dx')
ttestres_surfarea_residicv = ttestpvals(surfarea.melted,'ROI','SurfArea_ICVresid ~ Dx')                                           
ttestres_surfarea_residtotal = ttestpvals(surfarea.melted,'ROI','SurfArea_Totalresid ~ Dx') 

# write.csv(ttestres_surfarea_raw,'ttestres_surfarea_raw.csv')
# write.csv(ttestres_surfarea_residicv,'ttestres_surfarea_residicv.csv')
# write.csv(ttestres_surfarea_residtotal,'ttestres_surfarea_residtotal.csv')

#################the other stuff - volumes and total Thickness and Surface Area

volumes <- merge(subjects,LandRvolumes,by="SubjID")
volumes <- subset(volumes, !is.na(volumes$Dx))
volume.ROIs <-names(LandRvolumes)[2:(length(LandRvolumes)-2)]

## calculate some stuff to add
thickness$calcMeanThickness <- apply(thickness[,c(thickness.ROIs.notpooled)],1,mean)
thickness$TotalThickness <- thickness$RThickness + thickness$LThickness
surfarea$calcMeanSurfArea <- apply(surfarea[,c(surfarea.ROIs.notpooled)],1,mean)
surfarea$TotalSurfArea <- surfarea$RSurfArea + surfarea$LSurfArea

## merge the volumes with some surface into
volumesplussurf1 <- merge(volumes, 
                         thickness[ ,c("SubjID",thickness.ROIs.pooled,'calcMeanThickness','TotalThickness')], 
                         by="SubjID")
volumesplussurf <- merge(volumesplussurf1, 
                         surfarea[ ,c("SubjID",surfarea.ROIs.pooled,'calcMeanSurfArea','TotalSurfArea')], 
                         by="SubjID")
rm(volumesplussurf1)



## box plot thing of the ICV
ggplot(volumes, aes(x=site, y=ICV, shape=Dx, fill=Dx)) +
  geom_boxplot(outlier.shape = NA) + geom_point(alpha = 0.5, position=position_jitterdodge(jitter.width=0.01))


### now do the melting for both
volumessurf_ICVresid.melted <- melt(volumesplussurf,
                       id.vars=c(names(subjects),'ICV'),
                       measure.vars=c(volume.ROIs,
                                      thickness.ROIs.pooled,
                                      surfarea.ROIs.pooled,
                                      'calcMeanThickness','TotalThickness',
                                      'calcMeanSurfArea','TotalSurfArea'),
                       variable.name="ROI",
                       value.name="Volume")

volumessurf.melted <- melt(volumesplussurf,
                           id.vars=c(names(subjects)),
                           measure.vars=c(volume.ROIs,
                                          thickness.ROIs.pooled,
                                          surfarea.ROIs.pooled,
                                          'ICV',
                                          'calcMeanThickness','TotalThickness',
                                          'calcMeanSurfArea','TotalSurfArea'),
                           variable.name="ROI",
                           value.name="Volume")

### adding pretty labels for graphs
volumessurf.melted$Hemisphere <- NA
volumessurf.melted$Hemisphere[substring(volumessurf.melted$ROI,1,1)=="L"] <- "Left"
volumessurf.melted$Hemisphere[substring(volumessurf.melted$ROI,1,1)=="R"] <- "Right"

volumessurf.melted$Region <- NA
volumessurf.melted$Region[grepl('LatVent',volumessurf.melted$ROI)]<-'Lateral Ventricle'
volumessurf.melted$Region[grepl('thal',volumessurf.melted$ROI)]<-'Thalamus'
volumessurf.melted$Region[grepl('caud',volumessurf.melted$ROI)]<-'Caudate'
volumessurf.melted$Region[grepl('put',volumessurf.melted$ROI)]<-'Putamen'
volumessurf.melted$Region[grepl('pal',volumessurf.melted$ROI)]<-'Palidum'
volumessurf.melted$Region[grepl('hippo',volumessurf.melted$ROI)]<-'Hippocampus'
volumessurf.melted$Region[grepl('amyg',volumessurf.melted$ROI)]<-'Amygdala'
volumessurf.melted$Region[grepl('accumb',volumessurf.melted$ROI)]<-'Accumbens'

volumessurf_ICVresid.melted <- residualize.nested(volumessurf_ICVresid.melted,"site","ROI",'Volume ~ 1 + ICV','Volume_ICVresid') 

## run linear models on everybody - with wite, within ROI
results_volumes_raw <- run_nestedDxlm_and_plot(volumessurf.melted,"Volume", pdfwidth=9,pdfheight=7)
results_volumes_residicv <- run_nestedDxlm_and_plot(volumessurf.melted,"Volume_ICVresid", pdfwidth=9,pdfheight=7)

# run ttests on everybody  - pooling across sites 
ttestres_volumes_raw = ttestpvals(volumessurf.melted,'ROI','Volume ~ Dx')
ttestres_volumes_residICV2 = ttestpvals(volumessurf.melted,'ROI','Volume ~ Dx')
# ###plot only surface area and thickness stuff
# thick_raw_pooled <- subset(thickness.melted, 
#                                ROI == "LThickness" | ROI == "RThickness")
# ## box plot thing of the ICV
# ggplot(thick_raw_pooled, aes(x=Dx, y=Thickness)) +
#   geom_boxplot(outlier.shape = NA) + 
#   geom_jitter(position = position_jitter(width = .1), size=3, aes(color=site)) +
#   facet_wrap(~ROI) + labs(y="Thickness (raw)")
# ggsave('RawAllThickness_boxplot.png',width=6,height=5)
# 
# thick_ICV_pooled <- subset(thicknessICV.melted, 
#                            ROI == "LThickness" | ROI == "RThickness")
# ## box plot thing of the ICV
# ggplot(thick_ICV_pooled, aes(x=Dx, y=Thickness_ICVresid)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(position = position_jitter(width = .1), size=3, aes(color=site)) +
#   facet_wrap(~ROI) + labs(y="Thickness (resid. for ICV)")
# ggsave('ResidAllThickness_boxplot.png',width=6,height=5)
# 
# ##### for surface area
# sa_raw_pooled <- subset(surfarea.melted, 
#                            ROI == "LSurfArea" | ROI == "RSurfArea")
# ## box plot thing of the ICV
# ggplot(sa_raw_pooled, aes(x=Dx, y=surfarea)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(position = position_jitter(width = .1), size=3, aes(color=site)) +
#   facet_wrap(~ROI) + labs(y="Surface Area (raw)")
# ggsave('RawAllSurfaceArea_boxplot.png',width=6,height=5)
# 
# sa_ICV_pooled <- subset(surfareaICV.melted, 
#                         ROI == "LSurfArea" | ROI == "RSurfArea")
# ## box plot thing of the ICV
# ggplot(sa_ICV_pooled, aes(x=Dx, y=surfarea_ICVresid)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(position = position_jitter(width = .1), size=3, aes(color=site)) +
#   facet_wrap(~ROI) + labs(y="Surface Area (residualized for ICV)")
# ggsave('ResidAllSurfaceArea_boxplot.png',width=6,height=5)

## poole the Left and Right data
pooledVolumes <- poolLeftandRight(volumesplussurf,                  #df = volumes left and right
                                  c(names(subjects),'ICV'),         #idvars = names(subjects) AND 'ICV
                                  setdiff(names(df),names(subjects)))  #ROIvars is everything else..
pooledvars <- setdiff(names(pooledVolumes),names(subjects))

## 
pooledVolumes.melted <- melt(pooledVolumes,
                           id.vars=c(names(subjects)),
                           measure.vars=pooledvars,
                           variable.name="ROI",
                           value.name="Volume")

pooledVolumes.melted$Region <- NA
pooledVolumes.melted$Region[grepl('LatVent',pooledVolumes.melted$ROI)]<-'Lateral Ventricle'
pooledVolumes.melted$Region[grepl('thal',pooledVolumes.melted$ROI)]<-'Thalamus'
pooledVolumes.melted$Region[grepl('caud',pooledVolumes.melted$ROI)]<-'Caudate'
pooledVolumes.melted$Region[grepl('put',pooledVolumes.melted$ROI)]<-'Putamen'
pooledVolumes.melted$Region[grepl('pal',pooledVolumes.melted$ROI)]<-'Palidum'
pooledVolumes.melted$Region[grepl('hippo',pooledVolumes.melted$ROI)]<-'Hippocampus'
pooledVolumes.melted$Region[grepl('amyg',pooledVolumes.melted$ROI)]<-'Amygdala'
pooledVolumes.melted$Region[grepl('accumb',pooledVolumes.melted$ROI)]<-'Accumbens'
pooledVolumes.melted$Region[grepl('ICV',pooledVolumes.melted$ROI)]<-'ICV'

results_pvolumes_raw <- run_nestedDxlm_and_plot(pooledVolumes.melted,"Volume", save=F, pdfwidth=9,pdfheight=7)
ttestres_pvolumes_raw = ttestpvals(pooledVolumes.melted,'ROI','Volume ~ Dx')
ttestres_pvolumes_raw$site = 'pooled'
### do separately in all site and recombine
ttestres_pvolumes_raw_CMH = ttestpvals(subset(pooledVolumes.melted,site=='CMH'),'ROI','Volume ~ Dx')
ttestres_pvolumes_raw_CMH$site = 'CMH'
ttestres_pvolumes_raw_ZHH = ttestpvals(subset(pooledVolumes.melted,site=='ZHH'),'ROI','Volume ~ Dx')
ttestres_pvolumes_raw_ZHH$site = 'ZHH'
ttestres_pvolumes_raw_MRC = ttestpvals(subset(pooledVolumes.melted,site=='MRC'),'ROI','Volume ~ Dx')
ttestres_pvolumes_raw_MRC$site = 'MRC'

alleffects = rbind(ttestres_pvolumes_raw,ttestres_pvolumes_raw_CMH, ttestres_pvolumes_raw_ZHH, ttestres_pvolumes_raw_MRC)
alleffects$Region <- NA
alleffects$Region[grepl('LatVent',row.names(alleffects))]<-'Lateral Ventricle'
alleffects$Region[grepl('thal',row.names(alleffects))]<-'Thalamus'
alleffects$Region[grepl('caud',row.names(alleffects))]<-'Caudate'
alleffects$Region[grepl('put',row.names(alleffects))]<-'Putamen'
alleffects$Region[grepl('pal',row.names(alleffects))]<-'Palidum'
alleffects$Region[grepl('hippo',row.names(alleffects))]<-'Hippocampus'
alleffects$Region[grepl('amyg',row.names(alleffects))]<-'Amygdala'
alleffects$Region[grepl('accumb',row.names(alleffects))]<-'Accumbens'
alleffects$Region[grepl('ICV',row.names(alleffects))]<-'ICV'
alleffects$Region[grepl('Thickness',row.names(alleffects))]<-'Cortical Thickness'
alleffects$Region[grepl('SurfArea',row.names(alleffects))]<-'Cortical Surface Area'

alleffects$Region <- factor(alleffects$Region, 
                            levels = c('Cortical Thickness', 'Cortical Surface Area',
                                       'Hippocampus', 'Amygdala', 'Thalamus', 'Accumbens',
                                        'ICV', 'Caudate', 'Putamen', 'Palidum', 'Lateral Ventricle')[11:1])


ggplot(subset(alleffects,site=='pooled'), aes(y=cohen.D, x=Region)) +
  geom_errorbar(aes(ymin = (cohen.D - cohenD.se), ymax =  (cohen.D + cohenD.se)), width = 0.2) +
  geom_point(size = 8, shape = 18) + 
  geom_point(data=subset(alleffects,site!='pooled'), aes(y=cohen.D, x=Region, colour = site), size = 3) +
  geom_hline(yintercept=0) +
  scale_y_reverse() + 
  labs(x='',y="Group Effect Size (Hedges g)") +
  coord_flip()

ggsave("ENIGMAvolumes_horizpplot.pdf", width=9, height=7)

###################### Now for ENIGMA dti measures
FAresults <- read.csv("../dti/data/enigmaDTI-FA-results.csv")
FAvars <- names(FAresults)[2:length(names(FAresults))]
FA_average_col <- grep('Average', FAvars, value = T)
FA_LandRrois <- c(grep('.L_', FAvars, value = T),grep('.R_', FAvars, value = T))
PooledRois <- setdiff(FAvars,c(FA_LandRrois,FA_average_col))

ENIGMA_look_up_table <- read.delim("/projects/edickie/code/hcp_extras/templates/ENIGMA_look_up_table.txt", header=FALSE)
names(ENIGMA_look_up_table) <-c('key','name','X1','description')

FAvars_noFA <-gsub('_FA','',FAvars, fixed=T)
FAvars_noFA <-gsub('.','-',FAvars_noFA, fixed=T)

Skelvars <- intersect(FAvars_noFA,as.character(ENIGMA_look_up_table$name))
Skelvars_FAed <- paste0(gsub('-','.',Skelvars, fixed=T),'_FA')

## merge down the one's from the list
FAresults <- merge(subjects, FAresults, by.x = "SubjID", by.y = "id")

##melt
FAskel.melted <- melt(FAresults[,c(names(subjects),Skelvars_FAed)],
                             id.vars=c(names(subjects)),
                             measure.vars=Skelvars_FAed,
                             variable.name="ROI",
                             value.name="FA")

FAskel.melted <- residualize.nested(FAskel.melted,"site","ROI",'FA ~ 1','FA_resid') 
ttestres_FA_pooled = ttestpvals(FAskel.melted,'ROI','FA ~ Dx')
ttestres_FAresid_pooled = ttestpvals(FAskel.melted,'ROI','FA_resid ~ Dx')
ttestres_FAresid_pooled$ROI = gsub('_FA','',row.names(ttestres_FAresid_pooled), fixed=T)
ttestres_FAresid_pooled$ROI = gsub('.','-',ttestres_FAresid_pooled$ROI, fixed=T)
ttestres_FAresid_pooled <- ttestres_FAresid_pooled[ ,c("ROI",names(ttestres_FAresid_pooled)[1:9])]
write.csv(ttestres_FAresid_pooled,"ttestres_FAresid_allsites.csv",row.names = F)
