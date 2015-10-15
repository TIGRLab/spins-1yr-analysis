library("ggplot2", lib.loc="/opt/quarantine/R-extras/3.1.1/build")
library("reshape2", lib.loc="/opt/quarantine/R-extras/3.1.1/build")
library("effsize", lib.loc="~/R/x86_64-unknown-linux-gnu-library/3.1")

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
  yvar = gsub(' ','',(strsplit(ttestformula,'~')[[1]][1]))
  xvar = gsub(' ','',(strsplit(ttestformula,'~')[[1]][2]))
  coD = do.call("rbind", 
                lapply(split(df, df[ ,c(splitvar)]), 
                       function(x) cohen.d(myformula,x,hedges.correction=TRUE)))
  t.ttest = do.call("rbind", 
                    lapply(split(df, df[ ,c(splitvar)]), 
                           function(x) t.test(myformula,x)))
  n.ttest = do.call("rbind", lapply(split(df, df[ ,c(splitvar)]), nrow))
  p.ttest = cbind(coD[ ,c('estimate','var')],
                  t.ttest[ ,c('statistic', 'parameter', 'p.value')], 
                  p.adjust(t.ttest[ ,c('p.value')], method="fdr"), 
                  p.adjust(t.ttest[ ,c('p.value')], method="bonferroni"),
                  n.ttest)
  mycolnames = c('cohenD','cohenD.se','tstat','df',"raw", "fdr", "bonferroni","n")
  colnames(p.ttest) = mycolnames
  p.ttest = as.data.frame(p.ttest) 
  p.ttest$p0.05 = ifelse(p.ttest$fdr < 0.05, T, F)
  mycolnames <- c(mycolnames,'p0.05')
  p.ttest$cohenD = as.numeric(p.ttest$cohenD)
  p.ttest$cohenD.se = as.numeric(p.ttest$cohenD.se)
  
  ## now iterate through and find mean and se for subgroups
  for (iname in levels(df[ ,c(xvar)])) {
    h_df = subset(df,df[ ,c(xvar)]==iname)
    means = do.call("rbind", lapply(split(h_df[ ,c(yvar)], h_df[ ,c(splitvar)]), mean))
    SE = do.call("rbind", lapply(split(h_df[ ,c(yvar)], h_df[ ,c(splitvar)]), function(x) sd(x)/sqrt(length(x)) ))
    p.ttest = cbind(p.ttest, means, SE)
    p.ttest = as.data.frame(p.ttest)
    p.ttest$means = as.numeric(p.ttest$means)
    p.ttest$SE = as.numeric(p.ttest$SE)
    mycolnames = c(mycolnames,paste0(c('mean_','se_'),iname))
  }
  colnames(p.ttest) = mycolnames

  
  return(p.ttest)
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

###################### Now for ENIGMA dti measures
## read in the data
setwd("/projects/edickie/code/spins-1yr-analysis/dti")
FAresults <- read.csv("data/enigmaDTI-FA-results.csv")
FAvars <- names(FAresults)[2:length(names(FAresults))]
FA_average_col <- grep('Average', FAvars, value = T)
FA_LandRrois <- c(grep('.L_', FAvars, value = T),grep('.R_', FAvars, value = T))
PooledRois <- setdiff(FAvars,c(FA_LandRrois,FA_average_col,'CC_FA'))


## Ran previously...
# ENIGMA_look_up_table <- read.delim("/projects/edickie/code/hcp_extras/templates/ENIGMA_look_up_table.txt", header=FALSE)
# names(ENIGMA_look_up_table) <-c('KEY','NAME','X1','description')
# engimalut <- ENIGMA_look_up_table[ ,c('KEY',"NAME",'description')]
#write.csv(engimalut, "ENIGMA_look_up_table.csv", row.names = F)
## load the engima lookup table
ENIGMA_look_up_table <- read.csv("skeltemplates/ENIGMA_look_up_table.csv")


FAvars_noFA <-gsub('_FA','',FAvars, fixed=T)
FAvars_noFA <-gsub('.','-',FAvars_noFA, fixed=T)

Skelvars <- intersect(FAvars_noFA,as.character(ENIGMA_look_up_table$NAME))
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
#write.csv(ttestres_FAresid_pooled,"ttestres_FAresid_allsites.csv",row.names = F)

FApooled <- FAresults[ ,c(names(subjects),PooledRois)] 
FApooled_healthyavg <-cbind(apply(subset(FApooled, Dx == "Control")[ ,c(PooledRois)],2,mean),PooledRois)
FApooled_healthyavg <- as.data.frame(FApooled_healthyavg)
PooledRois_horder <- as.character(FApooled_healthyavg[order(FApooled_healthyavg$V1),'PooledRois'])

##melt
FApooled.melted <- melt(FAresults[,c(names(subjects),PooledRois)],
                        id.vars=c(names(subjects)),
                        measure.vars=PooledRois,
                        variable.name="ROI",
                        value.name="FA")

FApooled.melted <- residualize.nested(FApooled.melted,"site","ROI",'FA ~ 1','FA_resid') 

ttestres_pFA_raw = ttestpvals(FApooled.melted,'ROI','FA ~ Dx')
ttestres_pFA_raw$site = 'pooled'
### do separately in all site and recombine
ttestres_pFA_raw_CMH = ttestpvals(subset(FApooled.melted,site=='CMH'),'ROI','FA ~ Dx')
ttestres_pFA_raw_CMH$site = 'CMH'
ttestres_pFA_raw_ZHH = ttestpvals(subset(FApooled.melted,site=='ZHH'),'ROI','FA ~ Dx')
ttestres_pFA_raw_ZHH$site = 'ZHH'
ttestres_pFA_raw_MRC = ttestpvals(subset(FApooled.melted,site=='MRC'),'ROI','FA ~ Dx')
ttestres_pFA_raw_MRC$site = 'MRC'

alleffects_FA = rbind(ttestres_pFA_raw,ttestres_pFA_raw_CMH, ttestres_pFA_raw_ZHH, ttestres_pFA_raw_MRC)

alleffects_FA$Tract <- NA
for (i in 1:nrow(alleffects_FA)) {
  alleffects_FA$Tract[i] <- as.character(strsplit(row.names(alleffects_FA)[i],'_')[[1]][1])
}


tmp <- subset(alleffects_FA,site=='pooled')
alleffects_FA_order <- as.character(tmp[order(as.numeric(tmp$cohenD)),'Tract'])
alleffects_FA$Tract <- factor(alleffects_FA$Tract,levels = alleffects_FA_order)



ggplot(subset(alleffects_FA,site=='pooled'), aes(y=as.numeric(cohenD), x=Tract)) +
  geom_errorbar(aes(ymin = (cohenD - cohenD.se), ymax =  (cohenD + cohenD.se)), width = 0.2) +
  geom_point(size = 8, shape = 18) + 
  geom_point(data=subset(alleffects_FA,site!='pooled'), aes(y=cohenD, x=Tract, colour = site), size = 3) +
  geom_hline(yintercept=0) +
  scale_y_reverse() + 
  labs(x='',y="Group Effect Size (Hedges g)") +
  coord_flip()

alleffects_FA_horder <- gsub('_FA','',PooledRois_horder,fixed=T)
alleffects_FA$Tract <- factor(alleffects_FA$Tract,levels = alleffects_FA_horder)

ggplot(data=subset(alleffects_FA,site!='pooled'),aes(x=Tract, y=mean_Control, color=site, group = site)) +
  geom_errorbar(aes(ymin = mean_Control - se_Control, ymax = mean_Control + se_Control), width = 0.2, size=1.25) +
  geom_line(size=1.25) +
  labs(y = "FA Mean(SE) in Healthy Controls", color = "Site")

ggsave("FAmeaninControl.pdf", width = 9, height = 5)

## my new fav way to plot site by Dx
ggplot(data = FAresults, aes(y=UNC_FA, x=site, fill = Dx)) +
  geom_boxplot(outlier.shape = NA) + scale_fill_grey(start = 0.4, end = 1) +
  geom_point(position=position_jitterdodge(jitter.width=0.2), size=3, aes(color=site)) 







# ggplot(thick_raw_pooled, aes(x=Dx, y=Thickness)) +
#   geom_boxplot(outlier.shape = NA) + 
#   geom_jitter(position = position_jitter(width = .1), size=3, aes(color=site)) +
#   facet_wrap(~ROI) + labs(y="Thickness (raw)")
### run on resid values
ttestres_pFA_resid = ttestpvals(FApooled.melted,'ROI','FA_resid ~ Dx')
ttestres_pFA_resid$site = 'pooled'
### do separately in all site and recombine
ttestres_pFA_resid_CMH = ttestpvals(subset(FApooled.melted,site=='CMH'),'ROI','FA_resid ~ Dx')
ttestres_pFA_resid_CMH$site = 'CMH'
ttestres_pFA_resid_ZHH = ttestpvals(subset(FApooled.melted,site=='ZHH'),'ROI','FA_resid ~ Dx')
ttestres_pFA_resid_ZHH$site = 'ZHH'
ttestres_pFA_resid_MRC = ttestpvals(subset(FApooled.melted,site=='MRC'),'ROI','FA_resid ~ Dx')
ttestres_pFA_resid_MRC$site = 'MRC'

alleffects_FA_resid = rbind(ttestres_pFA_resid,ttestres_pFA_resid_CMH, ttestres_pFA_resid_ZHH, ttestres_pFA_resid_MRC)

alleffects_FA_resid$Tract <- NA
for (i in 1:nrow(alleffects_FA_resid)) {
  alleffects_FA_resid$Tract[i] <- as.character(strsplit(row.names(alleffects_FA_resid)[i],'_')[[1]][1])
}


tmp <- subset(alleffects_FA_resid,site=='pooled')
alleffects_FA_resid_order <- as.character(tmp[order(as.numeric(tmp$cohenD)),'Tract'])
alleffects_FA_resid$Tract <- factor(alleffects_FA_resid$Tract,levels = alleffects_FA_resid_order)

ggplot(subset(alleffects_FA_resid,site=='pooled'), aes(y=as.numeric(cohenD), x=Tract)) +
  geom_errorbar(aes(ymin = (cohenD - cohenD.se), ymax =  (cohenD + cohenD.se)), width = 0.2) +
  geom_point(size = 8, shape = 18) + 
  geom_point(data=subset(alleffects_FA,site!='pooled'), aes(y=cohenD, x=Tract, colour = site), size = 3) +
  geom_hline(yintercept=0) +
  scale_y_reverse() + 
  labs(x='',y="Group Effect Size (Hedges g)") +
  coord_flip()

ggsave("alleffects_FA_resid_effsiz.pdf",width = 8, height = 8)

alleffects_FA_resid_horder <- gsub('_FA','',PooledRois_horder,fixed=T)
#alleffects_FA_resid$Tract <- factor(alleffects_FA_resid$Tract,levels = alleffects_FA_resid_horder)

