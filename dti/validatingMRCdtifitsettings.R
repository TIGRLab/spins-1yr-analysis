## load libraries
library("reshape2", lib.loc="/opt/quarantine/R-extras/3.1.1/build")
library("ggplot2", lib.loc="/opt/quarantine/R-extras/3.1.1/build")

## load post settings change data from erin's projects
FA_postfix <- read.csv("/projects/edickie/analysis/dtifit_testing/SPINS/enigma-dti/enigmaDTI-FA-results.csv")
## load the pre settings change data from /archive/data-2.0
FA_prefix <- read.csv("/archive/data-2.0/SPINS/data/enigmaDTI/enigmaDTI-FA-results.csv")

## grab ROI names from the column names
ROIs <- names(FA_postfix)[2:ncol(FA_postfix)]

#change the column names to set up for merging
names(FA_postfix) <- c('id',paste0(ROIs,'post'))
names(FA_prefix) <- c('id',paste0(ROIs,'pre'))

# merge the two files
merged_FA <- merge(FA_prefix,FA_postfix)

# set up an empty dataframe to collect the differences
FA_diff <- FA_prefix
names(FA_diff) <- c('id',paste0(ROIs,'diff'))
FA_diff$id <- as.character(FA_diff$id)
FA_diff[ ,1:ncol(FA_diff) ] <- NA

## calculate the differences for each column
for (i in 1:length(ROIs)) {
  roi = ROIs[i]
  FA_diff$id[i] <- as.character(merged_FA$id[i])
  FA_diff[ ,i+1] <- merged_FA[ ,paste0(roi,'post')] - merged_FA[ ,paste0(roi,'pre')]
}

## make a site column for the plotting
FA_diff$site <- NA
FA_diff$site[grep('CMH',as.character(FA_diff$id))]<-"CMH"
FA_diff$site[grep('ZHH',as.character(FA_diff$id))]<-"ZHH"
FA_diff$site[grep('MRC',as.character(FA_diff$id))]<-"MRC"
FA_diff$site <- as.factor(FA_diff$site)

## print summary distributions for each site
summary(subset(FA_diff, site == "CMH"))
summary(subset(FA_diff, site == "ZHH"))
summary(subset(FA_diff, site == "MRC"))

## melt the data for ggplot
change.melted <- melt(FA_diff,
                          id.vars=c('id','site'),
                          measure.vars=paste0(ROIs,'diff'),
                          variable.name="ROIdiff",
                          value.name="Change")

## make prettier tick labels for ggplot
change.melted$ROI <- gsub('_FAdiff','',as.character(change.melted$ROIdiff))

## make the pretty plot
myplot <- ggplot(change.melted, aes(x=ROI, y=Change, colour=site)) +
          geom_boxplot(outlier.shape = NA) + geom_point(aes(colour=site)) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          labs(y="Change in FA", title="Change calculated FA after changing BET parameters for MRC site")

## save the plot
ggsave("~/code/spins-1yr-analysis/dti/ValidatingMRCFAthres.pdf",myplot,width=10,height=5)
