#!/usr/bin/env R
ggplot(plot, aes(x=sites, y=means)) + geom_point() + facet_wrap(~rois)
library(ggplot2)

# create a dataset with the diagnois
x <- read.csv('CorticalMeasuresENIGMA_ThickAvg.csv')
subjects <- read.table("list.csv", header=TRUE, quote="\"")
names(subjects) <-c("SubjID", "dx")
subjects$SubjID <- as.character(subjects$SubjID) # convert to string for later replacement

count <- 1
for (subj in subjects$SubjID) {
    subjects$SubjID[count] <- paste(subj, '_01', sep='')
    count <- count + 1
}
subjects$SubjID <- factor(subjects$SubjID)
data <- merge(x, subjects, by="SubjID")

# create a factor for site number
idx_cmh <- grepl('CMH', data$SubjID)
idx_mrc <- grepl('MRC', data$SubjID)
idx_zhh <- grepl('ZHH', data$SubjID)

data$site[idx_cmh] <- 1
data$site[idx_mrc] <- 2
data$site[idx_zhh] <- 3
data$site <- factor(data$site)

# get the data of inerest into the right format
cols <- colnames(data)
idxrois <- grepl('R_|L_', cols)
rois <- cols[idxrois]
keep <- c(rois, 'SubjID', 'dx', 'site')
data <- subset(data, select=keep) # remove non-raw data from dataframe
data <- melt(data, roi=rois) # reformat
data <- data[data$dx != 'x', ] # remove the to-be-ignored participants

# per ROI plot
ggplot(data, aes(x=site, y=value, fill=dx)) + geom_boxplot() + facet_wrap(~variable)

# overall site plot
ggplot(data, aes(x=site, y=value, fill=dx) + geom_boxplot()

# ttests: across all sites, and within sites.
data_cmh <- subset(data, site == 1)
data_mrc <- subset(data, site == 2)
data_zhh <- subset(data, site == 3)

ttest_all = do.call("rbind", lapply(split(data, data$variable), function(x) t.test(value ~ dx, x)$p.value))
ttest_cmh = do.call("rbind", lapply(split(data_cmh, data_cmh$variable), function(x) t.test(value ~ dx, x)$p.value))
ttest_mrc = do.call("rbind", lapply(split(data_mrc, data_mrc$variable), function(x) t.test(value ~ dx, x)$p.value))
ttest_zhh = do.call("rbind", lapply(split(data_zhh, data_zhh$variable), function(x) t.test(value ~ dx, x)$p.value))

ttest_all = cbind(ttest_all, p.adjust(ttest_all, method="fdr"), p.adjust(ttest_all, method="bonferroni"))
ttest_cmh = cbind(ttest_cmh, p.adjust(ttest_cmh, method="fdr"), p.adjust(ttest_cmh, method="bonferroni"))
ttest_mrc = cbind(ttest_mrc, p.adjust(ttest_mrc, method="fdr"), p.adjust(ttest_mrc, method="bonferroni"))
ttest_zhh = cbind(ttest_zhh, p.adjust(ttest_zhh, method="fdr"), p.adjust(ttest_zhh, method="bonferroni"))

colnames(ttest_all) = c("raw", "fdr", "bonferroni")
colnames(ttest_cmh) = c("raw", "fdr", "bonferroni")
colnames(ttest_mrc) = c("raw", "fdr", "bonferroni")
colnames(ttest_zhh) = c("raw", "fdr", "bonferroni")

ttest_all = as.data.frame(ttest_all)
ttest_cmh = as.data.frame(ttest_cmh)
ttest_mrc = as.data.frame(ttest_mrc)
ttest_zhh = as.data.frame(ttest_zhh)

sig_all = subset(ttest_all, fdr < 0.05)
sig_cmh = subset(ttest_cmh, fdr < 0.05)
sig_mrc = subset(ttest_mrc, fdr < 0.05)
sig_zhh = subset(ttest_zhh, fdr < 0.05)


