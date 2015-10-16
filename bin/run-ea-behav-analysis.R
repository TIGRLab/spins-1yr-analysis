#!/usr/bin/env R

# get a list of the available files
directory <- '/projects/jdv/data/spins/1yr/data/ea'
files <- Sys.glob(paste(directory, '/*corr-push*', sep=''))

# initalize data frame
df <- data.frame(correlation=NA, n.pushes.per.minute=NA, id=NA, site=NA, vid=NA)

# load in diagnosis from list.csv
dx <- read.csv('/projects/jdv/data/spins/1yr/outputs/bin/list.csv', sep=' ')
dx <- subset(dx, dx!='x')
dx$dx <- factor(dx$dx)

for (file in files) {

    # get subject name
    subj <- paste(strsplit(basename(file), '_')[[1]][1:3], collapse="_")

    # load in data
    x <- read.csv(paste(subj, '_01_corr-push.csv', sep=''))
    x$id <- subj
    x$site <- strsplit(subj, '_')[[1]][2]

    # add in video numbers
    x$vid <- c(1,2,3,4,5,6,7,8,9)

    # merge in this subject to the master dataframe
    df <- rbind(df, x)
}

# add in diagnosis
df <- merge(df, dx)

# save plots
plt <- ggplot(df, aes(x=site, y=correlation, fill=dx))
plt + geom_boxplot() + scale_color_brewer(palette="Paired") + scale_fill_brewer(palette="Paired") + theme_minimal() + facet_wrap(~vid)
ggsave(file='ea-behav-analysis_corrs.jpg', scale=1)

plt <- ggplot(df, aes(x=site, y=n.pushes.per.minute, fill=dx))
plt + geom_boxplot() + scale_color_brewer(palette="Paired") + scale_fill_brewer(palette="Paired") + theme_minimal() + facet_wrap(~vid)
ggsave(file='ea-behav-analysis_npush.jpg', scale=1)

df_cmh <- subset(df, site=='CMH')
df_mrc <- subset(df, site=='MRC')
df_zff <- subset(df, site=='ZHH')

ttest_all = do.call("rbind", lapply(split(df, df$vid), function(x) t.test(correlation ~ dx, x)$p.value))
ttest_cmh = do.call("rbind", lapply(split(df_cmh, df_cmh$vid), function(x) t.test(value ~ dx, x)$p.value))
ttest_mrc = do.call("rbind", lapply(split(df_mrc, df_mrc$vid), function(x) t.test(value ~ dx, x)$p.value))
ttest_zhh = do.call("rbind", lapply(split(df_zhh, df_zhh$vid), function(x) t.test(value ~ dx, x)$p.value))

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
