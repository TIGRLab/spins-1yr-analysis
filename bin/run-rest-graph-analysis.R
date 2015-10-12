#!/usr/bin/env R
library(ggplot2)
library(reshape2)

# load data, reformat
data <- read.csv('rest-graph-metrics.csv')
metrics<-c('conn.mean', 'ddist.var', 'conn.var', 'ddist.plaw', 'eff', 'ddist.mean', 'smworld', 'robtarget', 'robrandom')
data <- melt(data, metrics=metrics)

# ttests: across all sites, and within sites.
data_cmh <- subset(data, site == 'CMH')
data_mrc <- subset(data, site == 'MRC')
data_zhh <- subset(data, site == 'ZHH')

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

# per metric plot
ggplot(subset(data, variable == 'conn.mean'), aes(x=site, y=value, fill=dx)) + geom_boxplot()
ggplot(subset(data, variable == 'conn.var'), aes(x=site, y=value, fill=dx)) + geom_boxplot()
ggplot(subset(data, variable == 'ddist.var'), aes(x=site, y=value, fill=dx)) + geom_boxplot()
ggplot(subset(data, variable == 'ddist.plaw'), aes(x=site, y=value, fill=dx)) + geom_boxplot()
ggplot(subset(data, variable == 'eff'), aes(x=site, y=value, fill=dx)) + geom_boxplot()
ggplot(subset(data, variable == 'ddist.mean'), aes(x=site, y=value, fill=dx)) + geom_boxplot()
ggplot(subset(data, variable == 'smworld'), aes(x=site, y=value, fill=dx)) + geom_boxplot()
ggplot(subset(data, variable == 'robtarget'), aes(x=site, y=value, fill=dx)) + geom_boxplot()
ggplot(subset(data, variable == 'robrandom'), aes(x=site, y=value, fill=dx)) + geom_boxplot()

