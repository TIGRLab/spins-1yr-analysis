# An analysis of the enigma whitematter measures
library(reshape2)
library(ggplot2)

fa = read.csv('data/enigmaDTI-FA-results.csv')
fa$id = as.character(fa$id)
# extract site from scan id
fa$site = as.factor(unlist(lapply(fa$id, function (x) { unlist(strsplit(x,c("_"),fixed=T))[2]} )))

demog = read.csv('data/demog.csv')
demog$id = as.character(demog$id)
demog = subset(demog, dx != "x")
demog$dx = factor(demog$dx)
data = merge(fa,demog)

df = melt(data, id.vars = c("id", "site", "dx"), variable.name="ROI", value.name="FA")
# ggplot(data = df, aes(x=site, y=FA, colour=dx)) + geom_boxplot() + facet_wrap(~ROI)

# t-tests http://stats.stackexchange.com/questions/128894/p-value-correction-for-multiple-t-tests
p.ttest = do.call("rbind", 
                lapply(split(df, df$ROI), 
                       function(x) t.test(FA ~ dx, x)$p.value))
p.ttest = cbind(p.ttest, 
                 p.adjust(p.ttest, method="fdr"), 
                 p.adjust(p.ttest, method="bonferroni"))
colnames(p.ttest) = c("raw", "fdr", "bonferroni")
p.ttest = as.data.frame(p.ttest)
sig.ttest = subset(p.ttest, fdr < 0.05)

# residuals of FA ~ site
df$FA_resid =  melt(sapply(split(df, df$ROI), function(x) { resid(lm(FA ~ site, x)) }))$value
#ggplot(data = df, aes(x=dx, y=FA_resid, colour=dx)) + geom_boxplot() + facet_wrap(~ROI)
p.resid = do.call("rbind", 
                  lapply(split(df, df$ROI), 
                         function(x) t.test(FA_resid ~ dx, x)$p.value))
p.resid = cbind(p.resid, 
                p.adjust(p.resid, method="fdr"), 
                p.adjust(p.resid, method="bonferroni"))
colnames(p.resid) = c("raw", "fdr", "bonferroni")
p.resid = as.data.frame(p.resid)
sig.resid = subset(p.resid, fdr < 0.05)

df.sig.resid = subset(df, ROI %in% rownames(sig.resid))
ggplot(data = df.sig.resid, aes(x=site, y=FA_resid, colour=dx)) + geom_boxplot() + facet_wrap(~ROI)

# independant 
p.value = c()
beta = c()
n = c()
site_list = c()
roi = c()
index = 1
for (site_x in c("CMH", "ZHH", "MRC")) {
  site_df = subset(df, site == site_x)
  for (ROI_x in unique(as.character(df$ROI))) {
    roi_df = subset(site_df, ROI == ROI_x)
    mod = lm(FA ~ dx, roi_df)
    p.value[index] = anova(mod)[1,5]
    beta[index] = mod$coef[2]
    n[index] = nrow(roi_df)  # missing values counted, get this from model? 
    site_list[index] = site_x
    roi[index] = ROI_x
    index = index + 1
  }
}
df.ind = data.frame(site = as.factor(site_list), 
                    roi = as.factor(roi), 
                    beta = beta, 
                    p.value = p.value, 
                    n = n)

df.ind$p0.5 = ifelse(df.ind$p.value < 0.05, T, F)
ggplot(df.ind, aes(x = roi, 
                   y = beta, 
                   size = -log(p.value, base =10), 
                   shape = p0.5, 
                   color = as.factor(site))) + 
  geom_point() + 
  geom_hline(aes(yintercept=0)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

