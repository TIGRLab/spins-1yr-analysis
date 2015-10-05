# An analysis of the enigma whitematter measures
library(reshape2)

fa = read.csv('data/enigmaDTI-FA-results.csv')
fa$id = as.character(fa$id)
# extract site from scan id
fa$site = as.factor(unlist(lapply(fa$id, function (x) { unlist(strsplit(x,c("_"),fixed=T))[2]} )))

demog = read.csv('data/demog.csv')
demog$id = as.character(demog$id)
demog = subset(demog, dx != "x")
demog$dx = factor(demog$dx)
data = merge(fa,demog)

df = melt(data)
# t-tests http://stats.stackexchange.com/questions/128894/p-value-correction-for-multiple-t-tests
p.raw = do.call("rbind", 
                lapply(split(df, df$variable), 
                       function(x) t.test(value ~ dx, x)$p.value))
p.values = cbind(p.raw, 
                 p.adjust(p.raw, method="fdr"), 
                 p.adjust(p.raw, method="bonferroni"))
colnames(p.values) = c("raw", "fdr", "bonferroni")
