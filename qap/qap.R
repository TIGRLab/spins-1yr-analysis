library(ggplot2)
theme_set(theme_minimal(base_size = 22))
anat = read.csv('data/qap_anatomical_spatial.csv')
anat$subject = as.character(anat$subject)
anat = anat[-grep("_P", anat$subject),]
anat$site = as.factor(unlist(lapply(anat$subject, function (x) { unlist(strsplit(x,c("_"),fixed=T))[2]} )))

ggplot(anat, aes(x=site, fill = site, color=site, y=snr)) + geom_violin() + ggtitle("QAP T1 SNR")
ggsave(file = "qap_t1_snr.png")

ggplot(anat, aes(x=site, fill = site, color=site, y=cnr)) + geom_violin() + ggtitle("QAP T1 CNR")
ggsave(file = "qap_t1_cnr.png")

ggplot(anat, aes(x=site, fill = site, color=site, y=qi1)) + geom_violin() + ggtitle("QAP T1 Qi1 Artifact Detection (lower is better)")
ggsave(file = "qap_t1_qi1.png")

ggplot(anat, aes(x=site, fill = site, color=site, y=efc)) + geom_violin() + ggtitle("QAP T1 Head Motion (lower is better)") + 
  ylab("Entropy Focus Criterion (EFC)")
ggsave(file = "qap_t1_efc.png")


### Functional spatial
func = read.csv('data/qap_functional_spatial_temporal.csv')
func$subject = as.character(func$subject)
func = func[-grep("_P", func$subject),]
func$site = as.factor(unlist(lapply(func$subject, function (x) { unlist(strsplit(x,c("_"),fixed=T))[2]} )))

ggplot(func, aes(x=site, fill = site, color=site, y=snr)) + geom_violin() + ggtitle("QAP RST SNR")
ggsave(file = "qap_rst_snr.png")

ggplot(func, aes(x=site, fill = site, color=site, y=ghost_y)) + geom_violin() + ggtitle("QAP RST Ghost to Signal Ration (lower is better)")
ggsave(file = "qap_rst_ghost.png")

ggplot(func, aes(x=site, fill = site, color=site, y=efc)) + geom_violin() + ggtitle("QAP RST Head Motion") + 
  ylab("Entropy Focus Criterion (EFC)")
ggsave(file = "qap_rst_efc.png")

### Functional temporal
temp = read.csv('data/qap_functional_temporal.csv')
temp$subject = as.character(temp$subject)
temp = temp[-grep("_P", temp$subject),]
temp$site = as.factor(unlist(lapply(temp$subject, function (x) { unlist(strsplit(x,c("_"),fixed=T))[2]} )))

temp = subset(temp, subject != "SPN01_CMH_0031")
ggplot(temp, aes(x=site, fill = site, color=site, y=dvars)) + geom_violin() + ggtitle("QAP RST Standardized DVARS (low = better)")
ggsave(file = "qap_rst_dvars.png")

ggplot(temp, aes(x=site, fill = site, color=site, y=quality)) + geom_violin() + ggtitle("QAP RST Quality (lower is better)")
ggsave(file = "qap_rst_quality,png")

ggplot(temp, aes(x=site, fill = site, color=site, y=mean_fd)) + geom_violin() + ggtitle("QAP RST Head Motion") + 
  ylab("Mean Fractional Displacement - Jenkinson")
ggsave(file = "qap_rst_jenkinson.png")
