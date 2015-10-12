library(ggplot2)
theme_set(theme_classic(base_size = 22))
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