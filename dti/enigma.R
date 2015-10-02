# An analysis of the enigma whitematter measures
fa = read.csv('data/enigmaDTI-FA-results.csv')
fa$id = as.character(fa$id)
fa$site = as.factor(unlist(lapply(fa$id, function (x) { unlist(strsplit(x,c("_"),fixed=T))[2]} )))