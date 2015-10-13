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

# save plot
plt <- ggplot(df, aes(x=site, y=correlation, fill=dx))
plt + geom_boxplot() + scale_color_brewer(palette="Paired") + scale_fill_brewer(palette="Paired") + theme_minimal() + facet_wrap(~vid)
ggsave(file='ea-behavioural-analysis.jpg', scale=1)

