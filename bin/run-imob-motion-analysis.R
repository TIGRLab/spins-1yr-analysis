library(ggplot2)
library(stats)

directory <- '/archive/data-2.0/SPINS/data/imob/'
files <- Sys.glob(paste(directory, '/*/*motion.1D', sep=''))

# initialize data frame
df <- data.frame(fdmean=NA, site=NA, id=NA)

for (file in files) {

    # initialize MINIFRAME :D
    x <- data.frame(fdmean=NA, site=NA, id=NA)

    # get framewise displacement measures (mean, sum)
    data <- read.csv(file, sep="", header=FALSE) # any # of whitespace as delim
    data <- data.matrix(data) # convert to matrix
    data <- rowSums(abs(diff(data))) # sum of absoloute value of derivative

    # load in variables of interest
    x$fdmean <- mean(data)
    x$id <- paste(strsplit(basename(file), '_')[[1]][1:3], collapse="_")
    x$site <- strsplit(x$id, '_')[[1]][2]

    # append to master data frame
    df <- rbind(df, x)
}

# add diagnosis
dx <- read.csv('/projects/jdv/data/spins/1yr/outputs/bin/list.csv', sep=' ')
dx <- subset(dx, dx!='x')
dx$dx <- factor(dx$dx)
df <- merge(df, dx)

# site differences?
df_hc <- subset(df, dx=='h')
df_sz <- subset(df, dx=='p')
summary(aov(fdmean ~ site, df_hc)) # neither are significant
summary(aov(fdmean ~ site, df_sz)) #

# plot party
plt <- ggplot(df, aes(x=site, y=fdmean, fill=dx))
plt + geom_boxplot() + scale_color_brewer(palette="Paired") + scale_fill_brewer(palette="Paired") + theme_minimal()
ggsave(file='imob-headmotion.jpg', scale=1)

