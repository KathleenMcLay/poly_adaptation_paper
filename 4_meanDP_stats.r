library(data.table)

# histogram of read depths - tail is very long, so subset to < 100 to see where the decline starts 
system.time(data <- fread("/QRISdata/Q6656/avneet/snp_fil/pa_vo_sf4.ldepth.txt"))
subdata <- subset(data, MEAN_DEPTH < 100, select = c("CHROM", "POS", "MEAN_DEPTH","VAR_DEPTH"))
jpeg(file="/QRISdata/Q6656/avneet/snp_fil/pa_vo_sf4.ldepth.jpeg")
mean_depth <- hist(subdata$MEAN_DEPTH, col="purple", main="Distribution of Mean Read Depth", xlab="Read Depth", breaks=800) 
dev.off()

# calculate 90th quantile, d*sqrd and mode
statdata <- subset(data, select = "MEAN_DEPTH")

v <- round(statdata$MEAN_DEPTH)
d <- mean(v)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

sink(file = "/QRISdata/Q6656/avneet/snp_fil/pa_vo_sf4.ldepth_stats.txt")
quantile(v, probs = 0.9)
print("d squared")
d+3*sqrt(d)
print("mode")
getmode(v)
sink(file = NULL)

dev.off()