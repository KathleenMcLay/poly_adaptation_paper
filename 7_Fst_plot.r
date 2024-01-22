### Plot Fst with loess smoothing 

install.packages("fANCOVA")
library(fANCOVA)
library(ggplot2)
library(tools)

fst.files <- list.files(path = "/Users/kathleenmclay/Google Drive/Papers/Avneet_polygenic_adaptation/results/data_tables/fst", pattern=".fst", all.files=FALSE)
inversions <- read.csv("/Users/kathleenmclay/Google Drive/Papers/Avneet_polygenic_adaptation/methods/inversions.csv")
k = 1

pdf("/Users/kathleenmclay/Desktop/Fst_weighted.pdf", height=7, width=7) 
for (j in 1:length(fst.files)){
  # read in data, remove any na values, add index column set negative Fst values to 0
  fst.dat <- read.table(paste("/Users/kathleenmclay/Google Drive/Papers/Avneet_polygenic_adaptation/results/data_tables/fst/", fst.files[j], sep = ""), header=TRUE)
  fst.dat <- na.omit(fst.dat)
  fst.dat$index <- 1:nrow(fst.dat)
  fst.dat$WEIGHTED_FST[fst.dat$WEIGHTED_FST < 0] <- 0
  fst.dat$MEAN_FST[fst.dat$MEAN_FST < 0] <- 0
  
  # add window midpoint value to data for plotting
  for (i in 1:nrow(fst.dat)) {
    fst.dat$mid[i]=floor((fst.dat$BIN_START[i]+fst.dat$BIN_END[i])/2)
  }
  
  # determine optimal loess smoothing values
  fst.lo <- loess.as(fst.dat$index, fst.dat$WEIGHTED_FST, degree = 2, criterion ="aicc", user.span = NULL, plot = F)
  fst.lo.pred <- predict(fst.lo)
  
  #add loess smoothing values to the dataset 
  fst.dat <- cbind(fst.dat,fst.lo.pred)
  
  #get current inversion start and end, to add to plot
  current_inversion <- inversions[k,]
  inv_start <- current_inversion[1,3]
  inv_end <- current_inversion[1,4]
  k = k + 1 

  #plot
  title <- tools::file_path_sans_ext(fst.files[j])
  fst.plot <- ggplot(fst.dat, aes(x=mid/1000000, y=WEIGHTED_FST)) + theme_bw() +
    geom_point(col="grey") +
    geom_line(data=fst.dat, aes(x=mid/1000000, y=fst.lo.pred),col="navyblue") +
    annotate("rect", fill = "lightblue", alpha = 0.5, 
             xmin = inv_start/1000000, xmax = inv_end/1000000,
             ymin = 0, ymax = 1) +
    xlab("Positions (Mbp)") +
    ylab("Fst (weighted)") +
    ggtitle(title) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches"),
          axis.title.x = element_text(margin = margin(t = 20), face="bold", size = 12), 
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), face="bold", size = 12),
          plot.title = element_text(vjust = 4, size = 10), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
          axis.text = element_text(size = 12),
          legend.position = "none")
  print(fst.plot)
}
dev.off()