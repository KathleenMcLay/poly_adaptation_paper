install.packages("fANCOVA")
library(fANCOVA)
library(ggplot2)
library(tools)

dxy.files <- list.files(path = "/Users/kathleenmclay/Google Drive/Papers/Avneet_polygenic_adaptation/3_results/data_tables/dxy/", pattern=".dxy", all.files=FALSE)
inversions <- read.csv("/Users/kathleenmclay/Google Drive/Papers/Avneet_polygenic_adaptation/inversions.csv")
k = 1

pdf("/Users/kathleenmclay/Desktop/dxy_plot.pdf", height=7, width=7)
for (j in 1:length(dxy.files)){
  # read in data 
  dxy_dat <- read.table(paste("/Users/kathleenmclay/Google Drive/Papers/Avneet_polygenic_adaptation/3_results/data_tables/dxy/", dxy.files[j], sep = ""), header=TRUE)
  dxy_dat$pop1 <- as.character(dxy_dat$pop1)
  dxy_dat$pop2 <- as.character(dxy_dat$pop2)
  # remove windows with NA values for pi
  dxy_dat <- drop_na(dxy_dat) 
  # add window midpoint value to data for plotting
  for (i in 1:nrow(dxy_dat)) {
    dxy_dat$mid[i]=floor((dxy_dat$window_pos_1[i]+dxy_dat$window_pos_2[i])/2)
  }
  
  #subset to population, calculate loess smoothing values
  dxy_dat <- dxy_dat[(dxy_dat$pop1 == "0" & dxy_dat$pop2 == "2") | (dxy_dat$pop1 == "2" & dxy_dat$pop2 == "0"), ]
  dxy_dat$index <- 1:nrow(dxy_dat)
  dxy.lo <- loess.as(dxy_dat$index, dxy_dat$avg_dxy, degree = 0, criterion ="aicc", user.span = NULL, plot = F)
  dxy.lo.pred <- predict(dxy.lo)
  dxy_dat <- cbind(dxy_dat,dxy.lo.pred)

  #get current inversion start and end, to add to plot
  current_inversion <- inversions[k,]
  inv_start <- current_inversion[1,3]
  inv_end <- current_inversion[1,4]
  k = k + 1 
    
  # plot 
  title <- tools::file_path_sans_ext(dxy.files[j])
  dxy.plot <- ggplot(dxy_dat, aes(x=mid/1000000, y=dxy.lo.pred)) + theme_bw() +
    geom_line() +
    annotate("rect", fill = "lightblue", alpha = 0.5, 
             xmin = inv_start/1000000, xmax = inv_end/1000000,
             ymin = 0.025, ymax = 0.075) +
    xlab("Positions (Mbp)") +
    ylab("Dxy") +
    labs(color='0 - 2 Genotype\nComparison') +
    ggtitle(title) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches"),
          axis.title.x = element_text(margin = margin(t = 20), face="bold", size = 12), 
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), face="bold", size = 12),
          plot.title = element_text(vjust = 4, size = 10), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
          axis.text = element_text(size = 12))
  print(dxy.plot)
  
  
}
dev.off()