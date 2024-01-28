### Plot Pi with loess smoothing 

install.packages("fANCOVA")
library(fANCOVA)
library(ggplot2)
library(tools)

pi.files <- list.files(path = "/Users/kathleenmclay/Google Drive/Papers/Avneet_polygenic_adaptation/3_results/data_tables/pi/", pattern=".pi", all.files=FALSE)
inversions <- read.csv("/Users/kathleenmclay/Google Drive/Papers/Avneet_polygenic_adaptation/inversions.csv")
k = 1

pdf("/Users/kathleenmclay/Desktop/pi_plot.pdf", height=7, width=7)
for (j in 1:length(pi.files)){
  # read in data 
  pi_dat <- read.table(paste("/Users/kathleenmclay/Google Drive/Papers/Avneet_polygenic_adaptation/3_results/data_tables/pi/", pi.files[j], sep = ""), header=TRUE)
  pi_dat$pop <- as.character(pi_dat$pop)
  # remove windows with NA values for pi
  pi_dat <- drop_na(pi_dat) 
  # add window midpoint value to data for plotting
  for (i in 1:nrow(pi_dat)) {
    pi_dat$mid[i]=floor((pi_dat$window_pos_1[i]+pi_dat$window_pos_2[i])/2)
  }
  
  #subset to population, calculate loess smoothing values
  pi_dat_0 <- subset(pi_dat, pop == "0")
  pi_dat_0$index <- 1:nrow(pi_dat_0)
  pi.lo <- loess.as(pi_dat_0$index, pi_dat_0$avg_pi, degree = 1, criterion ="aicc", user.span = NULL, plot = F)
  pi.lo.pred <- predict(pi.lo)
  pi_dat_0 <- cbind(pi_dat_0,pi.lo.pred)
    
  pi_dat_1 <- subset(pi_dat, pop == "1")
  pi_dat_1$index <- 1:nrow(pi_dat_1)
  pi.lo <- loess.as(pi_dat_1$index, pi_dat_1$avg_pi, degree = 1, criterion ="aicc", user.span = NULL, plot = F)
  pi.lo.pred <- predict(pi.lo)
  pi_dat_1 <- cbind(pi_dat_1, pi.lo.pred)
  
  pi_dat_2 <- subset(pi_dat, pop == "2")
  pi_dat_2$index <- 1:nrow(pi_dat_2)
  pi.lo <- loess.as(pi_dat_2$index, pi_dat_2$avg_pi, degree = 1, criterion ="aicc", user.span = NULL, plot = F)
  pi.lo.pred <- predict(pi.lo)
  pi_dat_2 <- cbind(pi_dat_2,pi.lo.pred)
  
  # combine datasets 
  rm(pi_dat)
  pi_dat <- rbind(pi_dat_0, pi_dat_1, pi_dat_2)
  pi_dat$pop <- as.character(pi_dat$pop)
  
  #get current inversion start and end, to add to plot
  current_inversion <- inversions[k,]
  inv_start <- current_inversion[1,3]
  inv_end <- current_inversion[1,4]
  k = k + 1 
  
  # plot
  title <- tools::file_path_sans_ext(pi.files[j])
  pi.plot <- ggplot(pi_dat, aes(x=mid/1000000, y=pi.lo.pred)) + theme_bw() +
    geom_line(aes(colour=pop)) +
    scale_fill_manual(name="Genotype", values=c("red","purple","blue")) +
    annotate("rect", fill = "lightblue", alpha = 0.5, 
             xmin = inv_start/1000000, xmax = inv_end/1000000,
             ymin = 0, ymax = 0.1) +
    xlab("Positions (Mbp)") +
    ylab(expression(paste(pi))) +
    ggtitle(title) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches"),
          axis.title.x = element_text(margin = margin(t = 20), face="bold", size = 12), 
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), face="bold", size = 12),
          plot.title = element_text(vjust = 4, size = 10), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
          axis.text = element_text(size = 12),
          legend.position = "none")
  print(pi.plot)
}
dev.off()