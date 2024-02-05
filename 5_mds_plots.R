### Code adapted from Todesco et al. 2020 and Huang et al. 2020

library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggrepel)

###MDS PLOT###

#read in files 
mds_info <- read.table("/.../results/lpca_tables/8_mds_info.txt")
mds_info$mds <- tolower(mds_info$mds)
outliers <- read.table("/.../results/lpca_tables/3_outlier_windows.txt")
win_regions <- read.table("/.../results/lpca_tables/1_win_regions.txt")

#Filter the outliers and mds data(win_regions) by mds_info which lists the outlier clusters along the genome for all MDS. Then plot.  
pdf("/.../results/MDS_plots.pdf", height=7, width=7)
for (i in 1:nrow(mds_info)) {
  outliers %>%
    filter(mds_coord == mds_info$mds_coord[i]) -> outs
  win_regions %>% select(1:4, mds_info$mds[i], 45) %>%
      filter(chrom == mds_info$chromosome[i]) %>%
      mutate(outlier=case_when(n %in% outs$n ~ "Outlier", 
                              TRUE ~ "Non-outlier")) -> new
  colnames(new) <- c("chrom", "start", "end", "mid", "mds", "n", "outlier")
  genome_plot <- new %>%
    ggplot(., aes(x=mid/1000000, y=mds, color=outlier)) + geom_point() + theme_bw() +
    scale_color_manual(values=c("gray80","midnightblue")) +
    ggtitle(paste(mds_info$mds_coord[i], " ", mds_info$chromosome[i])) +
    xlab("Position (Mbp)") + ylab(toupper(mds_info$mds[i])) +
    theme(legend.position="none", 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches"),
          axis.title.x = element_text(margin = margin(t = 20), face="bold", size = 12), 
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), face="bold", size = 12),
          plot.title = element_text(vjust = 4, size = 10), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
          axis.ticks = element_line(colour="black"), panel.background=element_rect(colour="black"), panel.border=element_rect(colour="black", size=0.75))
  print(genome_plot)
}
dev.off()
