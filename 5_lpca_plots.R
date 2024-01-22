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

###PCA PLOTS###

#read in files 
out <- read.table("/.../results/lpca_tables/9_heterozygosity.txt")
out$genotype <- as.character(out$genotype)

#create empty data from to store filtered data 
new <- tibble(mds_coord=character(), name=character(), 
              PC1=numeric(), PC2=numeric(), genotype=character(), population=character())

#matching cluster data from 'out' to the outlier cluster list and producing the PCA plots 
pdf("/.../results/PCA_plots.pdf", height=7, width=7)
for (i in 1:nrow(mds_info)) {
  PC1_perc <- mds_info$PC1_perc[i]
  PC2_perc <- mds_info$PC2_perc[i]
  for (j in 1:nrow(out)) {
    if (mds_info$mds_coord[i] == out$mds_coord[j]) {
      new <- rbind(new, out[j,]) 
    }
  }
  pca_plot <- new %>%
    ggplot(., aes(x=PC1, y=PC2, col=genotype)) + geom_point(size = 2) + theme_bw() + 
    scale_color_manual(name="Genotype", values=c("red","purple","blue")) +
    ggtitle(paste(mds_info$mds_coord[i], " ", mds_info$chromosome[i])) +
    xlab(paste("PC",1," (",round(PC1_perc,2),"% PVE)",sep="")) +
    ylab(paste("PC",2," (",round(PC2_perc,2),"% PVE)",sep="")) +
    #geom_label_repel(aes(label = name), label.size = NA, fill = "NA", show.legend = FALSE, segment.color = "transparent") +
    theme(legend.title = element_text(size=18, face="bold"), 
          legend.text=element_text(size=15),
          plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches"), 
          axis.title.x = element_text(margin = margin(t = 20), face="bold", size = 12), 
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), face="bold", size = 12),
          plot.title = element_text(vjust = 4, size = 10), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
          axis.text = element_text(size = 12),
          legend.position = "none")
  print(pca_plot)
  new <- tibble(mds_coord=character(), name=character(), 
                PC1=numeric(), PC2=numeric(), genotype=character(), population=character())
}
dev.off()


###HET PLOTS###

heterozygosity <- read.table("/.../results/lpca_tables/9_heterozygosity.txt")
heterozygosity$genotype <- as.character(heterozygosity$genotype)

new <- tibble(mds_coord=character(), name=character(), 
              PC1=numeric(), PC2=numeric(), genotype=character(), het=character())

pdf("/.../results/HET_plots.pdf", height=7, width=7)
for (i in 1:nrow(mds_info)) {
  for (j in 1:nrow(heterozygosity)) {
    if (mds_info$mds_coord[i] == heterozygosity$mds_coord[j]) {
      new <- rbind(new, heterozygosity[j,]) 
    }
  }
  het_plot <- new %>%
    ggplot(., aes(x=as.character(genotype), y=het, fill=as.character(genotype))) + geom_boxplot() + theme_bw() +
    ggtitle(paste(mds_info$mds_coord[i], " ", mds_info$chromosome[i])) +
    scale_fill_manual(name="Genotype", values=c("red","purple","blue")) +
    xlab("Genotype") + ylab("Heterozygosity") +
    theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches"),
          axis.title.x = element_text(margin = margin(t = 20), face="bold", size = 12), 
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), face="bold", size = 12),
          plot.title = element_text(vjust = 4, size = 10), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
          axis.text = element_text(size = 12))
  print(het_plot)
  new <- tibble(mds_coord=character(), name=character(), 
                  PC1=numeric(), PC2=numeric(), genotype=character(), het=character())
}
dev.off()