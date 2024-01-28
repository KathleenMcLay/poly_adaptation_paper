### Combine data tables, calculate percentage ranks, plot HI-C data 
### Code taken from Todesco et al. 2020 and Huang et al. 2020 with modifications

library(data.table)
library(tidyverse)
library(grid)
library(gridExtra)
library(tools)
library(stringr)
library(ggplot2)
library(dplyr)

### combine the interaction tables for the two samples, for each chromosome with a region of interest

H01_files <- list.files(path = "/g/data/ht96/McLay_UQ/avneet_paper/hic/100kb/tables", pattern="H1", all.files=FALSE)
D01_files <- list.files(path = "/g/data/ht96/McLay_UQ/avneet_paper/hic/100kb/tables", pattern="D1", all.files=FALSE)

# read in the interaction tables for each chromosome 
for (i in 1:length(H01_files)) {
  print(H01_files[i])
  data_pop1 <- fread(paste("/g/data/ht96/McLay_UQ/avneet_paper/hic/100kb/tables/", H01_files[i], sep = ""))
  colnames(data_pop1) <- c("chr2","bin2","chr1","bin1","inv_links", "genotype") 
  for (j in 1:length(D01_files)){
      if (str_split(H01_files[i], "_", n = 2, simplify=TRUE)[2] == str_split(D01_files[j], "_", n = 2, simplify=TRUE)[2]) {
          data_pop2 <- fread(paste("/g/data/ht96/McLay_UQ/avneet_paper/hic/100kb/tables/", D01_files[j], sep = ""))
          colnames(data_pop2) <- c("chr2","bin2","chr1","bin1","noninv_links", "genotype")
          print(D01_files[j])
      }
  }   

  # combine the data from the two tables, and calc. link comparisons (difference in links between pop1 and pop2)
  data_pop1[,noninv_links:=data_pop2[,noninv_links]]
  data_pop1[,link_comparison:=inv_links - noninv_links] 
  data_pair <- data_pop1
  # write the new table to file
  chrom <- str_split(H01_files[i], "_", n = 3, simplify=TRUE)[2]
  table_file <- paste("/g/data/ht96/McLay_UQ/avneet_paper/hic/100kb/paired_tables/", chrom, ".txt", sep = "")
  write.table(data_pair, table_file)
}

Read in the file containing the 'inversion' details, create inversion id column 
inversions <- read.csv("/g/data/ht96/McLay_UQ/avneet_paper/hic/inversions.csv")
inversions$id <- paste0(inversions$chrom,".",inversions$MDS,".",inversions$outliers)

### Editing the data for the plot
hic_pairs <- list.files(path = "/g/data/ht96/McLay_UQ/avneet_paper/hic/100kb/paired_tables/",pattern=".txt", all.files=FALSE)

pdf("/g/data/ht96/McLay_UQ/avneet_paper/hic/HiC_plots.pdf", height=5, width=7)
for (i in 1:nrow(inversions)) {
  chosen_inversion <- inversions$id[i]
  print(paste0("preparing data: ", chosen_inversion))
  chosen_chr <- inversions$chrom[i]
  chosen_start <- inversions$start[i]
  chosen_end <- inversions$end[i]
  chosen_middle = ((floor(chosen_start/1000000) + floor(chosen_end/1000000))/2)*1000000 
  chosen_height = ((floor(chosen_end/100000) - floor(chosen_start/100000))/2) 
  print(paste0(chosen_inversion, " height is: ", chosen_height))
  view_start <- inversions$chr_start[i]
  view_end <- inversions$chr_end[i]
  
  for (j in 1:length(hic_pairs)){
      if (tools::file_path_sans_ext(hic_pairs[j]) == chosen_chr) {
          print(paste0("reading file for: ", hic_pairs[j]))
          pop_pair <- fread(paste("/g/data/ht96/McLay_UQ/avneet_paper/hic/100kb/paired_tables/", hic_pairs[j], sep = ""))
          colnames(pop_pair) <- c("n","chr2","bin2","chr1","bin1","inv_links", "genotype", "noninv_links","link_comparison")
          #pop_pair <- subset(pop_pair[,c(2:9)])
      }
  }
  print(head(pop_pair, n=5))
  window_size <- 100000
  #point_size = 100/((view_end - view_start)/10000)
  
  # prepare data for plot. 
  plegend <- as_tibble(pop_pair[chr2 ==chosen_chr & chr1 ==chosen_chr & bin2 >= view_start & bin2 <= view_end &
                                  bin1 >= view_start & bin1 <= view_end]) %>%  
    filter(bin1 >= bin2) %>%
    mutate(distance = bin1 - bin2) %>%
    group_by(distance) %>%
    mutate(percent_ranks = percent_rank(link_comparison)) %>%
    ungroup() %>%
    mutate(link_comparison = case_when(bin1 == bin2 ~ 99, #all bins that are the same are marked as 99
                                       link_comparison > 0.3 ~ 0.3, #all links > 0.3 are converted to 0.3
                                       link_comparison < -0.3 ~ -0.3, #all links < -0.3 are converted to 0.3
                                       link_comparison == 0 ~ 0.00001, #all links that are 0, are changed to 0.00001, presumably to allow for plotting
                                       TRUE ~ link_comparison)) %>% 
    mutate(link_comparison = na_if(link_comparison, link_comparison == 99))  %>% #bins that are the same are excluded (shown as grey triangles on plot)
    #filter(bin1 >= bin2) %>%
    mutate(x = (bin1 + bin2)/2,y=(bin1-bin2)/200000) %>% #CREATE X and Y variables 
    group_by(y,x,link_comparison,percent_ranks, distance) %>%
    expand(count = seq(1:4)) %>% #THIS STUFF IS TO USE GEOM_POLYGON! - it has four sides! 
    mutate(bin_id = paste0(x,"_",y)) %>%
    ungroup() %>%
    mutate(y = case_when(count == 2 ~ y+0.5,
                         count == 4 ~ y - 0.5,
                         TRUE ~ y),
           x = case_when(count == 1 ~ x -(window_size/2),
                         count == 3 ~ x + (window_size/2),
                         TRUE ~ x)) %>%
    mutate(y = case_when(y < 0 ~ 0,
                         TRUE ~ y)) 

  # Make the plot
  print(paste0("making plot for: ", chosen_inversion)) 
  hic_plot <- (ggplot(plegend,aes()) + geom_polygon(aes(x=x/1000000,y=y,fill=link_comparison,group = bin_id),size=0) +
                  scale_fill_gradient2(low = "#377EB8", mid = "white",
                                      high = "#E41A1C", midpoint = 0,name="HiC link\ndifference",limits=c(-0.3,0.3)) +
                  annotate("segment",x=floor(chosen_start/1000000),y=0,
                          yend=chosen_height,xend=chosen_middle/1000000,
                          alpha=1,color="black",size=0.5) +
                  annotate("segment",x=floor(chosen_end/1000000),y=0,
                          yend=chosen_height,xend=chosen_middle/1000000,
                          alpha=1,color="black",size=0.5) +
                  theme_linedraw() + ylab("") + xlab("Positions (Mbp)") +
                  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_blank(),
                        panel.border = element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank(),
                        legend.position="none") ) + 
                  ggtitle(chosen_inversion)
  print(hic_plot)
}
dev.off()


### ALL CHROMS 

# Read in the file containing the 'inversion' details, create inversion id column 
chromosomes <- read.csv("/g/data/ht96/McLay_UQ/avneet_paper/hic/chromosomes.csv")

### Editing the data for the plot
hic_pairs <- list.files(path = "/g/data/ht96/McLay_UQ/avneet_paper/hic/100kb/paired_tables/",pattern=".txt", all.files=FALSE)

pdf("/g/data/ht96/McLay_UQ/avneet_paper/hic/HiC_plots.pdf", height=5, width=7)
for (i in 1:nrow(chromosomes)) {
  chosen_chr <- chromosomes$chrom[i]
  view_start <- chromosomes$chr_start[i]
  view_end <- chromosomes$chr_end[i]
  
  for (j in 1:length(hic_pairs)){
      if (tools::file_path_sans_ext(hic_pairs[j]) == chosen_chr) {
          print(paste0("reading file for: ", hic_pairs[j]))
          pop_pair <- fread(paste("/g/data/ht96/McLay_UQ/avneet_paper/hic/100kb/paired_tables/", hic_pairs[j], sep = ""))
          colnames(pop_pair) <- c("n","chr2","bin2","chr1","bin1","inv_links", "genotype", "noninv_links","link_comparison")
          #pop_pair <- subset(pop_pair[,c(2:9)])
      }
  }
  print(head(pop_pair, n=5))
  window_size <- 100000
  #point_size = 100/((view_end - view_start)/10000)
  
  # prepare data for plot. 
  plegend <- as_tibble(pop_pair[chr2 ==chosen_chr & chr1 ==chosen_chr & bin2 >= view_start & bin2 <= view_end &
                                  bin1 >= view_start & bin1 <= view_end]) %>%  
    filter(bin1 >= bin2) %>%
    mutate(distance = bin1 - bin2) %>%
    group_by(distance) %>%
    mutate(percent_ranks = percent_rank(link_comparison)) %>%
    ungroup() %>%
    mutate(link_comparison = case_when(bin1 == bin2 ~ 99, #all bins that are the same are marked as 99
                                       link_comparison > 0.3 ~ 0.3, #all links > 0.3 are converted to 0.3
                                       link_comparison < -0.3 ~ -0.3, #all links < -0.3 are converted to 0.3
                                       link_comparison == 0 ~ 0.00001, #all links that are 0, are changed to 0.00001, presumably to allow for plotting
                                       TRUE ~ link_comparison)) %>% 
    mutate(link_comparison = na_if(link_comparison, link_comparison == 99))  %>% #bins that are the same are excluded (shown as grey triangles on plot)
    #filter(bin1 >= bin2) %>%
    mutate(x = (bin1 + bin2)/2,y=(bin1-bin2)/200000) %>% #CREATE X and Y variables 
    group_by(y,x,link_comparison,percent_ranks, distance) %>%
    expand(count = seq(1:4)) %>% #THIS STUFF IS TO USE GEOM_POLYGON! - it has four sides! 
    mutate(bin_id = paste0(x,"_",y)) %>%
    ungroup() %>%
    mutate(y = case_when(count == 2 ~ y+0.5,
                         count == 4 ~ y - 0.5,
                         TRUE ~ y),
           x = case_when(count == 1 ~ x -(window_size/2),
                         count == 3 ~ x + (window_size/2),
                         TRUE ~ x)) %>%
    mutate(y = case_when(y < 0 ~ 0,
                         TRUE ~ y)) 

  # Make the plot
  print(paste0("making plot for: ", chosen_chr)) 
  hic_plot <- (ggplot(plegend,aes()) + geom_polygon(aes(x=x/1000000,y=y,fill=link_comparison,group = bin_id),size=0) +
                  scale_fill_gradient2(low = "#377EB8", mid = "white",
                                      high = "#E41A1C", midpoint = 0,name="HiC link\ndifference",limits=c(-0.3,0.3)) +
                  theme_linedraw() + ylab("") + xlab("Positions (Mbp)") +
                  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_blank(),
                        panel.border = element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank(),
                        legend.position="none") ) + 
                  ggtitle(chosen_chr)
  print(hic_plot)
}
dev.off()