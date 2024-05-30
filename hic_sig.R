library(data.table)
library(tidyverse)
library(grid)
library(gridExtra)
library(tools)

hic_files <- list.files(path = "/.../data_tables/hic",pattern=".txt", all.files=FALSE)

# Read inversions 
inversions <- read.csv("/.../data_tables/hic/inversions.csv")
inversions$id <- paste0(inversions$chrom,".",inversions$MDS,".",inversions$outliers)

perc_peaks <- tibble(
  chr2 = character(),
  bin2 = integer(),
  chr1 = character(),
  bin1 = integer(),
  inv_links = numeric(),
  genotype = character(),
  noninv_links = numeric(),
  link_comparison = character(),
  distance = numeric(),
  percent_ranks = numeric(),
  perc = numeric()
)

for (i in 1:nrow(inversions)) {

  for (j in 1:length(hic_files)) {
    file <- tools::file_path_sans_ext(hic_files[j])
    if (file == inversions$chrom[i]) {
      pop_pair <- fread(paste("/.../data_tables/hic/", hic_files[j], sep = ""))
      pop_pair <- subset(pop_pair[,c(2:9)])
    }
  }
    # calculate percentage rank by distance 
    pop_pair %>%  
      filter(bin1 >= bin2) %>%
      mutate(distance = bin1 - bin2) %>%
      group_by(distance) %>%
      mutate(percent_ranks = percent_rank(link_comparison)) %>%
      ungroup() -> pop_pair

    # extract the inversion break point interaction percentage rank 
    peak <- subset(pop_pair, bin2 == inversions$start._window[i] & bin1==inversions$end_window[i]) 
    peak$perc <- floor(peak$percent_ranks*100)
    perc_peaks <- rbind(perc_peaks, peak)

}
write.csv(perc_peaks, "/.../data_tables/hic/hic_int_comparison_sig.csv")
