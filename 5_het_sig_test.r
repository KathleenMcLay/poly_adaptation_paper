### Mann-Whitney U test for heterozygosity~genotype - comparing presumed heterozygous and homozygous groups 

mds_info <- read.table("/.../lpca_tables/8_mds_info.txt")
het_data <- read.table("/.../lpca_tables/9_heterozygosity.txt")
het_data$genotype <- as.character(het_data$genotype)


# 0-1 comparison
het01 <- subset(het_data, genotype %in% c('0','1'))
new <- tibble(mds_coord=character(), name=character(), 
              PC1=numeric(), PC2=numeric(), genotype=character(), het=character())

man_u01 <- tibble(inversion=character(), group=character(), p_val=numeric())

for (i in 1:nrow(mds_info)) {
  for (j in 1:nrow(het01)) {
    if (mds_info$mds_coord[i] == het01$mds_coord[j]) {
      new <- rbind(new, het01[j,]) 
    }
  }
  test <- wilcox.test(new$het~new$genotype)
  p_val <- test$p.value
  inversion <- paste0(mds_info$chromosome[i], ":", mds_info$start[i], "-", mds_info$end[i])
  group <- "0-1"
  p_val <- cbind(p_val, inversion, group)
  
  man_u01 <- rbind(man_u01, p_val)
  
}

# 1-2 comparison
het12 <- subset(het_data, genotype %in% c('1','2'))
new <- tibble(mds_coord=character(), name=character(), 
              PC1=numeric(), PC2=numeric(), genotype=character(), het=character())
man_u12 <- tibble(inversion=character(), group=character(), p_val=numeric())

for (i in 1:nrow(mds_info)) {
  for (j in 1:nrow(het12)) {
    if (mds_info$mds_coord[i] == het12$mds_coord[j]) {
      new <- rbind(new, het12[j,]) 
    }
  }
  test <- wilcox.test(new$het~new$genotype)
  p_val <- test$p.value
  inversion <- paste0(mds_info$chromosome[i], ":", mds_info$start[i], "-", mds_info$end[i])
  group <- "1-2"
  p_val <- cbind(p_val, inversion, group)
  man_u12 <- rbind(man_u12, p_val)
  
}

het_sig <- rbind(man_u01, man_u12)
write.table(het_sig, "file_path/het_sig.txt")