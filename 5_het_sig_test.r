## Mann-W U test 
mds_info <- read.table("/.../8_mds_info_new.txt", header=TRUE)
het_data <- read.table("/.../all_het_results.txt")
het_data$genotype <- as.character(het_data$cluster)


#0-1
het01 <- subset(het_data, genotype %in% c('0','1'))
new <- tibble(mds_coord=character(), name=character(), 
              PC1=numeric(), PC2=numeric(), genotype=character(), het=character())

man_u01 <- tibble(inversion=character(), group=character(), p_val=numeric())

for (i in 1:nrow(mds_info)) {
  for (j in 1:nrow(het01)) {
    if (mds_info$inversion[i] == het01$inv[j]) {
      new <- rbind(new, het01[j,]) 
    }
  }
  test <- wilcox.test(new$het~new$genotype)
  p_val <- test$p.value
  inversion <- mds_info$inversion[i]
  group <- "0-1"
  p_val <- cbind(p_val, inversion, group)
  
  man_u01 <- rbind(man_u01, p_val)
  new <- tibble(mds_coord=character(), name=character(), 
                PC1=numeric(), PC2=numeric(), genotype=character(), het=character())
}

het12 <- subset(het_data, genotype %in% c('1','2'))
new <- tibble(mds_coord=character(), name=character(), 
              PC1=numeric(), PC2=numeric(), genotype=character(), het=character())
man_u12 <- tibble(inversion=character(), group=character(), p_val=numeric())

for (i in 1:nrow(mds_info)) {
  for (j in 1:nrow(het12)) {
    if (mds_info$inversion[i] == het12$inv[j]) {
      new <- rbind(new, het12[j,]) 
    }
  }
  
  unique_values <- unique(new$genotype)
  
  if (length(unique_values) == 2) {
    test <- wilcox.test(new$het~new$genotype)
    p_val <- test$p.value
    inversion <- mds_info$inversion[i]
    group <- "1-2"
    p_val <- cbind(p_val, inversion, group)
    man_u12 <- rbind(man_u12, p_val)
  } else {
    print(new)
  }
  new <- tibble(mds_coord=character(), name=character(), 
                PC1=numeric(), PC2=numeric(), genotype=character(), het=character())
}

het_sig <- rbind(man_u01, man_u12)
write.table(het_sig, "/.../het_significance.txt", sep = "\t")