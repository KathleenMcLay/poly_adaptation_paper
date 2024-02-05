### Plot heterozygosity results for each putative inversion 

pca_files <- list.files(path = "/.../pca", pattern="genotypes.txt", all.files=FALSE)
het_files <- list.files(path = "/.../het", pattern=".het", all.files=FALSE)

het_data <- data.frame()

pdf("/.../HET_plots_whole_region.pdf", height=7, width=7)
for (i in 1:length(pca_files)) {
  pca_file <- read.table(paste("/.../pca/", pca_files[i], sep = ""))

  for (j in 1:length(het_files)){
    if (str_split(pca_files[i], "-", n = 2, simplify=TRUE)[1] == str_split(het_files[j], "-", n = 2, simplify=TRUE)[1]) {
    het_file <- read.table(paste("/.../het/", het_files[j], sep = ""), header=TRUE)
    colnames(het_file) <- c("name", "observed_hom", "expected_hom", "sites", "F")
    het_file$het_sites <- het_file$sites - het_file$observed_hom
    het_file$het <- het_file$het_sites/het_file$sites
    het_file$cluster <- pca_file$cluster
    table.name <- paste("/.../het/", tools::file_path_sans_ext(het_files[j]), "_het.txt", sep = "")
    write.table(het_file, table.name, sep="\t")
    het_data <- rbind(het_data, het_file[c("name", "observed_hom", "expected_hom", "sites", "F", "het_sites", "het", "cluster", "inv")])
    write.table(het_data, "/.../all_het_results.txt", sep="\t")
    }
  }
  title <- tools::file_path_sans_ext(pca_files[i])  
  het_plot <- ggplot(het_file, aes(x=as.character(cluster), y=het, fill=as.character(cluster))) + geom_boxplot() + theme_bw() +
    ggtitle(title) +
    scale_fill_manual(name="Genotype", values=c("red","purple","blue")) +
    xlab("Genotype") + ylab("Heterozygosity") +
    theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches"),
          axis.title.x = element_text(margin = margin(t = 20), face="bold", size = 12), 
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), face="bold", size = 12),
          plot.title = element_text(vjust = 4, size = 10), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
          axis.text = element_text(size = 12))
  print(het_plot)
}
dev.off()