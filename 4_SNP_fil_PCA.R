### Create PCAs for SNP filtered dataset 

library(tidyverse)
library(SNPRelate)
library(Matrix)
library(ggrepel)
library(ggplot2)

samples <- read_tsv("/home/uqkmcla4/scripts/avneet/PCA_samples.tsv")
vcf_fn <- ("/QRISdata/Q6656/avneet/snp_fil/pa_as_sf5_10kbprune.vcf")

# Convert VCF file 
snpgdsVCF2GDS(vcf_fn, "senecio.gds", method="copy.num.of.ref") 
# Assign the new file to an object 
genofile <- snpgdsOpen("senecio.gds")
# Run PCA 
senecio_pca <- snpgdsPCA(genofile, autosome.only=FALSE)

# Extract PCs and make dataset with samples 
PC1 <- senecio_pca$eigenvect[,1]
PC2 <- senecio_pca$eigenvect[,2]
PCA <- cbind(PC1, PC2)
PCA <- cbind(samples, PCA)
write.table(PCA, "/QRISdata/Q6656/avneet/snp_fil/pa_vo_sf5_final.recode_PCA.txt", sep="\t")

# Plot the PCA 
tiff("/QRISdata/Q6656/avneet/snp_fil/pa_vo_sf5_final.recode_PCA.tiff", units="in", width=9, height=7, res=300)
ggplot(PCA, aes(x=PC1, y=PC2, col=population)) + geom_point(size = 2) + theme_bw() +
  scale_colour_discrete(name="Populations") +
  geom_label_repel(aes(label = sample), label.size = NA, fill = "NA", show.legend = FALSE, segment.color = "transparent") +
  xlab("PC1") + ylab("PC2") +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches"), axis.title = element_text(size = 12),axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)), text = element_text(family = "Helvetica") ) 
dev.off()
file.remove("/home/uqkmcla4/scripts/senecio.gds")
