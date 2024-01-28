## polygenic architecture of adaptation paper (Kaur et al.)

complete code for analysis contributed to Kaur et al. polygenic adaptation paper

### Data 

WGS data - H01 and D01 natural populations 
HiC data - H01 and D01 single individuals

### Alignment + SNP filtering 

- 1_data_prep.sh 
    
    Trim, align, clean, mark PCR duplicates for WGS data. Generate QC, alignment and coverage reports 

- 2_variant_calling.sh
    
    Call variants with GATK HaplotypeCaller by chromosome and concatenate output files using bcftools concat 

- 3_joint_calling.sh
    
    Combine all samples per population using GATK CombineGVCF, joint genotype using GATK GenotypeGVCF by chromosome and concatenate output files using bcftools concat 

- 4_SNP_filter.sh
    
    SNP filtering with GATK SelectVariants, VariantFiltration and Vcftools

- 4_meanDP_stats.r 
    
    Visualise, calculate 90th quantile, d*sqrd and mode for mean read depth to determine vcftools --max-meanDP filter 

- 4_SNP_fil_PCA.R
    
    Create a PCA of the final SNP filtered dataset 

### Local PCA + pop gen analysis 

- 5_LPCA_het_pca.R

    Local PCA (MDS), PCA of outliers and Heterozygosity of outliers. Code adapted from Todesco et al. 2020 and Huang et al. 2020 with minor modifications

- 5_het_sig_test.r

    Test significant difference in het between heterozygous and homozygous inversion clusters identified in PCA of local PCA results. Test: Mann-Whitney U test for heterozygosity~genotype

- 5_lpca_plots.R

    Plot MDS, PCA and Heterozygosity results. Code adapted from Todesco et al. 2020 and Huang et al. 2020 with minor modifications

- 7_Fst.sh

    Calculate fst between homozygous genotype clusters for each putative inversion identified by local PCA

- 7_Fst_plot.sh

    Plot Fst results 

- 8_pixy.sh

    Calculate fst between homozygous genotype clusters for each putative inversion identified by local PCA

- 8_pi_plot.sh

    Plot pi with loess smoothing 

 - 8_dxy_plot.sh

    Plot dxy with loess smoothing    


### Pop gen analysis 

- 9_pi_dxy_pop.sh 

    Calculte pi, dxy and fst between populations for every chromosome.


### HiC analysis 

- 6_hic_analysis.sh

    Trim, align, create Tag directories and interaction matricies with HiC data 

- 6_hic_plot.r

    Plot difference in HiC interactions between H01 and D01 samples. Code adapted from Todesco et al. 2020 with minor modifications



