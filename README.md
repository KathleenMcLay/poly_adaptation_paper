## Polygenic architecture of adaptation paper (James et al.)

Complete code for analysis contributed to Kaur et al. polygenic adaptation paper

### Data 

WGS data - H01 and D01 natural populations 
HiC data - H01 and D01 single individuals

### Alignment + SNP filtering 

- data_prep.sh 
    
    Trim, align, clean, mark PCR duplicates for WGS data. Generate QC, alignment and coverage reports 

- variant_calling.sh
    
    Call variants with GATK HaplotypeCaller by chromosome and concatenate output files using bcftools concat 

- joint_calling.sh
    
    Combine all samples per population using GATK CombineGVCF, joint genotype using GATK GenotypeGVCF by chromosome and concatenate output files using bcftools concat 

- SNP_filter.sh
    
    SNP filtering with GATK SelectVariants, VariantFiltration and Vcftools

- meanDP_stats.r 
    
    Visualise, calculate 90th quantile, d*sqrd and mode for mean read depth to determine vcftools --max-meanDP filter 

- SNP_fil_PCA.R
    
    Create a PCA of the final SNP filtered dataset 

### Local PCA + pop gen analysis 

- local_PCA.R

    Local PCA (MDS) Code adapted from Todesco et al. 2020 and Huang et al. 2020 with minor modifications

- lpca_plots.R

    Plot MDS, PCA and Heterozygosity results. Code adapted from Todesco et al. 2020 and Huang et al. 2020 with minor modifications

- inv_pca.R

    PCA and with kmeans clustering for all SNPs in each putative inversion and plot 

- heterozygosity.sh

    Calculate HET for each individual for all SNPs in the putative inversion with VCFtools 

- het_plot.R

    Plot heterozygosity results for each putative inversion 

- het_sig_test.r

    Test significant difference in het between heterozygous and homozygous inversion clusters identified in PCA of local PCA results. Test: Mann-Whitney U test for heterozygosity~genotype

- Fst.sh

    Calculate fst between homozygous genotype clusters for each putative inversion identified by local PCA

- Fst_plot.sh

    Plot Fst results 

- pixy.sh

    Calculate fst between homozygous genotype clusters for each putative inversion identified by local PCA

- pi_plot.sh

    Plot pi with loess smoothing 

 - dxy_plot.sh

    Plot dxy with loess smoothing    

### Pop gen analysis 

- _pvp_pop_gen.sh 

    Calculte pi, dxy and fst between populations for every chromosome.

- LD.sh

    Calculate LD within H01 and D01 populations, and for both populations collectively. 


### HiC analysis 

- hic_analysis.sh

    Trim, align, create Tag directories and interaction matricies with HiC data 

- hic_plot.r

    Plot difference in HiC interactions between H01 and D01 samples. Code adapted from Todesco et al. 2020 with minor modifications

- hic_sig.r 

    Calculate the percentage rank by distance for the difference in interactions for D01 and H01 between the inversion breakpoints

