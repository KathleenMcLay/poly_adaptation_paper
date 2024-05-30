### Calculate HET for each individual for all SNPs in the putative inversion

module load bcftools 
module load vcftools

file="/QRISdata/Q6684/for_paper/1_data/pa_vo_sf5_final.vcf.gz"
output="/QRISdata/Q6684/working_data/het/${1}"

# subset the vcf file to the inversion region
bcftools filter ${file} --threads 12 --regions ${1} -o ${output}

# calcuate heterozygosity per individual 
vcftools \
    --vcf ${output} \
    --het \
    --out /QRISdata/Q6684/working_data/het/${1}