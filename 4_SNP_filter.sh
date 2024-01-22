### SNP filtering
### for dataset including non-variant sites exclude non-variant site filter

module load gatk
module load vcftools 
module load bcftools 

JDIR="/QRISdata/Q6656/avneet/joint"
DIR="/QRISdata/Q6656/avneet/snp_fil"
ref="/QRISdata/Q6656/avneet/SLv141Asm_Ch20RN.fasta"

ds="pa_vo"

### merge joint called sample VCF files 
file_list="${DIR}/filelist.txt"
find "$JDIR" -type f -name "*jntcl.vcf.gz" >> "$file_list"

bcftools merge --threads 24 --file-list "$file_list" --output ${DIR}/${ds}_merged.vcf.gz 
bcftools index --threads 24 -t ${DIR}/${ds}_merged.vcf.gz 

### rename samples in VCF
bcftools reheader -s /home/uqkmcla4/scripts/avneet/rename_samples.tsv -o ${DIR}/${ds}_merged_final.vcf.gz ${DIR}/${ds}_merged.vcf.gz
bcftools index --threads 24 -t ${DIR}/${ds}_merged_final.vcf.gz 

### Remove invariant sites and indels 
gatk SelectVariants \
    -R ${ref} \
    -V ${DIR}/${ds}_merged_final.vcf.gz \
    --exclude-non-variants \
    --select-type-to-include SNP \
    -O ${DIR}/${ds}_sf0.vcf.gz

### gatk filtering - adds PASS to the filter field, otherwise, if failed adds the name/s of the failed filter - does not remove sites 
gatk VariantFiltration \
    -V ${DIR}/${ds}_sf0.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O ${DIR}/${ds}_sf1.vcf.gz \
    > ${DIR}/VF_prog.out 2> ${DIR}/VF_prog.err

### Filter for max missing per site. The --remove-filtered-all option removes all sites that have been filtered out by the previous step
vcftools \
    --gzvcf ${DIR}/${ds}_sf1.vcf.gz \
    --max-missing 0.8 \
    --remove-filtered-all \
    --recode --recode-INFO-all --stdout | gzip -c > ${DIR}/${ds}_sf2.vcf.gz

### Remove samples with too much missing data 

# generate list of missing data per sample
vcftools \
    --gzvcf ${DIR}/${ds}_sf2.vcf.gz \
    --missing-indv \
    --out ${DIR}/${ds}_sf2

# create a list of samples with more than x missing data 
awk '$5 > 0.5' ${DIR}/${ds}_sf2.imiss | cut -f1 > ${DIR}/${ds}_sf2_lowDP.indv

# filter out the samples 
vcftools \
     --gzvcf ${DIR}/${ds}_sf2.vcf.gz \
     --remove ${DIR}/${ds}_sf2_lowDP.indv \
     --recode --recode-INFO-all --stdout | gzip -c > ${DIR}/${ds}_sf3.vcf.gz

### remove sites with too much missing data for any population 

# check missing site data - per population
populations="/home/uqkmcla4/scripts/avneet/populations.txt"

while IFS= read -r pop; do
    vcftools \
        --gzvcf ${DIR}/${ds}_sf3.vcf.gz \
        --keep /home/uqkmcla4/scripts/avneet/${pop}.txt \
        --missing-site \
        --out ${DIR}/${ds}_sf3_${pop}
done < "$populations"

wait

# concatenate the population lists 
touch "${DIR}/${ds}_sf3_popmiss.txt"

for file in "${DIR}"/*${ds}*.lmiss; do
        cat "$file" >> "${DIR}/${ds}_sf3_popmiss.txt"
done

# filter the list for those sites where the missing data for any population is > 50%  
awk '!/CHR/ && $6 > 0.5 {print $1, $2}' ${DIR}/${ds}_sf3_popmiss.txt > ${DIR}/${ds}_sf3_badloci.txt

# filter out sites where the missing data for any population is > 50%  
vcftools \
    --gzvcf ${DIR}/${ds}_sf3.vcf.gz \
    --exclude-positions ${DIR}/${ds}_sf3_badloci.txt \
    --recode --recode-INFO-all --stdout | gzip -c > ${DIR}/${ds}_sf4.vcf.gz
 
### Filter for mean depth -- here dataset is split into high coverage and low coverage samples, which require different mean depth filter thresholds

# visualise the distribution of mean depth per site
vcftools \
    --gzvcf ${DIR}/${ds}_sf4.vcf.gz \
    --site-mean-depth \
    --out ${DIR}/${ds}_sf4

# Filter for max-meanDP and min-meanDP, biallelic sites only 
vcftools \
    --gzvcf ${DIR}/${ds}_sf4.vcf.gz \
    --max-meanDP 20 \
    --min-meanDP 2 \
    --max-alleles 2 \
    --recode --recode-INFO-all --out ${DIR}/${ds}_sf5_final

# calculate missing data per sample for the final dataset 
vcftools \
    --vcf ${DIR}/${ds}_sf5_final.recode.vcf \
    --missing-indv \
    --out ${DIR}/${ds}_sf5_final.recode

# create bcf file 
bcftools convert -O b ${DIR}/${ds}_sf5_final.recode.vcf > ${DIR}/${ds}_sf5_final.recode.bcf
bcftools index -f ${DIR}/${ds}_sf5_final.recode.bcf