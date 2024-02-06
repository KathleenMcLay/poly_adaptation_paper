### Calculate Pi and Dxy for each putative inversion
### GNU parallel by chromosome/inverison

#activate the miniconda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bio

module load bcftools 

DIR="/g/data/ht96/McLay_UQ/avneet_paper/pi_dxy"
FILE="/g/data/ht96/McLay_UQ/avneet_paper/pa_as_sf5_final.rm_reheader.vcf.gz"

# filter to the chromosome
bcftools filter ${FILE} --regions ${1} --output ${DIR}/${2}.vcf 

# Zip and index the inversion files 
bgzip -@ 12 ${DIR}/${2}.vcf
tabix ${DIR}/${2}.vcf.gz

#create a file to allocate each sample to based on its genotype format: D01203 G0
POPS="${DIR}/${2}_POPS.txt"

#read in inversion genotype data 
genotypes="/home/564/km6006/Scripts/avneet_paper/inversion_genotypes.txt"

# for each inversion check the genotype, and assign the genotype value as the 'population' for that individual
while IFS=$'\t' read -r inversion name genotype population; do
    echo "current $inversion has "$genotype" genotype"
    if [ "$inversion" == "$2" ] && [ "$genotype" -eq 0 ]; then
        echo -e "$name\t$genotype" >> "$POPS"
    elif [ "$inversion" == "$2" ] && [ "$genotype" -eq 1 ]; then
        echo -e "$name\t$genotype" >> "$POPS"
    elif [ "$inversion" == "$2" ] && [ "$genotype" -eq 2 ]; then
        echo -e "$name\t$genotype" >> "$POPS"        
    else
        echo "no match"
    fi
done < "$genotypes"

# Run pixy - calculate pi and dxy values
pixy --stats pi dxy --vcf ${DIR}/${2}.vcf.gz --populations ${POPS} --window_size 10000 --output_prefix ${2} --output_folder ${DIR} --n_cores 12