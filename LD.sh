module load vcftools

FILE="/QRISdata/Q6684/for_paper/1_data/pa_vo_sf5_final.vcf.gz"
DIR="/QRISdata/Q6684/for_paper/1_data/data_tables/pvp_popgen/LD"

# calculate LD and convert to table for all samples
vcftools --gzvcf ${FILE} --chr ${1} --geno-r2 --maf 0.05 --max-missing 0 --out ${DIR}/${1}
cat ${DIR}/${1}.geno.ld | perl /QRISdata/Q6656/reformat/emerald2windowldcounts.pl > ${DIR}/${1}_ld.txt
rm ${DIR}/${1}.geno.ld

# calculate LD and convert to table for D01
D01_list="/home/uqkmcla4/scripts/avneet/D01.txt"
vcftools --gzvcf ${FILE} --keep ${D01_list} --chr ${1} --geno-r2 --maf 0.05 --max-missing 0 --out ${DIR}/D01_${1}

rm ${DIR}/D01_${1}.geno.ld

# calculate LD and convert to table for H01
H01_list="/home/uqkmcla4/scripts/avneet/H01.txt"
vcftools --gzvcf ${FILE} --keep ${H01_list} --chr ${1} --geno-r2 --maf 0.05 --max-missing 0 --out ${DIR}/H01_${1}

rm ${DIR}/H01_${1}.geno.ld