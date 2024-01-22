### Combine all samples per population using GATK CombineGVCF, joint call using GATK GenotypeGVCF by chromosome and concatenate output files using bcftools concat 
### GNU parallel by population

module load gatk
module load bcftools

# create directory variables 
VCDIR="/QRISdata/Q6656/avneet/variant"
CBDIR="/QRISdata/Q6656/avneet/joint/combined"
SPDIR="/QRISdata/Q6656/avneet/joint/split"
JCDIR="/QRISdata/Q6656/avneet/joint"

# create a list of all sample files for the current population
ls -d -1 "${VCDIR}/${1}"*.g.vcf.gz > "${1}_gvcf.list"

# combine all samples into a single vcf 
gatk CombineGVCFs -R /QRISdata/Q6656/avneet/SLv141Asm_Ch20RN.fasta --variant ${1}_gvcf.list -O ${CBDIR}/${1}_cohort.g.vcf.gz 

# create an array of scaffold names to use as interval values with GenotypeGVCFs
declare -a arr=("SLv141Ch1" "SLv141Ch2" "SLv141Ch3" "SLv141Ch4" "SLv141Ch5" "SLv141Ch6" "SLv141Ch7" "SLv141Ch8" "SLv141Ch9" "SLv141Ch10" "SLv141Ch11" "SLv141Ch12" "SLv141Ch13" "SLv141Ch14" "SLv141Ch15" "SLv141Ch16" "SLv141Ch17" "SLv141Ch18" "SLv141Ch19" "SLv141Ch20")

# joint call with GenotypeGVCFs
for i in "${arr[@]}";
do
    # create a variable to isolate the scaffold no. for file naming  
    if [ ${#i} -eq 9 ]; then
        modified_i="${i:0:8}0${i:8:1}"
    else
        modified_i="${i}"
    fi
    chr="${modified_i: -2}"
    # run GenotypeGVCFs per scaffold using -intervals 
    gatk --java-options "-Xmx8g" GenotypeGVCFs \
       -R /QRISdata/Q6656/avneet/SLv141Asm_Ch20RN.fasta \
       -V ${CBDIR}/${1}_cohort.g.vcf.gz \
       -O ${SPDIR}/${1}_${chr}_cohort_jntcl.vcf.gz \
       --include-non-variant-sites \
       --intervals $i &
done

wait

# concatenate the chromosomes for each sample, remove the chromsome files once done  

#create list variables 
file_list="${JCDIR}/${1}_file_list.txt"
sort_list="${JCDIR}/${1}_sorted_file_list.txt"

# find all files for the current population in the output directory 
find "$SPDIR" -type f -name "${1}_*.vcf.gz" >> "$file_list"

# print the file to the file_list variable, sort it in the sort list variable and then concatenate and index the file 
if [ -s "$file_list" ]; then
    awk '{print $1}' "$file_list" | sort -t'_' -k5 -n > "$sort_list" 
    # Concatenate and index 
    bcftools concat --threads 12 --file-list "$sort_list" --output "${JCDIR}/${1}_jntcl.vcf.gz" 
    bcftools index --threads 12 -t "${JCDIR}/${1}_jntcl.vcf.gz" 
    
else
    echo "No match for: $1" 
fi  
