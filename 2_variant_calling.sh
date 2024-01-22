### Call variants with GATK HaplotypeCaller by chromosome and concatenate output files using bcftools concat 
### GNU parallel by sample 

module load gatk 
module load bcftools

# create directory variables 
CLDIR="/g/data/ht96/McLay_UQ/avneet_paper/clean"
SPDIR="/g/data/ht96/McLay_UQ/avneet_paper/variant/split"
VCDIR="/g/data/ht96/McLay_UQ/avneet_paper/variant"

# create an array of scaffold names to use as interval values with HaplotypeCaller
declare -a arr=("SLv141Ch1" "SLv141Ch2" "SLv141Ch3" "SLv141Ch4" "SLv141Ch5" "SLv141Ch6" "SLv141Ch7" "SLv141Ch8" "SLv141Ch9" "SLv141Ch10" "SLv141Ch11" "SLv141Ch12" "SLv141Ch13" "SLv141Ch14" "SLv141Ch15" "SLv141Ch16" "SLv141Ch17" "SLv141Ch18" "SLv141Ch19" "SLv141Ch20")

# run HaplotypeCaller 
for i in "${arr[@]}";
do
    # create a variable to isolate the scaffold no. for file naming  
    if [ ${#i} -eq 9 ]; then
        modified_i="${i:0:8}0${i:8:1}"
    else
        modified_i="${i}"
    fi
    chr="${modified_i: -2}"
    # run HaplotypeCaller per chromosome using -intervals
    gatk --java-options "-Xmx8g" HaplotypeCaller \
        --input ${CLDIR}/${1}_PCRm_cln_srt.bam \
        --output ${SPDIR}/${1}_${chr}_vrnt.g.vcf.gz \
        --reference /scratch/ht96/km6006/SLv141Asm_Ch20RN.fasta \
        --emit-ref-confidence GVCF \
        --intervals $i &
done

wait

# concatenate the scaffolds for each sample

#create list variables 
file_list="${VCDIR}/${1}_file_list.txt"
sort_list="${VCDIR}/${1}_sorted_file_list.txt"

# find all files for the current sample in the output directory 
find "$SPDIR" -type f -name "${1}_*.vcf.gz" >> "$file_list"

# print the file to the file_list variable, sort it in the sort list variable and then concatenate 
if [ -s "$file_list" ]; then
    awk '{print $1}' "$file_list" | sort -t'_' -k5 -n > "$sort_list" 
    # Concatenate and index 
    bcftools concat --threads 6 --file-list "$sort_list" --output "${VCDIR}/${1}_vrnt.g.vcf.gz" 
    bcftools index --threads 6 -t "${VCDIR}/${1}_vrnt.g.vcf.gz" 

else
    echo "No match for: $1" 
fi  

