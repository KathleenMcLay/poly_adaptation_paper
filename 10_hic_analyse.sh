dir="/g/data/ht96/McLay_UQ/avneet_paper/hic"

# array of chroms 
declare -a arr=("SLv141Ch1" "SLv141Ch2" "SLv141Ch3" "SLv141Ch4" "SLv141Ch5" "SLv141Ch6" "SLv141Ch7" "SLv141Ch8" "SLv141Ch9" "SLv141Ch10" "SLv141Ch11" "SLv141Ch12" "SLv141Ch13" "SLv141Ch14" "SLv141Ch15" "SLv141Ch16" "SLv141Ch17" "SLv141Ch18" "SLv141Ch19" "SLv141Ch20")

# create Hi-C interaction matrix and convert it to a table 
for i in "${arr[@]}";
do
    /home/564/km6006/homer/bin/analyzeHiC ${dir}/1_tags/${1}_HiC-TAG -chr ${i} -res 100000 -coverageNorm > ${dir}/${1}_${i}_hic_matrix.txt
    cat ${dir}/${1}_${i}_hic_matrix.txt | perl /home/564/km6006/Scripts/inversion_paper/hic/hicmatrix2table.pl ${1} > ${dir}/${1}_${i}_hic_table.txt
done

wait
