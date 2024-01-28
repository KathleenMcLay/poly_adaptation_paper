### Calculate Pi, Dxy and Fst between populations, for all chromosomes

source ~/miniconda3/etc/profile.d/conda.sh
conda activate bio

module load bcftools 

FST_IN="/g/data/ht96/McLay_UQ/avneet_paper/pa_vo_sf5_final.vcf.gz"
D01="/home/564/km6006/Scripts/avneet_paper/D01.txt"
H01="/home/564/km6006/Scripts/avneet_paper/H01.txt"

PIXY_IN="/g/data/ht96/McLay_UQ/avneet_paper/pa_as_sf5_final.vcf.gz"
PIXY_DIR="/g/data/ht96/McLay_UQ/avneet_paper/pvp_popgen"
POPS="/home/564/km6006/Scripts/avneet_paper/D01H01.txt"

# Calculate Pi and Dxy between populations 
tabix ${PIXY_IN}
pixy --stats pi dxy --vcf ${PIXY_IN} --populations ${POPS} --window_size 10000 --output_prefix D01H01 --output_folder ${PIXY_DIR} --n_cores 12

# Calculate Fst between populations 
/home/564/km6006/bin/vcftools \
    --gzvcf ${FST_IN} \
    --weir-fst-pop ${D01} \
    --weir-fst-pop ${H01} \
    --fst-window-size 10000 \
    --fst-window-step 10000 \
    --out /g/data/ht96/McLay_UQ/avneet_paper/pvp_popgen/D01H01
