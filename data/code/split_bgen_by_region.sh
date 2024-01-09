
CHROMOSOME=$1
STARTPOS=$2
ENDPOS=$3
BGEN_DIR="/home/home01/scrb/nobackup/data/ukbb/genotypes/imputed"
BGEN=${BGEN_DIR}/full/ukb22828_c${CHROMOSOME}_b0_v3.bgen
BGI=${BGEN_DIR}/full/ukb_imp_chr${CHROMOSOME}_v3.bgen.bgi
OFILE=${BGEN_DIR}/by_region/ukb_imp_chr${CHROMOSOME}_${STARTPOS}_${ENDPOS}.bgen


# bgenix -g $BGEN -i $BGI -incl-range $CHROMOSOME:$STARTPOS-$ENDPOS > $OFILE 
