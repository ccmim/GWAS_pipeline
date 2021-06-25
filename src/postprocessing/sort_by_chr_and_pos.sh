FILE=$1
TMPFILE=${FILE}.tmp

awk 'NR == 1; NR > 1 {print $0 | "sort -k1 -n -k3 -n"}' $FILE > $TMPFILE

mv $TMPFILE $FILE
