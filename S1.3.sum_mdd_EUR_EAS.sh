trait='mdd'

population1='EUR'

population2='EAS'

cd /home/baihanw/psych/data/pgs/${trait}_${population1}_${population2}

awk 'BEGIN {OFS="\t"} FNR>1 {sum[$2] += $6} END {for (i in sum) print i, sum[i]}' ${trait}_META_chr*.profile > ${trait}_${population1}_${population2}.profile