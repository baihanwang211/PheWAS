cd /home/baihanw/psych/data/pgs

trait='mdd'

population1='EUR'

population2='EAS'

# SLURM_ARRAY_TASK_ID=22

plink --bfile /ckb/prs/data/chr${SLURM_ARRAY_TASK_ID}_rs \
--score ./${trait}_${population1}_${population2}/${trait}_META_pst_eff_a1_b0.5_phiauto_chr${SLURM_ARRAY_TASK_ID}.txt 2 4 6 sum \
--out ./${trait}_${population1}_${population2}/${trait}_META_chr${SLURM_ARRAY_TASK_ID}


