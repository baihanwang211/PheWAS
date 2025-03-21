cd /home/baihanw/psych/data/pgs

trait='mdd'

population1='EUR'
population2='EAS'

n_gwas1=1309348
n_gwas2=51910

mkdir ${trait}_${population1}_${population2}

# SLURM_ARRAY_TASK_ID=22

python /home/baihanw/software/PRScsx/PRScsx.py \
--ref_dir=/home/baihanw/reference/1kg \
--bim_prefix=/ckb/prs/data/chr${SLURM_ARRAY_TASK_ID}_rs \
--sst_file=/home/baihanw/psych/data/gwas/harmonised/${trait}_${population1}_prscs.txt,/home/baihanw/psych/data/gwas/harmonised/${trait}_${population2}_prscs.txt \
--n_gwas=${n_gwas1},${n_gwas2} \
--pop=${population1},${population2} \
--chrom=${SLURM_ARRAY_TASK_ID} \
--out_dir=./${trait}_${population1}_${population2}/ \
--out_name=${trait} \
--meta=TRUE