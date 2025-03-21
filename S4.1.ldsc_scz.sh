# Define an array of phenos
phenos=("bmi" "height" "ea" "smkinit" "cigday")

pheno="${phenos[$SLURM_ARRAY_TASK_ID-1]}"

# Run ldsc for EAS
munge_sumstats.py \
--sumstats "${pheno}_EAS.hapmap3.ldsc.tsv.gz" \
--snp SNP \
--N-col N \
--a1 A1 \
--a2 A2 \
--p P \
--signed-sumstats Beta,0 \
--out "${pheno}_EAS" \
--merge-alleles w_hm3.snplist \
--chunksize 500000

ldsc.py \
--rg "scz_EAS.sumstats.gz","${pheno}_EAS.sumstats.gz" \
--ref-ld-chr eas_ldscores/ \
--w-ld-chr eas_ldscores/ \
--out "scz_${pheno}_EAS"

# Run ldsc for EUR
munge_sumstats.py \
--sumstats "${pheno}_EUR.hapmap3.ldsc.tsv.gz" \
--snp SNP \
--N-col N \
--a1 A1 \
--a2 A2 \
--p P \
--signed-sumstats Beta,0 \
--out "${pheno}_EUR" \
--merge-alleles w_hm3.snplist \
--chunksize 500000

ldsc.py \
--rg "scz_EUR.sumstats.gz","${pheno}_EUR.sumstats.gz" \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out "scz_${pheno}_EUR"

# between scz and mdd
pheno="mdd"

ldsc.py \
--rg "scz_EAS.sumstats.gz","${pheno}_EAS.sumstats.gz" \
--ref-ld-chr eas_ldscores/ \
--w-ld-chr eas_ldscores/ \
--out "scz_${pheno}_EAS"

ldsc.py \
--rg "scz_EUR.sumstats.gz","${pheno}_EUR.sumstats.gz" \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out "scz_${pheno}_EUR"
