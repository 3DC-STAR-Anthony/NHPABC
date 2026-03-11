#########################################################################################################
###############--------------------------LD Score Estimation Tutorial--------------------------##########
#########################################################################################################
conda activate liftover

# Compute LD Score
# Step 1: Transform single cell ATAC-seq peaks bed files from Mac5 to hg19
# Note: Cell types are numbered for convenience (e.g., subtype1.bed corresponds to Ast_PFC.bed)
filedir=~/01.cCRE/01.raw/
outdir=~/01.cCRE/02.mac5Tohg38/
    
for i in {1..166};
do
    liftOver \
    ${filedir}/subtype${i}.bed \
    GCF_037993035.2ToHg38.over.chain.gz \
    ${outdir}${i}.MacTohg38.bed \
    ${outdir}${i}.unlifted.bed
done

filedir=~/01.cCRE/02.mac5Tohg38/
outdir=~/01.cCRE/03.hg38Tohg19/

for i in {1..166};
do
    liftOver \
    ${filedir}/${i}.MacTohg38.bed \
    hg38ToHg19.over.chain.gz \
    ${outdir}${i}.hg38Tohg19.bed \
    ${outdir}${i}.unlifted.bed
done
    
# Step 2: Create annotation files
source activate ldsc

outdir=~/01.cCRE/03.hg38Tohg19/
for i in {1..166};
do
    for j in {1..22}
    do 
        python make_annot.py \
        --bed-file ${outdir}${i}.hg38Tohg19.bed \
        --bimfile ~01.1000G_EUR_Phase3_plink/1000G.EUR.QC.${j}.bim \
        --annot-file ~/02.LDSC_annot_file/${i}.${j}.annot.gz
    done
done

# Step 3: Compute LD scores with annotation files
for i in {1..166};
do
    for j in {1..22}
    do 
        python ldsc.py \
        --l2 \
        --bfile ~01.1000G_EUR_Phase3_plink/1000G.EUR.QC.${j} \
        --ld-wind-cm 1 \
        --annot ~/02.LDSC_annot_file/${i}.${j}.annot.gz \
        --thin-annot \
        --print-snps ~02.hapmap3_snps/hm.${j}.snp \
        --out ~/03.LDSC_LD_scores/${i}.${j}
    done
done

# Step 4: Execute enrichment analysis using Perl script
#!/usr/bin/perl -w
use strict;
open IN, "~/phenotype.tsv" || die "$!";
while(<IN>){
    next if($. == 1);
    chomp;
    my($a, $b) = split(/\s+/, $_);
    `python ldsc.py --h2-cts $b --ref-ld-chr ~03.1000G/1000G_EUR_Phase3_baseline/baseline. --out ~/04.Enrichment/$a --ref-ld-chr-cts ~/subtype_cCREs_ID.ldcts --w-ld-chr ~03.1000G/weights_hm3_no_hla/weights.`;
}

perl ~/Enrichment_ID.pl
