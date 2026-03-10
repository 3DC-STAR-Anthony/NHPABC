#############################################################################
#################################----Liftover----############################
#############################################################################
# Convert monkey bed files to hg38

conda activate liftover

filedir=~/01.Project/NHPABC/Figure6/06.Conservation/02.NHPABC_cCRE_bed
outdir=~/01.Project/NHPABC/Figure6/06.Conservation/02.NHPABC_cCRE_bed/hg38

# All subtypes except DG Ex
for f in `ls ~/01.Project/NHPABC/Figure6/06.Conservation/02.NHPABC_cCRE_bed/*.bed`;
  do
   n=`basename $f | sed 's/.bed//g'`
   liftOver \
     ${filedir}/${n}.bed \
     ~/GCF_037993035.2ToHg38.over.chain.gz \
     ${outdir}/${n}.MacTohg38.bed \
     ${outdir}/${n}.unlifted.bed
    done

# Convert DG Ex to hg19
filedir=~/01.Project/NHPABC/Figure6/06.Conservation/02.NHPABC_cCRE_bed/hg38/
outdir=~/01.Project/NHPABC/Figure6/06.Conservation/02.NHPABC_cCRE_bed/hg19

   liftOver \
     ${filedir}/Hf_DG_Ex.Ex.nhpabc.MacTohg38.bed \
     ~/hg38ToHg19.over.chain.gz \
     ${outdir}/Hf_DG_Ex.Ex.nhpabc.hg38Tohg19.bed \
     ${outdir}/Hf_DG_Ex.Ex.nhpabc.unlifted.bed

#############################################################################
#################################----Bedtool----#############################
#############################################################################
# Use bedtools to find overlaps between monkey and human peaks for corresponding cell types
cd "~/01.Project/NHPABC/Figure6/06.Conservation/03.NHPABC_cCRE_conservation/"

# Sort and merge human peaks
for f in `ls Hip.WXQ_human_brain.bed`
do
name=`echo $f | sed 's/.bed//g'`
sort -k1,1 -k2,2n $f > ${name}_sort.bed
bedtools merge -i ${name}_sort.bed > ${name}_unique.bed
done

cd ~/01.Project/NHPABC/Figure6/06.Conservation/03.NHPABC_cCRE_conservation/

# Excitatory neurons
for f in `ls ~/01.Project/NHPABC/Figure6/06.Conservation/02.NHPABC_cCRE_bed/hg38/*.Ex.nhpabc.MacTohg38.bed`
do
n=`basename $f | sed 's/.Ex.nhpabc.MacTohg38.bed//g'`
bedtools intersect \
 -a $f \
 -b ~/01.Project/NHPABC/Figure6/06.Conservation/01.human_maincelltype_cCRE/Ex.renbing_human_brain_unique.bed \
 -wa -wb > ${n}_intersect.bedtool.bed
done

# Inhibitory neurons
for f in `ls ~/01.Project/NHPABC/Figure6/06.Conservation/02.NHPABC_cCRE_bed/hg38/*.Inh.nhpabc.MacTohg38.bed`
do
n=`basename $f | sed 's/.Inh.nhpabc.MacTohg38.bed//g'`
bedtools intersect \
 -a $f \
 -b ~/01.Project/NHPABC/Figure6/06.Conservation/01.human_maincelltype_cCRE/In.renbing_human_brain_unique.bed \
 -wa -wb > ${n}_intersect.bedtool.bed
done

# Astrocytes
for f in `ls ~/01.Project/NHPABC/Figure6/06.Conservation/02.NHPABC_cCRE_bed/hg38/*.Ast.nhpabc.MacTohg38.bed`
do
n=`basename $f | sed 's/.Ast.nhpabc.MacTohg38.bed//g'`
bedtools intersect \
 -a $f \
 -b ~/01.Project/NHPABC/Figure6/06.Conservation/01.human_maincelltype_cCRE/Astrocyte.renbing_human_brain_unique.bed \
 -wa -wb > ${n}_intersect.bedtool.bed
done

# Microglia
for f in `ls ~/01.Project/NHPABC/Figure6/06.Conservation/02.NHPABC_cCRE_bed/hg38/*.Mic.nhpabc.MacTohg38.bed`
do
n=`basename $f | sed 's/.Mic.nhpabc.MacTohg38.bed//g'`
bedtools intersect \
 -a $f \
 -b ~/01.Project/NHPABC/Figure6/06.Conservation/01.human_maincelltype_cCRE/Microglia.renbing_human_brain_unique.bed \
 -wa -wb > ${n}_intersect.bedtool.bed
done

# OPC (Oligodendrocyte Precursor Cells)
for f in `ls ~/01.Project/NHPABC/Figure6/06.Conservation/02.NHPABC_cCRE_bed/hg38/*.OPC.nhpabc.MacTohg38.bed`
do
n=`basename $f | sed 's/.OPC.nhpabc.MacTohg38.bed//g'`
bedtools intersect \
 -a $f \
 -b ~/01.Project/NHPABC/Figure6/06.Conservation/01.human_maincelltype_cCRE/OPC.renbing_human_brain_unique.bed \
 -wa -wb > ${n}_intersect.bedtool.bed
done

# Oligodendrocytes
for f in `ls ~/01.Project/NHPABC/Figure6/06.Conservation/02.NHPABC_cCRE_bed/hg38/*.ODC.nhpabc.MacTohg38.bed`
do
n=`basename $f | sed 's/.ODC.nhpabc.MacTohg38.bed//g'`
bedtools intersect \
 -a $f \
 -b ~/01.Project/NHPABC/Figure6/06.Conservation/01.human_maincelltype_cCRE/Oligodendrocyte.renbing_human_brain_unique.bed \
 -wa -wb > ${n}_intersect.bedtool.bed
done

# VS (Vascular cells)
for f in `ls ~/01.Project/NHPABC/Figure6/06.Conservation/02.NHPABC_cCRE_bed/hg38/*.VS.nhpabc.MacTohg38.bed`
do
n=`basename $f | sed 's/.VS.nhpabc.MacTohg38.bed//g'`
bedtools intersect \
 -a $f \
 -b ~/01.Project/NHPABC/Figure6/06.Conservation/01.human_maincelltype_cCRE/VS.renbing_human_brain_unique.bed \
 -wa -wb > ${n}_intersect.bedtool.bed
done

# DG Ex neurons (special case, using hippocampus data)
for f in `ls ~/01.Project/NHPABC/Figure6/06.Conservation/02.NHPABC_cCRE_bed/hg19/Hf_DG_Ex.Ex.nhpabc.hg38Tohg19.bed`
do
n=`basename $f | sed 's/.Ex.nhpabc.hg38Tohg19.bed//g'`
bedtools intersect \
 -a $f \
 -b ~/01.Project/NHPABC/Figure6/06.Conservation/01.human_maincelltype_cCRE/Hip.WXQ_human_brain_unique.bed \
 -wa -wb > ${n}_intersect.bedtool.bed
done
