#!/bin/bash
#SBATCH --job-name=Phasing
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time 12:00:00
#SBATCH -o %A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH --mem=15G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=


#######Load modules

module load gcc12-env
module load miniconda3
module load singularity


######Create environments

conda create -n minimap2 -c bioconda minimap2
conda create -n whatshap -c bioconda whatshap
conda create -n bcftools -c bioconda bcftools
conda create -n picard -c bioconda picard



#define file locations and names

pacbioR=/your/path/5816_A.run725.hifi_reads.sub.fastq.gz
assembly=your/path/R3-I4.fasta
strain_name=R3-I4
vcf_file=your/path/Mguizhouense977_on-R3-I.vcf ##VCF file generated using M.guizhouense illumina files mapped onto R3-I4 assembly 


##create output directory

outdir=/work_beegfs/sunbo356/data/Metarhizium_2023/Revision_horizontalTransfer/phasing/
mkdir -p $outdir

################Phasing of SNPs and Small Indeels##########################################

#Mapping of pacbio reads onto the on R3-I4 assembly using minimap2 (version 2.24-r1122), samtools (version 1.3.1):

conda activate minimap2
minimap2 -L -ax map-hifi  \
	$assembly  $pacbioR \
	> ${outdir}MguizhouenseARSEF977_on_${strain_name}.sam
conda deactivate

#####reformating and sorting using samtools and picard (for changing the readgroup information to match the sample information in the VCF file)

conda activate samtools
samtools view -S -b ${outdir}MguizhouenseARSEF977_on_${strain_name}.sam  > ${outdir}MguizhouenseARSEF977_on_${strain_name}.bam
samtools sort ${outdir}MguizhouenseARSEF977_on_${strain_name}.bam -o ${outdir}MguizhouenseARSEF977_on_${strain_name}.sorted.bam
samtools index ${outdir}MguizhouenseARSEF977_on_${strain_name}.sorted.bam
rm ${outdir}MguizhouenseARSEF977_on_${strain_name}.sam ${outdir}MguizhouenseARSEF977_on_${strain_name}.bam
conda deactivate

##changing of readgroup information

conda activate picard
picard AddOrReplaceReadGroups INPUT=${outdir}MguizhouenseARSEF977_on_${strain_name}.sorted.bam  \
	OUTPUT=${outdir}MguizhouenseARSEF977_on_${strain_name}.sorted.RG.bam  \
	RGID=Mguizhouense977 RGLB=lib1 RGPL=ILLUMINA RGPU=unknown RGSM=Mguizhouense977 SORT_ORDER=coordinate
conda deactivate

conda activate samtools
samtools index ${outdir}MguizhouenseARSEF977_on_${strain_name}.sorted.RG.bam 
conda deactivate

#Filter of vcf produced by published illumina data of M.guizhouense ARSEF977 mapped on R3-I4 

conda activate bcftools

vcf_name=$(basename -- "${vcf_file%.vcf}")
bcftools filter -o ${outdir}${vcf_name}.filtered.vcf  -i 'QUAL>50 && DP>10' $vcf_file

conda deactivate

#Using WhatsHap to phase the SNPS and small InDels (version1.6)
conda activate whatshap
whatshap phase -o ${outdir}${vcf_name}.filtered.phased.vcf \
		--reference=$assembly \
		${outdir}${vcf_name}.filtered.vcf  \
		${outdir}MguizhouenseARSEF977_on_${strain_name}.sorted.RG.bam 
conda deactivate

#################################

