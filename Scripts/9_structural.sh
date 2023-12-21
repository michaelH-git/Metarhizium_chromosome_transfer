#!/bin/bash
#SBATCH --job-name=struct
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time 12:00:00
#SBATCH -o %A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH --mem=15G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=

module load gcc12-env
module load miniconda3

#create environments

conda create -n sniffels -c bioconda sniffles
conda creare -n samtools -c bioconda samtools
conda create -n bcftools -c bioconda bcftools

##define input

Ref=/your/path/R3I4.fasta
NanoR=/your/path/R3-I4.nanopore.fastq.gz

##define output
DIR=/your/path/Structural/
mkdir -p $DIR
cd $DIR

##generate names

refname=$(basename -- "${Ref%.fasta}")
readname=$(basename -- "${NanoR%.nanopore.fastq.gz}")


#map Nanopore reads using minimap2 (version 2.26-r1175)

conda activate minimap2 

minimap2 -L -ax map-ont $Ref  $NanoR > ${DIR}${readname}_on-${refname}.sam

conda deactivate


##some reformating and sorting using samtools (v 1.17)

conda activate samtools

samtools view -bS ${DIR}${readname}_on-${refname}.sam > ${DIR}${readname}_on-${refname}.bam
samtools sort ${DIR}${readname}_on-${refname}.bam -o ${DIR}${readname}_on-${refname}.sorted.bam
samtools index ${DIR}${readname}_on-${refname}.sorted.bam

conda deactivate

##calling larger strutural variants using sniffles (version 2.2)


conda activate sniffels

sniffles -i ${DIR}${readname}_on-${refname}.sorted.bam -v ${DIR}${readname}_on-${refname}.sorted.vcf

conda deactivate

##filtering VCF to those strutural variants that have an AF of > 0.7 (haploid genome) using bcftools (v 1.17)

conda activate bcftools

bcftools filter -o ${DIR}${readname}_on-${refname}.sorted.filtered.vcf  -i'AF>0.7' ${DIR}${readname}_on-${refname}.sorted.vcf

conda deactivate
