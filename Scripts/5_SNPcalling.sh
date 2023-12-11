#!/bin/bash
#SBATCH --job-name=NanoAsseml
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time 12:00:00
#SBATCH -o %A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH --mem=15G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=


#####SNP calling using the IlluminaReads of ancestral and evolved strain using the Nanopore-based assemblies


### load modules

module load gcc12-env
module load miniconda3
module load singularity


### create environments

conda create -n trimmomatic -c bioconda trimmomatic
conda create -n fastqc -c bioconda fastqc
conda create -n bowtie2_samtools -c bioconda bowtie2 samtools
conda create -n picard -c bioconda picard
conda create -n bcftools -c bioconda bcftools



### input locations and names

assembly=/your/path/R3-I4.fasta
assembly_name=R3-I4


### input dir with the Illuminareads 

read_dir=/your/path/Reads/


### define and create outdir

out_dir=/your/path/SNP/

mkdir -p  $out_dir


##Mapping and Calling

for read1 in "$read_dir"/*_R1.fastq.sub.gz; do
	read1_stamm=${read1%.fastq.sub.gz}
	read2=${read1/_R1/_R2}
	read2_stamm=${read2%.fastq.sub.gz}
    sample=$(basename -- "${read1_stamm%_R1}")
	echo $read1_stamm 
	echo $read2_stamm
	echo $sample
	
	#Use trimmomatic (V0.39) to remove adapter
	conda activate trimmomatic
	trimmomatic PE \
		${read1} ${read2} \
		${read1_stamm}.paired.fq.gz ${read1_stamm}.unpaired.fq.gz \
		${read2_stamm}.paired.fq.gz ${read2_stamm}.unpaired.fq.gz \
		ILLUMINACLIP:TruSeq2-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
	conda deactivate
	
	#Check quality using FastQC (v0.11.5)
	conda activate fastqc
	fastqc_dir=${read_dir}fastqc/
	mkdir -p ${fastqc_dir}
	fastqc \
		${read1_stamm}.paired.fq.gz \
		${read2_stamm}.paired.fq.gz \
		--outdir $fastqc_dir
	conda deactivate
	
	#Mapping reads und sorting them using  bowtie2 and samtools
	#create BAM folder
	bam_dir=${out_dir}/BAM/
	mkdir -p $bam_dir
	
	conda activate bowtie2_samtools
	assembly_index="${assembly%.*}"
	bowtie2 -p 11 --sensitive -x $assembly_index \
		-1 ${read1_stamm}.paired.fq.gz \
		-2 ${read2_stamm}.paired.fq.gz \
		-U ${read1_stamm}.unpaired.fq.gz,${read2_stamm}.unpaired.fq.gz | samtools view -@11 -h -bS \
		| samtools sort -@11 -o ${bam_dir}${sample}_on-${assembly_name}.sorted.bam
	conda deactivate
	
	#Some formatting and marking of duplicates of generated BAM-files using picard (version 2.18.29)
	conda activate picard
	picard AddOrReplaceReadGroups \
		I=${bam_dir}${sample}_on-${assembly_name}.sorted.bam \
		O=${bam_dir}${sample}_on-${assembly_name}.sorted.RG.bam \
		RGID=${sample} RGLB=lib1 RGPL=ILLUMINA RGPU=unknown RGSM=${sample} SORT_ORDER=coordinate
	
	picard MarkDuplicates \
		I=${bam_dir}${sample}_on-${assembly_name}.sorted.RG.bam \
		O=${bam_dir}${sample}_on-${assembly_name}.sorted.RG.Dedup.bam \
		M=${bam_dir}${sample}_on-${assembly_name}.sorted.RG.Dedup.txt
		
	rm ${bam_dir}${sample}_on-${assembly_name}.sorted.RG.bam
	rm ${bam_dir}${sample}_on-${assembly_name}.sorted.RG.Dedup.txt
	conda deactivate
	
	##indexining bam file 
	conda activate samtools
	samtools index ${bam_dir}${sample}_on-${assembly_name}.sorted.RG.Dedup.bam
	conda deactivate
	
	
	#Variant calling using bcftools mpileup (version= 1.14)
	#create VCF folder
	vcf_dir=${out_dir}/VCF/
	mkdir -p $vcf_dir
	
	#calling SNPs
	conda activate bcftools
	bcftools mpileup \
		-E -C50 -Q20 -q20 \
		-f $assembly ${bam_dir}${sample}_on-${assembly_name}.sorted.RG.Dedup.bam \
		| bcftools call --ploidy 1 -vc -Ou \
		-o ${vcf_dir}${sample}_on-${assembly_name}.raw.bcf

	bcftools view \
		-v snps \
		${vcf_dir}${sample}_on-${assembly_name}.raw.bcf \
		-o ${vcf_dir}${sample}_on-${assembly_name}.raw.snps.vcf
		
	bcftools view \
		-v indels \
		${vcf_dir}${sample}_on-${assembly_name}.raw.bcf \
		-o ${vcf_dir}${sample}_on-${assembly_name}.raw.indels.vcf
	
	bgzip -c ${vcf_dir}${sample}_on-${assembly_name}.raw.snps.vcf  \
		> ${vcf_dir}${sample}_on-${assembly_name}.raw.snps.vcf.gz && tabix -p vcf ${vcf_dir}${sample}_on-${assembly_name}.raw.snps.vcf.gz
		
	bgzip -c ${vcf_dir}${sample}_on-${assembly_name}.raw.indels.vcf \
		> ${vcf_dir}${sample}_on-${assembly_name}.raw.indels.vcf.gz && tabix -p vcf ${vcf_dir}${sample}_on-${assembly_name}.raw.indels.vcf.gz
	
	
	#Filter to SNPs/Indels in with high quality in regions that are callable 
	bcftools filter \
		-o ${vcf_dir}${sample}_on-${assembly_name}.filtered.snps.vcf.gz \
		-i'QUAL>50 && DP>6 && AF1>0.7' \
		${vcf_dir}${sample}_on-${assembly_name}.raw.snps.vcf.gz
		
	bcftools filter \
		-o ${vcf_dir}${sample}_on-${assembly_name}.filtered.indels.vcf.gz  \
		-i'IDV>12 && IMF>0.7' \
		${vcf_dir}${sample}_on-${assembly_name}.raw.indels.vcf.gz
		
	tabix -p vcf \
		${vcf_dir}${sample}_on-${assembly_name}.filtered.snps.vcf.gz
		
	tabix -p vcf \
		${vcf_dir}${sample}_on-${assembly_name}.filtered.indels.vcf.gz
		
	bcftools concat -a \
		-o ${vcf_dir}${sample}_on-${assembly_name}.filtered.ALL.vcf.gz \
		${vcf_dir}${sample}_on-${assembly_name}.filtered.snps.vcf.gz \
		${vcf_dir}${sample}_on-${assembly_name}.filtered.indels.vcf.gz
		
	tabix -p vcf \
		${vcf_dir}${sample}_on-${assembly_name}.filtered.ALL.vcf.gz
		
	bcftools sort \
		${vcf_dir}${sample}_on-${assembly_name}.filtered.ALL.vcf.gz \
		-o ${vcf_dir}${sample}_on-${assembly_name}.filtered.ALL.sorted.vcf.gz
		
	conda deactivate
	
done

 