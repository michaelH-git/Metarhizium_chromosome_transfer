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


###########Analysis of published Assemblies & Reads

#Load modules

module load gcc12-env
module load miniconda3
module load singularity


# create environments

conda create -n mummer -c bioconda mummer
conda create -n minimap2 -c bioconda minimap2
conda create -n bedtools -c bioconda bedtools
conda create -n trimmomatic -c bioconda trimmomatic
conda create -n fastqc -c bioconda fastqc
conda create -n bowtie2_samtools -c bioconda bowtie2 samtools
conda create -n samtools -c bioconda samtools
conda create -n picard -c bioconda picard
conda create -n bcftools -c bioconda bcftools
conda create -n whatshap -c whatshap


#Folder with the published assemblies
publa=/your/path/publishedAssemblies/

#Folder with published reads
publreads=/your/path/publishedReads/

#location of the reformat.sh of the bbmap package
reformat=/zfshome/sunbo356/bbmap/reformat.sh

#input locations and names
assembly=your/path/R3-I4.fasta
assembly_name=R3-I4
outdir=your/path/Comparison/

mkdir -p  $outdir

######Analysis of published Assemblies



conda activate mummer

for asm in "$publa"/*.fna; do
	#use nucmer to make alignment
    filename=$(basename -- "$asm")
    filename_no_ext="${filename%.*}"
	findir=${outdir}/nucmer/
	mkdir -p $findir
	out=${findir}${filename_no_ext}
	nucmer --maxmatch -p $out  $assembly $asm

	#Filter to those matches of min 1000bp length with a minimum identity of 90%
	delta-filter -i 90 -l 1000 ${out}.delta > ${out}.i90.l1000.delta

	#generate bed file for coverage analysis 
	show-coords ${out}.i90.l1000.delta > ${out}.i90.l1000.coords
	awk -v OFS='\t' '{if ($10 >= 90 && $7 >= 90) print $12,$1-1,$2,$13-$10}'  ${out}.i90.l1000.coords > ${out}.i90.l1000.bed
	
	conda activate bedtools
	sortBed -i ${out}.i90.l1000.bed > ${out}.i90.l1000.sorted.bed 
	bedtools merge -i  ${out}.i90.l1000.sorted.bed > ${out}.i90.l1000.sorted.merged.bed
	conda deactivate
	
	#Determining the number of SNPs in regions with coverage
	show-snps -T -Clr  ${out}.i90.l1000.delta > ${out}.i90.l1000.snps
	all2vcf mummer --snps ${out}.i90.l1000.snps --reference $assembly --type SNP  --input-header --output-header > ${out}.i90.l1000.vcf
	conda activate minimap2
	paftools.js delta2paf ${out}.i90.l1000.delta  > ${out}.i90.l1000.paf
	conda deactivate
done
conda deactivate
 
###############################################
 
 
########Analysis of published Reads



for reads in "$publreads"/*.fastq; do
    rname=$(basename -- "$reads")
    rname_no_ext="${rname%.*}"
	readdir=${outdir}/sra/
	mkdir -p $readdir
	outreads=${readdir}${rname_no_ext}
	
	#Deinterleaving the using the reformat.sh utility of bbmap (version 39.01)
	conda activate picard
	$reformat in=$reads \
	out1=${outreads}_1.fastq \
	out2=${outreads}_2.fastq
	conda deactivate
	
	#Use trimmomatic (V0.39) to remove adapter
	conda activate trimmomatic
	trimmomatic PE \
	${outreads}_1.fastq ${outreads}_2.fastq \
	${outreads}_1_paired.fq.gz ${outreads}_1_unpaired.fq.gz \
	${outreads}_2_paired.fq.gz ${outreads}_2_unpaired.fq.gz \
	ILLUMINACLIP:TruSeq2-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
	conda deactivate
	
	#Check quality using FastQC (v0.11.5)
	conda activate fastqc
	fastqc_dir=${readdir}fastqc/
	mkdir -p ${fastqc_dir}
	fastqc ${outreads}_1_paired.fq.gz ${outreads}_2_paired.fq.gz --outdir $fastqc_dir
	conda deactivate
	
	#Mapping reads und sorting them using  bowtie2 and samtools
	conda activate bowtie2_samtools
	assembly_index="${assembly%.*}"
	bowtie2-build $assembly $assembly_index
	bowtie2 -p 11 --fast -x $assembly_index \
	-1 ${outreads}_1_paired.fq.gz \
	-2 ${outreads}_2_paired.fq.gz \
	-U  ${outreads}_1_unpaired.fq.gz,${outreads}_2_unpaired.fq.gz | samtools view -@11 -h -bS \
	| samtools sort -@11 -o ${readdir}${rname_no_ext}_on-${assembly_name}.sorted.bam
	conda deactivate
	
	#Some formatting and marking of duplicates of generated BAM-files using picard (version 2.18.29)
	conda activate picard
	picard AddOrReplaceReadGroups \
	I=${readdir}${rname_no_ext}_on-${assembly_name}.sorted.bam \
	O=${readdir}${rname_no_ext}_on-${assembly_name}.sorted.RG.bam \
	RGID=${rname_no_ext} RGLB=lib1 RGPL=ILLUMINA RGPU=unknown RGSM=${rname_no_ext} SORT_ORDER=coordinate
	
	picard MarkDuplicates \
	I=${readdir}${rname_no_ext}_on-${assembly_name}.sorted.RG.bam  \
	O=${readdir}${rname_no_ext}_on-${assembly_name}.sorted.RG.Dedup.bam M=${readdir}${rname_no_ext}_on-${assembly_name}.sorted_RG_Dedup.txt
	rm ${readdir}${rname_no_ext}_on-${assembly_name}.sorted.RG.bam
	rm ${readdir}${rname_no_ext}_on-${assembly_name}.sorted_RG_Dedup.txt

	#Variant calling using bcftools mpileup (version= 1.14)
	conda activate bcftools
	bcftools mpileup -A -f $assembly ${readdir}${rname_no_ext}_on-${assembly_name}.sorted.RG.Dedup.bam | bcftools call --ploidy 1 -vc -Ou -o ${readdir}${rname_no_ext}_on-${assembly_name}.raw.bcf   ##please note that M.guizhouense ARSEF977 the polidy was set to two to account for the disomic chrA 
	
	bcftools view ${readdir}${rname_no_ext}_on-${assembly_name}.raw.bcf -o ${readdir}${rname_no_ext}_on-${assembly_name}.raw.vcf 
	bgzip -c ${readdir}${rname_no_ext}_on-${assembly_name}.raw.vcf  > ${readdir}${rname_no_ext}_on-${assembly_name}.raw.vcf.gz && tabix -p vcf ${readdir}${rname_no_ext}_on-${assembly_name}.raw.vcf.gz 
	
	#Filter to SNPs in with high quality in regions that are callable 
	bcftools filter -o ${readdir}${rname_no_ext}_on-${assembly_name}.Q50DP6.vcf.gz  -i'QUAL>50 && DP>6' ${readdir}${rname_no_ext}_on-${assembly_name}.raw.vcf.gz 
	bcftools view -v snps ${readdir}${rname_no_ext}_on-${assembly_name}.Q50DP6.vcf.gz  > ${readdir}${rname_no_ext}_on-${assembly_name}.Q50DP6.snps.vcf
	conda deactivate
done

############################################################