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

####Generating a Nanopore-based_Assemblie


#Load modules

module load gcc12-env
module load miniconda3


#create conda environments required for the analysis.

conda create -n bwa-mem2 -c bioconda bwa-mem2
conda create -n nanofilt -c bioconda nanofilt
conda create -n trimmomatic -c bioconda trimmomatic
conda create -n fmlrc2 -c bioconda fmlrc2  -c conda-forge rust
conda create -n canu -c bioconda canu
conda create -n flye -c bioconda flye
conda create -n racon -c biodonda racon
conda create -n medaka -c bioconda medaka
conda create -n minimap2 -c bioconda minimap2
conda create -n samtools -c bioconda samtools
conda create -n pilon -c bioconda pilon
conda create -n tapestry -c bioconda tapestry


###########################

#creating directory & defining file locations

DIR=Nanopore_Assembly/
mkdir -p $DIR


#Define the location of the Nanopore reads and the Illuminareads used for correction

#Nanopore reads

NReads=/your/path/${strain_name}.nanopore.fastq.gz


#Illuminareads

IReadsR1=/your/path/${strain_name}_R1.fastq.gz
IReadsR2=/your/path/${strain_name}_R2.fastq.gz


#Strain-name

strain_name=R3-A


#location of the reformat.sh from bbmap

reformat=/your/path/bbmap/reformat.sh


##############################


###filtering Nanopore reads to those excieding 5000 bases using NanoFilt (v2.3.0)

conda activate nanofilt
gunzip -c $NReads | NanoFilt  -l 5000 | gzip > ${DIR}${strain_name}.5000.fq.gz
conda deactivate


#Trimming of Illumina PE reads (trimmed using Trimmomatic V0.39:

conda activate trimmomatic
trimmomatic PE $IReadsR1 $IReadsR2 \
	${DIR}${strain_name}_paired_R1.fq ${DIR}${strain_name}_unpaired_R1.fq \
	${DIR}${strain_name}_paired_R2.fq ${DIR}${strain_name}_unpaired_R2.fq \
	ILLUMINACLIP:TruSeq2-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
conda deactivate


#interleave the trimmed paired reads

conda activate picard
$reformat  \
	in1=${DIR}${strain_name}_paired_R1.fq \
	in2=${DIR}${strain_name}_paired_R2.fq  \
	out=${DIR}${strain_name}_trimmed.fq
gzip ${DIR}${strain_name}_trimmed.fq > ${DIR}${strain_name}_trimmed.fq.gz 
conda deactivate 


#Correction of Nanopore reads by trimmed Illumina reads using FMLRC 2 (0.1.4)

conda activate fmlrc2
gunzip -c ${DIR}${strain_name}_trimmed.fq.gz  | awk 'NR % 4 == 2' | sort | tr NT TN | ropebwt2 -LR | tr NT TN | fmlrc2-convert ${DIR}comp_msbwt.npy
fmlrc2 \
	${DIR}comp_msbwt.npy ${DIR}${strain_name}.5000.fq.gz \
	${DIR}${strain_name}.5000.corrected.fq
conda deactivate


#Further trimming of corrected reads using Canu (v2.1.1)

conda activate canu
canu -trim -p ${strain_name} \
	-d ${strain_name} \
	genomeSize=38m  minInputCoverage=0  \
	-nanopore -corrected ${DIR}${strain_name}.5000.corrected.fq
conda deactivate


#Assembly of corrected and trimmed Nanopore reads using flye (V 2.8.3) and two rounds of polishing

conda activate flye
flye \
	--nano-corr ${Dir}${strain_name}.5000.corrected.trimmed.fq.gz \
	--out-dir $DIR --threads 11 --trestle -i 2
conda deactivate
mv ${DIR}assembly.fasta ${DIR}${strain_name}.fasta


#Assembly polishing by uncorrected Nanopore reads using Racon (V1.4.20) twice. 

cp  ${DIR}${strain_name}.fasta ${DIR}${strain_name}.rancon1.fasta
for i in {1..2}; do
	#Mapping the uncorrected Nanopore reads onto the assembly
	conda activate minimap2
	minimap2 -ax map-ont -t 6 -L ${DIR}${strain_name}.rancon${i}.fasta ${DIR}${strain_name}.5000.fq.gz > ${DIR}5000Reads_on_${strain_name}.rancon${i}.sam
	conda deactivate
	
	#Using rancon to polish the assembly (twice)
	conda activate racon
	e=$((i + 1))
	racon -t 6 --no-trimming ${DIR}${strain_name}.5000.fq.gz \
		${DIR}5000Reads_on_${strain_name}.sam \
		${DIR}${strain_name}.rancon${i}.fasta > ${DIR}${strain_name}.rancon${e}.fasta	
	conda deactivate


#Assembly polishing by uncorrected Nanopore reads using Medaka (V 1.4.3): 

conda activate medaka
medaka_consensus -i ${DIR}${strain_name}.5000.fq.gz \
	-d ${DIR}${strain_name}.rancon3.fasta \
	-o ${DIR} \
	-t 4 -m  r941_min_high_g360
conda deactivate

cp ${DIR}consensus.fasta ${DIR}${strain_name}.rancon3.medaka.pilon1.fasta


#Assembly polishing by Illumina reads using Pilon (V 1.24) after mapping with BWA-mem2 (V 2.2.1) and samtools (1.12) (four times):

for i in {1..4}; do
	#Mapping Illumina reads on Assembly using bwa-mem2
	conda activate bwa-mem2
	bwa-mem2 index -p ${DIR}${strain_name}.rancon3.medaka.pilon${i} ${DIR}${strain_name}.rancon3.medaka.pilon${i}.fasta 
	bwa-mem2 mem -t 11 ${DIR}${strain_name}.rancon3.medaka.pilon${i} ${DIR}${strain_name}_trimmed.fq.gz  > ${DIR}${strain_name}_trimmed_on-${strain_name}.rancon3.medaka.pilon${i}.sam 
	conda deactivate

	conda activate samtools
	samtools sort ${DIR}${strain_name}_trimmed_on-${strain_name}.rancon3.medaka.pilon${i}.sam  -o ${DIR}${strain_name}_trimmed_on-${strain_name}.rancon3.medaka.pilon${i}.bam
	samtools index ${DIR}${strain_name}_trimmed_on-${strain_name}.rancon3.medaka.pilon${i}.bam
	conda deactivate

	#use pilon to improve the assmebly
	conda activate pilon
	e=$((i + 1))
	java -Xmx16G -jar pilon.jar --genome ${DIR}${strain_name}.rancon3.medaka.pilon${i}.fasta  --bam ${DIR}${strain_name}_trimmed_on-${strain_name}.rancon3.medaka.pilon${i}.bam --output ${DIR}${strain_name}.rancon3.medaka.pilon${e} --outdir ${DIR} --changes
	#remove the "_pilon" from the fasta file
	sed -i 's/_pilon//g' "${DIR}${strain_name}.rancon3.medaka.pilon${e}.fasta"
done


# The final assembly was analyzed using Tapestry (V 1.0.0) using the telomeric sequence  TTAGGG:

conda activate Tapestry
weave -a ${DIR}${strain_name}.rancon3.medaka.pilon5.fasta \
	-r ${DIR}${strain_name}.5000.fq.gz \
	-t TTAGGG -o ${SIR} -c 10
 
#visually inspect te tapestry report and indicate in the csv which contigs should be removed (coverage lower than 30 or higher than60)

clean -a ${DIR}${strain_name}.rancon3.medaka.pilon5.fastaa -c ${DIR}${strain_name}.rancon3.medaka.pilon5.filtered.csv -o ${DIR}${strain_name}.rancon3.medaka.pilon5.clean.fasta

conda deactivate
