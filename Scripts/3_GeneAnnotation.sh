#!/bin/bash
#SBATCH --job-name=GeneAnno
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time 12:00:00
#SBATCH -o %A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH --mem=15G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=

##load modules

module load gcc12-env
module load miniconda3
module load singularity


##creating environments

conda create -n gffread -c bioconda gffread
conda create -n busco -c bioconda  busco


#run braker as a singulrity container

SINGULARITY_IMAGE=/your/path/braker2.sif
GENOME=/your/path/R3-I4.softmasked_TE.fasta
PROT_SEQ=/your/path/odb10_fungi_fasta.tar.gz
SPECIES="R3-I4"
AUG_PATH=/your/path/augustus_config/


singularity exec --bind /your_path/:/your_path/  ${SINGULARITY_IMAGE} braker.pl \
       --cores 4 \
        --species=${SPECIES} \
		--useexisting \
       --genome=${GENOME} \
        --prot_seq=${PROT_SEQ} \
		--softmasking \
		--fungus \
		--AUGUSTUS_CONFIG_PATH=$AUG_PATH
		

# Use of the Augustus Gene prediction: augustus.hints.gtf

conda create -n gffread -c bioconda gffread

conda activate gffread
gffread -o augustus.hints.gff3  /Braker/augustus.hints.gtf

#extract protein sequences using gffread

gffread -g ${GENOME} -y ${GENOME}.proteins.fasta augustus.hints.gff3

conda deactivate

#Busco Analysis (V5.2.2)

conda create -n busco -c bioconda  busco

busco -i ${GENOME}.proteins.fasta  -l sordariomycetes_odb10 -o buscoR3-I4 -m protein -c 10 -f --offline --download_path /path/to/datasets

