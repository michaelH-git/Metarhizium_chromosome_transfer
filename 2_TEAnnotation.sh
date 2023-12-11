#!/bin/bash
#SBATCH --job-name=repet
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --partition=standard 

##### Transposable Element Annotation using REPET3 
###IMPORTANT: You need access to a MySQL server. REPET creates an SQL table for each project. These SQL table need to be removed before rerunning the project again - otherwise the REPET pipeline fill fail.

##make directory for 

mkdir R1A
cd R1A


#copying the setEnv.sh and set the environment to working directory (specifc the Repet host and dependencies - adjust for your local system) 

cp /home/setEnv.sh ./


#Sourcing the environment (adjust setEnv to your system)

. setEnv.sh


#create new folder for TEannot

mkdir TEannot
cd TEannot


#Copy the configuration files TEdenovo.cfg TEannot.cfg into your working directory and the link needed files (genome fasta). Adappt TEdenovo to your system.

cp /home/TEannot.cfg ./

#linking assembly into workfolder
ln s /home/R3-A.assembly.fasta ./R1A.fa

#Running the denovo pipeline

launch_TEdenovo.py -P R1A -f MCL >& denovo.txt &

#Output: R1A_sim_denovoLibTEs_filtered.fa  (filtered consensus sequences of TEs)

#Second, this filtered consensus sequences of TEs identified by the denovo annotation were used to annote the TEs in the genome. Here, two rounds of TEannot were conducted. 

TEannot.py -P R1A_annot -C TEannot.cfg -S 1
TEannot.py -P  R1A_annot  -C TEannot.cfg -S 2 -a BLR  -v 2
TEannot.py -P  R1A_annot  -C TEannot.cfg -S 2 -a RM -v 2
TEannot.py -P  R1A_annot  -C TEannot.cfg -S 2 -a CEN -v 2
TEannot.py -P  R1A_annot  -C TEannot.cfg -S 2 -a BLR -r -v 2
TEannot.py -P  R1A_annot -C TEannot.cfg -S 2 -a RM -r -v 2
TEannot.py -P  R1A_annot  -C TEannot.cfg -S 2 -a CEN -r -v 2
TEannot.py -P  R1A_annot -C TEannot.cfg -S 3 -c BLR+RM+CEN  -v 2 
TEannot.py -P  R1A_annot -C TEannot.cfg -S 7 -v 2

#Output:  R1A_annot_chr_allTEs_nr_join_path.annotStatsPerTE_FullLengthFrag.fa

#Second Round of TEannot using  the  R1A_annot_chr_allTEs_nr_join_path.annotStatsPerTE_FullLengthFrag.fa as FLF_refTE.fa

TEannot.py -P  R1A_FLF  -C TEannot.cfg -S 2 -a BLR  -v 2
TEannot.py -P  R1A_FLF  -C TEannot.cfg -S 2 -a RM -v 2
TEannot.py -P  R1A_FLF  -C TEannot.cfg -S 2 -a CEN -v 2
TEannot.py -P  R1A_FLF  -C TEannot.cfg -S 2 -a BLR -r -v 2
TEannot.py -P  R1A_FLF -C TEannot.cfg -S 2 -a RM -r -v 2
TEannot.py -P  R1A_FLF  -C TEannot.cfg -S 2 -a CEN -r -v 2
TEannot.py -P  R1A_FLF -C TEannot.cfg -S 3 -c BLR+RM+CEN  -v 2 
TEannot.py -P  R1A_FLF -C TEannot.cfg -S 7 -v 2
PostAnalyzeTELib.py -a 3 -g 42768942 -p   R1A_FLF_chr_allTEs_nr_noSSR_join_path -s R1A_FLF_refTEs_seq

#Concatenate .gff3 files
cat  R1A_FLF_GFF3/*.gff3 >>  R1A_TEs_final_annotation.gff3

#Softmasking using bedtools (V2.30.0)


#load modules and create conda environemn 
module load gcc12-env
module load miniconda3

conda create -n bedtools -c bioconda bedtools

conda activate bedtools
bedtools maskfasta -soft \
	-fi /home/R3-A.assembly.fasta \
	-bed TEs_final_annotation.gff3 \
	-fo /home/R3-A.assembly.softmasked_TE.fasta

conda deactivate
