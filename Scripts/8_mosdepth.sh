#!/bin/bash
#SBATCH --job-name=mosdepth
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time 12:00:00
#SBATCH -o _%A_%a.out
#SBATCH --error=_%A_%a.err
#SBATCH --array=0-1
#SBATCH --mem=15G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=



module load gcc12-env
module load miniconda3


#activate create conda environment

conda create -n mosdepth -c bioconda mosdepth # mosdepth (v 0.3.4)


input_dir=/your/path/04_reformat/
output_dir=/your/path/mosdepth/

bed_file=/your/path/R3-I4_TE-complement_50kbwindows.sorted.bed

mkdir -p $output_dir




for bam in $input_dir/*.bam; do
	conda activate mosdepth
	outbase=${output_dir}"$(basename "$bam" on-R3I4.RG.Dedup.bam)"	
	mosdepth -t 10 -n -Q 10 -T 5 -b $bed_file $outbase $bam 
	conda deactivate
done

