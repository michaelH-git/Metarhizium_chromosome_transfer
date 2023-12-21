#!/bin/bash

module load gcc12-env
module load miniconda3
module load singularity
module load R/4.2.1 gcc/12.1.0


DIR=/your/Path/DNDS/


cd $DIR


#make all called SNPs in Mguizhouense homozygoutous

conda activate bcftools #bcftools version 1.17

VCF_IN="./phased.not_in_TEs.snps.inCDS.vcf"
VCF_OUT=".//phased.not_in_TEs.snps.inCDS.HOM.vcf"

cp $VCF_IN $VCF_OUT

sed -i -e 's/0|1/1\/1/g' -e 's/1|0/1\/1/g' -e 's/0\/1/1\/1/g' -e 's/;AF1=[^;]*;/;AF1=1;/g' -e 's/AC1=1/AC1=2/g' "$VCF_OUT"
bgzip -k $VCF_OUT && tabix -p vcf $VCF_OUT.gz

# generating a alternative genome file for Mguizhouense using the SNPs located in cds called on Mguizhouense reads mapped on Mrobertsii R3-I4


cat ./R3-I4.fasta | bcftools consensus ${VCF_OUT}.gz > ./Mguizhouense.consensus_from_R3-I4.fasta

conda deactivate


# Extracting cds for Mguizhouense.consensus_from_R3-I4 and Mrobertsii R3-I4 using AGAT (v1.0.0)


conda activate AGAT # AGAT version 1.0.0

agat_sp_extract_sequences.pl -g ./R3-I4.gff3 -f ./R3-I4.fasta -t cds -o R4-I4.cds.fasta -t cds -o ./R3-I4.cds.fasta
agat_sp_extract_sequences.pl -g ./R3-I4.gff3 -f ./Mguizhouense.consensus_from_R3-I4.fasta -t cds -o ./Mguizhouense.consensus_from_R3-I4.cds.fasta

conda deactivate


## alignement of the corresponing cds from Mrobertsii R3-I4 and Mguizhouense.consensus_from_R3-I4 using clustalO (version 1.2.4)


file1="./R3-I4.cds.fasta"
file2="./Mguizhouense.consensus_from_R3-I4.cds.fasta"

OUT_Dir=./CLUSTALO/

mkdir -p $OUT_Dir

conda activate seqtk

# Extract gene IDs from both files
ids_file1=$(grep "^>" "$file1" | sed 's/^>\([^[:space:]]*\).*/\1/')
ids_file2=$(grep "^>" "$file2" | sed 's/^>\([^[:space:]]*\).*/\1/')

# Find common gene IDs
common_ids=$(comm -12 <(echo "$ids_file1" | tr ' ' '\n' | sort) <(echo "$ids_file2" | tr ' ' '\n' | sort))

# Create temporary files for storing sequences
	list_fi=./list.txt
	fast_fi=${OUT_Dir}/tmp.fa
# Extract sequences for common gene IDs
for id in $common_ids; do
	echo $id > $list_fi
	seqtk subseq  $file1 $list_fi > $fast_fi
	seqtk subseq  $file2 $list_fi >> $fast_fi
	conda activate clustalo
	clustalo --outfmt clu -i $fast_fi -o ${OUT_Dir}/${id}-aligned_sequences.aln
	conda deactivate
done

conda deactivate



# Calculating dnds for each alignment using the ape package (version 5.7-1) in R 

export R_LIBS=$HOME/R_libs:$R_LIBS

Rscript - <<EOF

# Load the required library
library(ape)

# Set the directory path
directory_path <- "/CLUSTALO/"

# Get a list of all files in the directory
file_list <- list.files(path = directory_path, pattern = ".aln", full.names = TRUE)

# Create an empty data frame to store results
result_df <- data.frame(File = character(), dnds = numeric(), stringsAsFactors = FALSE)

# Loop through each file
for (file_path in file_list) {
  # Extract the file name without extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  file_name <- gsub("-aligned_sequences", "", file_name)
  dnds <- "NaN"
  
  tryCatch({
    # Read the sequence alignment using read.dna
    alignment <- read.dna(file_path, format = "clustal")
    
    # Calculate dN/dS
    dnds <- dnds(alignment, code = 1, quiet = TRUE, return.categories = FALSE)
    
    # Extract dN/dS ratio and disreagrd errors (if error occurs there are no SNPs and NaN is exported)
  }, error = function(e) {
    # If an error occurs, set dN/dS ratio to "NaN"
    dnds <- "NaN"
  })
  
  # Append results to the data frame
  result_df <- rbind(result_df, data.frame(File = file_name, dnds = dnds[1]))
}

# Write the results to a single file

write.table(result_df, file = "output_file.txt", sep = "\t", quote = FALSE, row.names = FALSE)

EOF



#generating a file that for each of the transcripts determines on which contig this transcript is located


# input files
gff=./R3-I4.gff3
dnds_file=./output_file.txt

# output files
output_file=./transcriptID_contigs.txt
merged_file=./dnds_transcriptID_contigs.txt
out_file=./dnds_final.txt


# extract "contig-ID" and "transcript-ID" for lines containing "transcript"
grep "transcript" "$gff" | awk '{split($9, a, ";"); sub(/ID=/, "", a[1]); print a[1], $1}' > "$output_file"

# Use join to merge based on the matching expression (first column)
join -1 1 -2 1 -o 1.1,1.2,2.2 <(sort -k1,1 "$output_file") <(sort -k1,1 "$dnds_file") > "$merged_file"

##filter out all genes with NaN or Inf or negative Value

awk -F' ' '!($3 == "NaN" || $3 == "Inf" || $3 < 0)' "$merged_file" >  "$out_file"

