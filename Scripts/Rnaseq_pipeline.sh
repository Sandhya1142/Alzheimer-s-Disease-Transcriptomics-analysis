#!/bin/bash
#SBATCH --job-name=rnaseq_AD
#SBATCH --partition=large
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --output=rnaseq_%j.log

module load anaconda/2024.10
source /iitgn/apps/packages/anaconda/2024.10/etc/profile.d/conda.sh
conda activate rnaseq_env

set -e 


#*************************************
#STEP 1 - Quality check using Fastqc
#*************************************
mkdir -p /iitgn/homedirs/sandhya/Rush_AD/fastqc_results
echo "Files are being checked for quality 🌻"
fastqc -t 40 /iitgn/homedirs/sandhya/Rush_AD/URLS_fastq/*.fastq.gz \
	-o /iitgn/homedirs/sandhya/Rush_AD/fastqc_results
echo "Quality check done!🌻"



#**************************************
#STEP 2 - Alignment
#**************************************
mkdir -p /iitgn/homedirs/sandhya/Rush_AD/aligned_Rush
echo "Alignment to reference genome has begun 🌻"

for file in /iitgn/homedirs/sandhya/Rush_AD/URLS_fastq/*.fastq.gz
do
	filename=$(basename ${file})
	sample=${filename%.fastq.gz}
	
	echo "processing $sample 🌻"
	
	STAR \
	--runThreadN 16 \
	--genomeDir /iitgn/homedirs/sandhya/reference/star_index_GRCh38 \
	--readFilesIn ${file} \
	--readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix /iitgn/homedirs/sandhya/Rush_AD/aligned_Rush/${sample}_

done
echo "Alignment done! 🌻"

#******************************************
#STEP 3 - Checking integrity of bam files
#******************************************
for sample in /iitgn/homedirs/sandhya/Rush_AD/aligned_Rush/*.bam
do
	
	base=$(basename ${sample})
	echo "processing $base"

	if samtools quickcheck "$sample" && samtools index "$sample"
	then
		echo "🌟 BAM passed validation 🌟"
	else
		echo "❌ ERROR !! BAM failed validation ❌"
	exit 1
	fi
done
echo "Files checked for integrity🌻"

#*******************************************
#STEP 4 - Quantifying the aligned reads
#*******************************************

mkdir -p /iitgn/homedirs/sandhya/Rush_AD/counts

echo "Quantification started!🌻"

featureCounts \
-T 40 \
-t exon \
-g gene_id \
-a /iitgn/homedirs/sandhya/annotation/gencode.v45.annotation.gtf \
-o /iitgn/homedirs/sandhya/Rush_AD/counts/gene_counts.txt \
/iitgn/homedirs/sandhya/Rush_AD/aligned_Rush/*_Aligned.sortedByCoord.out.bam

echo "Quantification done!🌻"
















