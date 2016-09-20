#this script works to prepare and conduct the initial mapping of bisulfite sequencing reads for the diverse accessions

#bismark and bowtie1 will need to be installed and in your PATH
#samtools, fastqc, trim_galore, cutadapt as well

#index and prepare genomes
bismark_genome_preparation rawdata/genomes/Bd21/ &
bismark_genome_preparation rawdata/genomes/Bd21-3/ &
bismark_genome_preparation rawdata/genomes/Bd3-1/ &
bismark_genome_preparation rawdata/genomes/Bd30-1/ &
bismark_genome_preparation rawdata/genomes/Bd1-1/ &
bismark_genome_preparation rawdata/genomes/BdTR12c/ &
bismark_genome_preparation rawdata/genomes/Koz-3/ &

#perform alignments
#note that this takes some time (days). It could likely be done better given the power of your machine and running things in parallel

mkdir 0_initial_mapping
cd r0_initial_mapping

cp ../rawdata/reference_methylomes/Bd21_R1.fastq.gz .
cp ../rawdata/reference_methylomes/Bd21_R2.fastq.gz .
../wgbs_alignment_script.sh -pese Bd21_R1.fastq.gz Bd21_R2.fastq.gz ../rawdata/genomes/Bd21/ Bd21


cp ../rawdata/reference_methylomes/Bd21-3_R1.fastq.gz .
cp ../rawdata/reference_methylomes/Bd21-3_R2.fastq.gz .
../wgbs_alignment_script.sh -pese Bd21-3_R1.fastq.gz Bd21-3_R2.fastq.gz ../rawdata/genomes/Bd21-3/ Bd21-3


cp ../rawdata/reference_methylomes/Bd3-1_R1.fastq.gz .
cp ../rawdata/reference_methylomes/Bd3-1_R2.fastq.gz .
../wgbs_alignment_script.sh -pese Bd3-1_R1.fastq.gz Bd3-1_R2.fastq.gz ../rawdata/genomes/Bd3-1/ Bd3-1


cp ../rawdata/reference_methylomes/Bd30-1_R1.fastq.gz .
cp ../rawdata/reference_methylomes/Bd30-1_R2.fastq.gz .
../wgbs_alignment_script.sh -pese Bd30-1_R1.fastq.gz Bd30-1_R2.fastq.gz ../rawdata/genomes/Bd30-1/ Bd30-1


cp ../rawdata/reference_methylomes/Bd1-1_R1.fastq.gz .
cp ../rawdata/reference_methylomes/Bd1-1_R2.fastq.gz .
../wgbs_alignment_script.sh -pese Bd1-1_R1.fastq.gz Bd1-1_R2.fastq.gz ../rawdata/genomes/Bd1-1/ Bd1-1


cp ../rawdata/reference_methylomes/BdTR12c_R1.fastq.gz .
cp ../rawdata/reference_methylomes/BdTR12c_R2.fastq.gz .
../wgbs_alignment_script.sh -pese BdTR12c_R1.fastq.gz BdTR12c_R2.fastq.gz ../rawdata/genomes/BdTR12c/ BdTR12c


cp ../rawdata/reference_methylomes/Koz-3_R1.fastq.gz .
cp ../rawdata/reference_methylomes/Koz-3_R2.fastq.gz .
../wgbs_alignment_script.sh -pese Koz-3_R1.fastq.gz Koz-3_R2.fastq.gz ../rawdata/genomes/Koz-3/ Koz-3

#pull out all the files that are useful

#cleanup (which you can do earlier)

