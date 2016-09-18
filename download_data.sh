mkdir rawdata
cd rawdata

#get the seven reference methylomes (single replicates) from NCBI SRA
mkdir reference_methylomes
cd reference_methylomes
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR197/SRR1972494/SRR1972494.sra #Bd21
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR197/SRR1972495/SRR1972495.sra #Bd21-3
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR197/SRR1972496/SRR1972496.sra #Bd3-1
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR197/SRR1972497/SRR1972497.sra #Bd30-1
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR197/SRR1972498/SRR1972498.sra #Bd1-1
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR197/SRR1972499/SRR1972499.sra #BdTR12c
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR197/SRR1972500/SRR1972500.sra #Koz-3

for FILE in *.sra
do
fastq-dump --gzip --split-files --readids --skip-technical --dumpbase --clip $FILE
rm $FILE
done

mv SRR1972494_1.fastq.gz Bd21_R1.fastq.gz
mv SRR1972494_2.fastq.gz Bd21_R2.fastq.gz

mv SRR1972495_1.fastq.gz Bd21-3_R1.fastq.gz
mv SRR1972495_2.fastq.gz Bd21-3_R2.fastq.gz

mv SRR1972496_1.fastq.gz Bd3-1_R1.fastq.gz
mv SRR1972496_2.fastq.gz Bd3-1_R2.fastq.gz

mv SRR1972497_1.fastq.gz Bd30-1_R1.fastq.gz
mv SRR1972497_2.fastq.gz Bd30-1_R2.fastq.gz

mv SRR1972498_1.fastq.gz Bd1-1_R1.fastq.gz
mv SRR1972498_2.fastq.gz Bd1-1_R2.fastq.gz

mv SRR1972499_1.fastq.gz BdTR12c_R1.fastq.gz
mv SRR1972499_2.fastq.gz BdTR12c_R2.fastq.gz

mv SRR1972500_1.fastq.gz Koz-3_R1.fastq.gz
mv SRR1972500_2.fastq.gz Koz-3_R2.fastq.gz


cd ../
#get 14 replicate samples of Bd1-1, Bd21, and Bd3-1 from NCBI SRA
mkdir replicate_methylomes
cd replicate_methylomes
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR247/SRR2475390/SRR2475390.sra #1-1 rep1
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR247/SRR2475391/SRR2475391.sra #1-1 rep2
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR247/SRR2475392/SRR2475392.sra #1-1 rep3
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR247/SRR2475393/SRR2475393.sra #1-1 rep4
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR247/SRR2475394/SRR2475394.sra #1-1 rep5
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR247/SRR2475395/SRR2475395.sra #21 rep1
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR247/SRR2475396/SRR2475396.sra #21 rep2
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR247/SRR2475397/SRR2475397.sra #21 rep3
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR247/SRR2475398/SRR2475398.sra #21 rep4
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR247/SRR2475399/SRR2475399.sra #21 rep5
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR247/SRR2475400/SRR2475400.sra #3-1 rep2
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR247/SRR2475401/SRR2475401.sra #3-1 rep3
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR247/SRR2475402/SRR2475402.sra #3-1 rep4
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR247/SRR2475403/SRR2475403.sra #3-1 rep5

for FILE in *.sra
do
fastq-dump --gzip --split-files --readids --skip-technical --dumpbase --clip $FILE
rm $FILE
done

mv SRR2475390_1.fastq.gz Bd1-1_rep1.fastq.gz
mv SRR2475391_1.fastq.gz Bd1-1_rep2.fastq.gz
mv SRR2475392_1.fastq.gz Bd1-1_rep3.fastq.gz
mv SRR2475393_1.fastq.gz Bd1-1_rep4.fastq.gz
mv SRR2475394_1.fastq.gz Bd1-1_rep5.fastq.gz
mv SRR2475395_1.fastq.gz Bd21_rep1.fastq.gz
mv SRR2475396_1.fastq.gz Bd21_rep2.fastq.gz
mv SRR2475397_1.fastq.gz Bd21_rep3.fastq.gz
mv SRR2475398_1.fastq.gz Bd21_rep4.fastq.gz
mv SRR2475399_1.fastq.gz Bd21_rep5.fastq.gz
#note Bd3-1_rep1 does not exist (failed sample, no library, not sequenced, no SRA run)
mv SRR2475400_1.fastq.gz Bd3-1_rep2.fastq.gz
mv SRR2475401_1.fastq.gz Bd3-1_rep3.fastq.gz
mv SRR2475402_1.fastq.gz Bd3-1_rep4.fastq.gz
mv SRR2475403_1.fastq.gz Bd3-1_rep5.fastq.gz





cd ../


#download the SNP data from Gordon et al., 2014 Brachypodium genomic resequencing paper
mkdir reference_snp_data
cd reference_snp_data

#TBD as brachypodium.org is down now....

wget https://www.dropbox.com/s/dol2ykdpyunw5y1/Bdistachyon_283.vcf.gz

cd ../../

#confirm data files are as expected
shasum -c rawdata_shasums.sha




