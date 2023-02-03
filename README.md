# RNA-Seq-kallisto
Here we describe a pipeline to analyse RNA sequencing data (NovaSeq 6000 PE150). Overall, we need to 1) check the quality of raw data (FASTQC analyisis), 2) filter te raw data (FASTP), and 3) estimate the gene expression of the RNAseq data (KALLISTO). Kallisto is a program for quantifying abundances of transcripts from bulk RNA-Seq data. It is based on the pseudoalignment of reads against the transcriptome of a species for rapid estimation of transcript expression. Of note, this pipeline is described considering that the genome/transcriptome of the species that we are working with is publicly available.

## 0. Download files
The first step is downloading the raw sequencing files into our server of computer. Most sequencing providers will upload the raw files into an FTP server. The files can usually be downloaded using the Unix command "wget":
```
wget https://...
```
Paired-end sequencing generates two output fastq files for each sample. These usually end in "_1.fastq.gz" and "_2.fastq.gz" (but "_R1_001.fastq.gz" and "_R2_001.fastq.gz" are also frequent). The files are usually compressed using gunzip - there is no need to unzip them, most software will work with gunzipped files. The fastq files will usually come with corresponding MD% files to check that the integrity of the downloaded files.


## 1. Quality control with FastQC 
Once we have downloaded the fastq files, it is good practice to check their quality. It can be done using the software [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). The command to run FastQC in a sample (a pair of fastq files) is:

```
fastqc sampleA_1.fastq.gz
```
Usually, it will be more convinient to run FastQC for all samples at once:
```
fastqc *fastq.gz
```
#You will obtain 2 output files from each sample: fastqc.zip and fastqc.html, which you might organize in a new folder called 'QC_samples'

DUDA: Cómo harías para meter los archivos de salida directamente en una nueva carpeta, en vez de crearla después y mover los archivos a ella? Probé a crear la carpeta qc_samples dentro de raw_data y nombrarla al final del comando, pero no funcionó

#4 DOWNLOAD THE TRANSCRIPTOME OF THE SPECIES YOU ARE WORKING WITH
#You can create a new folder called 'transcriptome' 
```
mkdir transcriptome
```
#How to obtain the transcriptome link to download? NCBI -> Genome -> Transcript -> Right mouse botton -> copy link address
```
wget LINK
```
#5 CREATE KALLISTO INDEX. You will need to create a Kallisto index from the transcriptome, which means that you are giving the transcriptome the instrunctions of how to use the Kallisto to match the transcriptome with the RNAseq data
#You should do this from the folder 'transcriptome'
```
module load cesga/2018  gcccore/6.4.0 #This will be only if you are using 'Cesga-Galicia'#You might need other programs depending on your system and you can check it with: module spider kallisto
module load kallisto/0.46.1
kallisto index -i transcriptome_rabbit_index_kallisto.idx rabbit_transcriptome.fna.gz
```
#6 FILTERING RAW DATA - FASP
#We will use FASP to filter the data. This will remove low quality reads, contaminating sequence, low complexity reads (repeats), short reads, etc. We will do this with the 'raw_data', NOT with the QC data.
#We can do this step from a general folder called 'RABBIT' for example:
```
module load cesga/2018 gcccore/6.4.0 #This will be only if you are using 'Cesga-Galicia'
module load fastp/0.20.1
```
#You can create a new folder called inside 'RABBIT' called 'filtered'. In 'RABBIT' you will also have the folder 'raw_data' which contain all the RNAseq raw data files (fastq.gz)
```
mkdir filtered
```
#You have to indicate both files (forward and reverse)/each sample (--in 1 and --in 2), and also the folder (create it in advance) and names of the output files that are still paired (--out1, --out2). --unpaired1, --unpaired2 indicate the output files for the reads that are not 'paired'. -l 35: this specifies that if a read is shorter than 50 basepairs after all filters, it should be removed. -q: quality threshold per base required (default: 15, which means that a Phred quality score of at least 15 is required) -h: specifies name for the html file with plots showing the read quality before and after filtering 
```
fastp --in1 raw_data/E283_1.fq.gz --in2 raw_data/E283_2.fq.gz --out1 filtered/E283_1.filtered.fastq.gz --out2 filtered/E283_2.filtered.fastq.gz --unpaired1 filtered/E283_1.unpaired.fastq.gz --unpaired2 filtered/E283_2.unpaired.fastq.gz -l 35 -q 20 -h filtered/E283.html &> filtered/E283.log
```
#You will have to do this for each of your samples. If you have many samples, you can run a '.sh file' for all samples
#How to create a .sh file? #Write the following in a text document:

```
#!/bin/sh
#SBATCH -p thinnodes 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 15:00:00
#SBATCH --mail-type=begin 
#SBATCH --mail-type=end 
#SBATCH --mail-user=paularodriguez.villamayor@usc.es #Direccion

module load cesga/2018 gcccore/6.4.0
module load fastp/0.20.1

fastp --in1 raw_data/P03_1.fq.gz --in2 raw_data/P03_2.fq.gz --out1 filtered/P03_1.filtered.fastq.gz --out2 filtered/P03_2.filtered.fastq.gz --unpaired1 filtered/P03_1.unpaired.fastq.gz --unpaired2 filtered/P03_2.unpaired.fastq.gz -l 35 -q 20 -h filtered/P03.html &> filtered/P03.log
fastp --in1 raw_data/P04_1.fq.gz --in2 raw_data/P04_2.fq.gz --out1 filtered/P04_1.filtered.fastq.gz --out2 filtered/P04_2.filtered.fastq.gz --unpaired1 filtered/P04_1.unpaired.fastq.gz --unpaired2 filtered/P04_2.unpaired.fastq.gz -l 35 -q 20 -h filtered/P04.html &> filtered/P04.log
fastp --in1 raw_data/P05_1.fq.gz --in2 raw_data/P05_2.fq.gz --out1 filtered/P05_1.filtered.fastq.gz --out2 filtered/P05_2.filtered.fastq.gz --unpaired1 filtered/P05_1.unpaired.fastq.gz --unpaired2 filtered/P05_2.unpaired.fastq.gz -l 35 -q 20 -h filtered/P05.html &> filtered/P05.log
fastp --in1 raw_data/P06_1.fq.gz --in2 raw_data/P06_2.fq.gz --out1 filtered/P06_1.filtered.fastq.gz --out2 filtered/P06_2.filtered.fastq.gz --unpaired1 filtered/P06_1.unpaired.fastq.gz --unpaired2 filtered/P06_2.unpaired.fastq.gz -l 35 -q 20 -h filtered/P06.html &> filtered/P06.log
fastp --in1 raw_data/P41_1.fq.gz --in2 raw_data/P41_2.fq.gz --out1 filtered/P41_1.filtered.fastq.gz --out2 filtered/P41_2.filtered.fastq.gz --unpaired1 filtered/P41_1.unpaired.fastq.gz --unpaired2 filtered/P41_2.unpaired.fastq.gz -l 35 -q 20 -h filtered/P41.html &> filtered/P41.log
fastp --in1 raw_data/P42_1.fq.gz --in2 raw_data/P42_2.fq.gz --out1 filtered/P42_1.filtered.fastq.gz --out2 filtered/P42_2.filtered.fastq.gz --unpaired1 filtered/P42_1.unpaired.fastq.gz --unpaired2 filtered/P42_2.unpaired.fastq.gz -l 35 -q 20 -h filtered/P42.html &> filtered/P42.log
fastp --in1 raw_data/P43_1.fq.gz --in2 raw_data/P43_2.fq.gz --out1 filtered/P43_1.filtered.fastq.gz --out2 filtered/P43_2.filtered.fastq.gz --unpaired1 filtered/P43_1.unpaired.fastq.gz --unpaired2 filtered/P43_2.unpaired.fastq.gz -l 35 -q 20 -h filtered/P43.html &> filtered/P43.log
```

#You convert the .txt document into .sh 
#I created the .txt in Windows, then converted into .sh and then transfered to Linux #You need to use the command dos2unix to transfer from windows to linux
```
dos2unix name.sh
```
#Run the file
```
sbatch name.sh #In this case it is important to put the .sh file and to run it from the 'RABBIT' folder but it will depend on how you organize your folders
```
#Results: you will obtain different files for each sample but the 2 important ones to run Kallisto are: _1.filtered.fastq.gz, and _2.filtered.fastq.gz


#ESTIMATE THE EXPRESSION - KALLISTO
#You will now match the transcriptome with the FILTERED data to obtain the gene expression of your RNAseq. Do this from the folder 'RABBIT'. Kallisto will create the '-o folder' inside 'RABBIT', so don't create it in advance, and you will 'call' your transcriptome with '-i'. The 'transcriptome' folder should be inside 'RABBIT'. -b 100 indicates the maximum lenght of a read (not sure about this?¿). Finally you have to indicate the directory of the filtered data (both paired reads)
#You can do that with a compute session if you have few samples. Else, you can queue them with .sh and 'sbatch' (same procedure as for FASTP)
#How to create a .sh file to run kallisto? You have to do a different script for each sample and run them independently. See below an example for ONE SAMPLE.
```
#!/bin/sh
#SBATCH -p thinnodes
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 2:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=paularodriguez.villamayor@usc.es

module load cesga/2018 gcccore/6.4.0 kallisto/0.46.1

kallisto quant -i transcriptome/transcriptome_rabbit_index_kallisto.idx -o kallisto_results_E283 -b 100 /mnt/lustre/scratch/home/usc/ge/prv/fetos/filtered/E283_1.filtered.fastq.gz /mnt/lustre/scratch/home/usc/ge/prv/fetos/filtered/E283_2.filtered.fastq.gz
```

#You convert the .txt document into .sh 
#I created the .txt in Windows, then converted into .sh and then transfered to Linux #You need to use the command dos2unix to transfer from windows to linux
```
dos2unix name.sh
```
#Run the file
```
sbatch name.sh #it is important to do this from the folder 'RABBIT' but it will depend on how you organize your folders
```
#Results: Each sample will have its own folder (i.e. kallisto_results_E283) and inside that folder you will find 3 files: abundance.h5, abundance.tsv, and run_info.json





