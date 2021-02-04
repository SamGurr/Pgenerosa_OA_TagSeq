#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=samuel_gurr@uri.edu
#SBATCH --output=../../../sgurr/Geoduck_TagSeq/output/clean/"%x_out.%j"
#SBATCH --error=../../../sgurr/Geoduck_TagSeq/output/clean/"%x_err.%j"
#SBATCH -D /data/putnamlab/KITT/hputnam/20201217_Geoduck_TagSeq/

# load modules needed
module load fastp/0.19.7-foss-2018b
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.7-foss-2018b-Python-2.7.15

# Make an array of sequences to trim
array1=($(ls *.fastq.gz)) 

# fastp loop; trim the Read 1 TruSeq adapter sequence 
for i in ${array1[@]}; do
	fastp --in1 ${i} --out1 ../../../sgurr/Geoduck_TagSeq/output/clean/clean.${i} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
        fastqc ../../../sgurr/Geoduck_TagSeq/output/clean/clean.${i}
done 

echo "Read trimming of adapters complete." $(date)

# Quality Assessment of Trimmed Reads

cd ../../../sgurr/Geoduck_TagSeq/output/clean #The following command will be run in the /clean directory

multiqc ./ #Compile MultiQC report from FastQC files

echo "Cleaned MultiQC report generated." $(date)
