# Geoduck TagSeq Bioinformatics

## <span style="color:blue">**Table of Contents**</span>
  - [Upon upload to HPC...](#Initial-diagnostics-upon-sequence-upload-to-HPC)
      - [Count raw reads](#Count-the-number-of-read-files)
      - [Digital fingerprint md5sums](#run-checksum) (HPC script)
  - [Initial quality check](#Quality-check-of-raw-reads) (HPC script)
  - [Trim and post-trim quality check](#Trimming-and-post-trim-quality-check-of-'clean'-reads) (HPC script)
  - [HISAT2: Alignment of cleaned reads to reference](#Alignment-of-cleaned-reads-to-reference)
	- [Upload reference genome](#Reference-genome-upload-to-HPC)
	- [About HISAT2](#HISAT2-alignment)
	- [About samtools](#samtools)
	- [Run index and alignment](#HPC-Job-HISAT2-Index-Reference-and-Alignment) (HPC script)
  - [StringTie: Assembly](#Assembly-and-quantification)  
	- [About StringTie](#StringTie)
## Initial diagnostics upon sequence upload to HPC
--------------------------------------------
### Count the number of read files
ls -1 | wc -l

should equal 141 seq samples *2lanes per sample + 1 md5 = 283

### run checksum

```
nano transfer_checks.sh
```
```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH -D /data/putnamlab/KITT/hputnam/20201217_Geoduck_TagSeq/
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

# load modules needed

# generate md5
md5sum /data/putnamlab/KITT/hputnam/20201217_Geoduck_TagSeq/*.gz > /data/putnamlab/KITT/hputnam/20201217_Geoduck_TagSeq/URI.md5

# Count the number of sequences per file

zcat /data/putnamlab/KITT/hputnam/20201217_Geoduck_TagSeq/*fastq.gz | echo $((`wc -l`/4)) > /data/putnamlab/KITT/hputnam/20201217_Geoduck_TagSeq/rawread.counts.txt
```

```
sbatch transfer_checks.sh
```

- check the digital fingerprint of the files with md5sum 
- compare md5sum of our output URI.md5 file to the UT Library Information pdf; okay the upload - move forward


## Quality check of raw reads
-------------------------------------------
### RUN FASTQC and MUTLIQC


### shell script: <span style="color:green">**fastqc_raw.sh**<span>
```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH -D /data/putnamlab/KITT/hputnam/20201217_Geoduck_TagSeq/
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

# load modules needed
module load all/FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

for file in /data/putnamlab/KITT/hputnam/20201217_Geoduck_TagSeq/*gz
do
fastqc $file
done

multiqc /data/putnamlab/KITT/hputnam/20201217_Geoduck_TagSeq/

```

```
sbatch fastqc_raw.sh
```
- view the multiqc.html report 

## Trimming and post-trim quality check of 'clean' reads
-------------------------------------------
- notes from the previous raw data multiqc.html report:
	- *Per Sequence GC Content*: double peak of ~10% GC and 30-40% GC
	- *Mean Quality Scores*: > 30
	- *Per Base Sequence Content*: 15-25% C and G, 25-28% T, 35-40+% A
	- *Adapter Content*: high adapters present, 1-6% of sequences 
- trim adapter sequence
- output mutliqc report of the 'clean' reads
### shell script: <span style="color:green">**fastp.sh**<span>
```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=samuel_gurr@uri.edu
#SBATCH --output=../../../sgurr/Geoduck_TagSeq/output/clean/"%x_out.%j"
#SBATCH --error=../../../sgurr/Geoduck_TagSeq/output/clean/"%x_err.%j"
#SBATCH -D /data/putnamlab/KITT/hputnam/20201217_Geoduck_TagSeq/
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

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
```

### EXPORT MUTLIQC REPORT 
*exit bluewaves and run from terminal*
- save to gitrepo as multiqc_clean.html
```
scp samuel_gurr@bluewaves.uri.edu:/data/putnamlab/sgurr/Geoduck_TagSeq/output/clean/multiqc_report.html  C:/Users/samjg/Documents/My_Projects/Geoduck_TagSeq
```

## Alignment of cleaned reads to reference
-------------------------------------------

### Reference genome upload to HPC 

*exit bluewaves and run from terminal*

```
scp C:/Users/samjg/Documents/Bioinformatics/genomes/Panopea-generosa-v1.0.fa samuel_gurr@bluewaves.uri.edu:/data/putnamlab/sgurr/refs/
```

-  file name: Panopea-generosa-genes.fna
	- Roberts, S., White, S., Trigg, S. A., and Putnam, H.M. (2020). Panopea_generosa_genome. OSF. June 24. doi:10.17605/OSF.IO/YEM8N.
	 [Roberts et al. 2020; link to genome](doi:10.17605/OSF.IO/YEM8N)
-  reference genome file size: 368MB
-  uploaded to sgurr/refs/

### HISAT2 alignment

**About:** HISAT2 is a sensitive alignment program for mapping sequencing reads.
In our case here, we will use HISAT2 to **(1)** index our *P. generosa* reference genome **(2)** align our clean TagSeq reads to the indexed reference genome. 

More information on HISAT2 can be read [here](http://daehwankimlab.github.io/hisat2/manual/)!

**Main arguments used below**: 

**(1)** *Index the reference genome*

``` hisat2-build ``` =
builds a HISAT2 index from a set of DNA sequences. Outputs 6 files that together consitute the index. 
ALL output files are needed to slign reads to the reference genome and the original sequence FASTA file(s) 
are no longer used for th HISAT2 alignment 

``` -f <reads.fasta> ``` =
the reads (i.e <m1>, <m2>, <m100>)
FASTA files usually have extension .fa, .fasta, .mfa, .fna or similar. 
FASTA files do not have a way of specifying quality values, so when -f is set, 
the result is as if --ignore-quals is also set

Note: other options for your reads are ```-q ``` for FASTQ files, 
 ```--qseq ``` for QSEQ  files, etc. - check [here](http://daehwankimlab.github.io/hisat2/manual/) for more file types
 
**(2)** *align reads to reference*

``` -x <hisat2-indx> ``` =
base name for the index reference genome, first looks in the current directory then in the directory specified in HISAT_INDEXES environment variable

``` --dta ``` =
 **important!** reports the alignments tailored for *StringTie* assembler.
 With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. 
 This leads to fewer alignments with short-anchors, which helps transcript assemblers improve 
 significantly in computation and memory usage. 
 
 This is important relative to Cufflinks ```--dta-cufflinks``` 
 in which HISAT2 looks for novel splice sites with three signals (GT/AG, GC/AG, AT/AC)

``` -U <r> ``` =
Comma-separated list of files contained unparied reads to be aligned. 
*In other words...*, our array of post-trimmed 'clean' TagSeq reads

``` -p NTHREADS``` =
Runs on separate processors, increasing -p increases HISAT2's memory footprint, increasing -p from 1 to 8 
increased the footprint by a few hundred megabytes 

Note: Erin Chile from Putnam Lab ran -p 8 for HISAT2 alignment of *Montipora* [here](https://github.com/echille/Montipora_OA_Development_Timeseries/blob/master/Amil/amil_RNAseq-analysis.sh)

``` -S <hit> ``` =
file to write SAM alignments to. By default, alignments are written to the “standard out” or “stdout” filehandle (i.e. the console).

**Note:** HISAT2 also has several criteria to trim such as the phred score and base count 5' (left) and 3' (right) 
Since we already assessed quality and trimmed, we will not use these commands


### samtools 
**About:** used to manipuate alignments to a binary BAM format - files contain the spliced reads alignemnts sorted by the reference position
with a tag to indicate the genomic strand that produced the RNA from which the read was sequenced. 
samtools quickly extracts alignments overlapping particular genomic regions - outputs allow viewers to quickly display alignments in each genomic region
note: the SAM format output by HISAT2 must be sorted and converted to BAM format using the samtools program

``` sort ``` =
sorts the alignments by the leftmost coordinates

```-o <out.bam``` = 
outputs as a specified .bam file 
	 
```-@ threads``` = 
Set number of sorting and compression threads. By default, operation is single-threaded
	 
more information on samtools commands [here](http://www.htslib.org/doc/1.1/samtools.html)
	 
## HPC Job: HISAT2 Index Reference and Alignment 

- create directory output\hisat2  

``` mkdir hisat2 ```

- index reference and alignment 

**input**
- Panopea-generosa-genes.fna *= reference genome*
- clean/*.fastq.gz *= all clean TagSeq reads*

**ouput**
- Pgenerosa_ref *= indexed reference by hisat2-build; stored in the output/hisat2 folder as 1.hy2, 2.ht2... 8.ht2*
- <clean.fasta>.sam *=hisat2 output, readable text file; removed at the end of the script*
- <clean.fasta>.bam *=converted binary file complementary to the hisat sam files*

### shell script: <span style="color:green">**hisat2.sh**<span>

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=samuel_gurr@uri.edu
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"
#SBATCH -D /data/putnamlab/sgurr/Geoduck_TagSeq/output/hisat2

#load packages
module load HISAT2/2.1.0-foss-2018b #Alignment to reference genome: HISAT2
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

# symbolically link 'clean' reads to hisat2 dir
ln -s ../clean/*.fastq.gz ./

# index the reference genome for Panopea generosa output index to working directory
hisat2-build -f ../../../refs/Panopea-generosa-v1.0.fa ./Pgenerosa_ref
echo "Referece genome indexed. Starting alingment" $(date)

# This script exports alignments as bam files
# sorts the bam file because Stringtie takes a sorted file for input (--dta)
# removes the sam file because it is no longer needed
array=($(ls *.fastq.gz)) # call the symbolically linked sequences - make an array to align
for i in ${array[@]}; do
        hisat2 -p 8 --dta -x Pgenerosa_ref -U ${i} -S ${i}.sam
        samtools sort -@ 8 -o ${i}.bam ${i}.sam
                echo "${i} bam-ified!"
        rm ${i}.sam
done
```
- HISAT2 complete with format prepared for StringTie assembler!

## Assembly and quantification'

### Upload annotation reference .gff or .gff3 to HPC 

- Sam White and Steven Roberts completed the annotation with open resources available for download

*exit bluewaves and run from terminal*

- "one file to rule them all" (- Sam W.) with the CDS, mRNA, exon, and genes
```
scp C:/Users/samjg/Documents/Bioinformatics/genomes/Pgenerosa_annotations/Panopea-generosa-vv0.74.a3-merged-2019-09-03-6-14-33.gff3 samuel_gurr@bluewaves.uri.edu:/data/putnamlab/sgurr/refs/
```

- genes only
```
scp C:/Users/samjg/Documents/Bioinformatics/genomes/Pgenerosa_annotations/Panopea-generosa-v1.0.a4.gene.gff3 samuel_gurr@bluewaves.uri.edu:/data/putnamlab/sgurr/refs/
```

- CDS only
```
scp C:/Users/samjg/Documents/Bioinformatics/genomes/Pgenerosa_annotations/Panopea-generosa-v1.0.a4.CDS.gff3 samuel_gurr@bluewaves.uri.edu:/data/putnamlab/sgurr/refs/
```

### StringTie

**About:** StringTie is a fast and efficient assembler of RNA-Seq alignments to assemble and quantitate full-length transcripts.  putative transcripts. 
For our use, we will input short mapped reads from HISAT2. StringTie's ouput can be used to identify DEGs in programs such as DESeq2 and edgeR

More information on StringTie can be read [here](https://ccb.jhu.edu/software/stringtie/)!

[Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/) is another assembler option we will not use. [Read this](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4643835/#:~:text=On%20a%20simulated%20data%20set,other%20assembly%20software%2C%20including%20Cufflinks.)

**Main StringTie arguments used below**: 

``` -p ``` =
specify the number of threads (CPUs). Note that a single node has at least 8 CPUs or 'threads' to call here

``` -A ``` =
Gene abundances will be reported (tab delimited format) in the output file with the given name.

``` -e ``` = 
Limits the processing of read alignments to only estimate and output the assembled transcripts matching the 
reference transcripts given with the -G option (*requires -G, recommended for -B/-b*). **With this option, 
read bundles with no reference transcripts will be entirely skipped**, which may provide a considerable 
speed boost when the given set of reference transcripts is limited to a set of target genes, for example.

``` -G <ref_ann.gff>``` =
Use the reference annotation file (in GTF or GFF3 format) to guide the assembly process. 
The output will include expressed reference transcripts as well as any novel transcripts that are assembled. 
*This option is required by options -B, -b, -e, -C*.

``` -o ``` =
Sets the name of the output GTF file where StringTie will write the assembled transcripts

``` --merge ``` =
StringTie takes as input a list of GTF/GFF files and merges/assembles these transcripts into a non-redundant set of transcripts. This mode is used in the new differential analysis pipeline to generate a global, unified set of transcripts (isoforms) across multiple RNA-Seq samples.
If the -G option (reference annotation) is provided, StringTie will assemble the transfrags from the input GTF files with the reference transcripts.

- merge commands

 ```-G <guide_gff>```	reference annotation to include in the merging (GTF/GFF3)
 
 ```-o <out_gtf>```	output file name for the merged transcripts GTF (default: stdout)

**gtf_list.txt**  ``` --merge``` requires a list file to call each of the gtf files 

*from within the stringtie directory with all .gtf outputs of the for loop, run the following to build txt file* **gtf_list.txt**

**``` ls *.gtf > gtf_list.txt```**

**Main gffcompare arguments used below**:

#### gffcompare
**About:** gffcompare is used to estimate the accuracy of one or more GFF query files compared to a reference annotation
``` -r ``` = 
reference annotation GFF file - each file is matched against it and tagged as overlapping, matching, or novel where appropriate.
outputs .refmap and .tmap files 

``` -G ``` =

``` -o <outprefix>``` = give a prefix for output files created by gffcompare (i.e. merged)

#### Python step prepDE.py**

**About:** in preparation for differential expression analysis using Bioconductor packages in R (i.e. DESeq2, edgeR, WGCNA), 
we will need to first aquire a matrix of read counts to particular genomic features (i.e. genes). prepDE.py is a Python call to extract 
hypothetical read counts for each transcript from the files generated by ```-e``` using StringTie (here as ${i}.gtf from .bam files by HISAT2)

```prepDE.py``` =
.py found [here](https://github.com/gpertea/stringtie/blob/master/prepDE.py) - add to scripts folder to call using python 
prepDE.py creats the matrix of read counts 

``` -g ``` =
where to ouput the gene count matrix (defaults as gene_count_matrix.csv)
```-i INPUT``` = 
a folder containing all sample sub-directories; Alternatively can provide a text file (i.e. sample_list.tct) with sample ID and path to its 
GTF file on each line (default '.' meaning the working directory is the subdirectory with all GTF files) 
*Alternatively call a **listGTF.txt** file...* this file has two columns with the sample ID and the <Path to sample> for the .gtf files 
run the following **```for filename in *.gtf; do echo $filename $PWD/$filename; done > listGTF.txt```** to make this text file and call in ``` -i ``` for this python step

## HPC job: Assembly 

### shell script: <span style="color:green">**stringtie.sh**<span>
```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=200GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=samuel_gurr@uri.edu
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"
#SBATCH -D /data/putnamlab/sgurr/Geoduck_TagSeq/output/stringtie

#load packages
module load StringTie/2.1.1-GCCcore-7.3.0 #Transcript assembly: StringTie

array=($(ls ../hisat2/*.bam)) #Make an array of sequences to assemble

for i in ${array[@]}; do #Running with the -e option to compare output to exclude novel genes. Also output a file with the gene abundances
        stringtie -p 8 -e -G ../../../refs/Panopea-generosa-vv0.74.a3-merged-2019-09-03-6-14-33.gff3 -A ../stringtie/${i}.gene_abund.tab -o ../stringtie/${i}.gtf
        echo "StringTie assembly for seq file ${i}" $(date)
done
echo "StringTie assembly COMPLETE, starting assembly analysis" $(date)
```

## HPC job: Merge and Build Read Count Matrix for DEG analysis

### shell script: <span style="color:green">**merge_prepDEpy.sh**<span>
```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=200GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=samuel_gurr@uri.edu
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"
#SBATCH -D /data/putnamlab/sgurr/Geoduck_TagSeq/output/stringtie

#load packages
module load Python/2.7.15-foss-2018b #Python
module load StringTie/2.1.1-GCCcore-7.3.0 #Transcript assembly: StringTie
module load gffcompare/0.11.5-foss-2018b #Transcript assembly QC: GFFCompare

stringtie --merge -p 8 -G ../../../refs/Panopea-generosa-vv0.74.a3-merged-2019-09-03-6-14-33.gff3 -o Pgenerosa_merged.gtf gtf_list.txt #Merge GTFs to form full GTF for analysis of assembly accuracy and precision
echo "Stringtie merge complete" $(date)

gffcompare -r ../../../refs/Panopea-generosa-vv0.74.a3-merged-2019-09-03-6-14-33.gff3 -G -o merged Pgenerosa_merged.gtf #Compute the accuracy and precision of assembly
echo "GFFcompare complete, Starting gene count matrix assembly..." $(date)

python ../../scripts/prepDE.py -g Pgenerosa_gene_count_matrix.csv -i listGTF.txt #Compile the gene count matrix
echo "Gene count matrix compiled." $(date)
```
*exit bluewaves and run from terminal*

- save the read count matrix to local pc
```
scp  samuel_gurr@bluewaves.uri.edu:/data/putnamlab/sgurr/Geoduck_TagSeq/output/stringtie/transcript_count_matrix.csv C:/C:/Users/samjg/Documents/My_Projects/Geoduck_TagSeq/Geoduck_TagSeq
```