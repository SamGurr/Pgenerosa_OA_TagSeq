#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=200GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=samuel_gurr@uri.edu
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"
#SBATCH -D /data/putnamlab/sgurr/refs

date

module load BLAST+/2.7.1-foss-2018a

# Pass the database you want as the first argument
# Pass the query you want as the second argument

database=$1
query=$2

mkdir ./db/

# Make the database
makeblastdb -in $1 -dbtype nucl -out ./db/database


#runs blast on the P generosa genome files against the Pacific oyster KEGG database we created above
blastn -query $2 -db ./db/database \
  evalue 20 # default evalue is 10 if you do not run this  
  -out ./${2}_out.tsv -outfmt 6 # qseqid sseqid pident evalue length qlen slen qstart qend sstart send sseq

echo "Done"
date

