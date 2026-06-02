Bioinformatics workflow scripts.
SPAdes comandos: spades.py \
-1 reads_R1.fastq \
-2 reads_R2.fastq \
-o spades_output
ABySS: abyss-pe \
k=64 \
name=abyss \
in='reads_R1.fastq reads_R2.fastq'
VelvetOptimisier: VelvetOptimiser.pl \
-s 21 \
-e 99 \
-f '-shortPaired -fastq -separate reads_R1.fastq reads_R2.fastq'
QUAST: quast.py \
spades.fasta \
abyss.fasta \
velvet.fasta \
-o quast_results
