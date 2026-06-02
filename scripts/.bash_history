cd /home/robertolede/bioinformatica_2026-1/Practica_6/spades/spades.log
/home/robertolede/bioinformatica_2026-1/Practica_6/spades/spades.log
ls spades
pwd
ls
find /home/robertolede -type d -name "spades"
cd /home/robertolede/bioinformatica_2026-1/Practica_6
ls
find /home/robertolede -type d -name "spades"
ls /home/robertolede/bioinformatica_2026-1/Practica_6/spades
find /home/robertolede -type d -name "spades"
find /home/robertolede -name "contigs.fasta"
cd /home/robertolede/bioinformatica_2026-1/Practica_6
rm -rf spades
spades.py --careful --isolate -1 Erw2_1.fastq -2 Erw2_2.fastq -k 21,33,55,77,99,105 -o spades
ls
cd /home/robertolede/bioinformatica_2026-1/Practica_6
pwd
ls
df -h .
spades.py --careful --isolate -1 Erw2_1.fastq -2 Erw2_2.fastq -k 21,33,55,77,99,105 -o spades 2>&1 | tee spades_run.log
cd /home/robertolede/bioinformatica_2026-1/Practica_6
spades.py --isolate -1 Erw2_1.fastq -2 Erw2_2.fastq -k 21,33,55,77,99,105 -o spades
rm -rf spades
spades.py --isolate -1 Erw2_1.fastq -2 Erw2_2.fastq -k 21,33,55,77 -o spades
free -h
for k in 55 65 75 85 95; do   echo "=== Ejecutando ABySS con k=$k ===";      mkdir -p abyss-k$k;      abyss-pe -C abyss-k$k   B=2G   k=$k   name=abyss   in='sub_1.fastq sub_2.fastq';      echo "=== Terminado k=$k ==="; done
cd /home/robertolede/bioinformatica_2026-1/Practica_6
seqkit sample -p 0.2 Erw2_1.fastq > sub_1.fastq
seqkit sample -p 0.2 Erw2_2.fastq > sub_2.fastq
sudo apt update
sudo apt install seqkit -y
seqkit sample -p 0.2 Erw2_1.fastq > sub_1.fastq
seqkit sample -p 0.2 Erw2_2.fastq > sub_2.fastq
seqkit sample -p 0.2 -s 42 Erw2_1.fastq > sub_1.fastq
seqkit sample -p 0.2 -s 42 Erw2_2.fastq > sub_2.fastq
cd /home/robertolede/bioinformatica_2026-1/Practica_6
rm -rf spades
spades.py --isolate -1 sub_1.fastq -2 sub_2.fastq -k 21,33,55,77 -m 10 -o spades
cd spades
ls
grep -v ">" contigs.fasta | wc -c
grep ">" contigs.fasta | wc -l
cd /bioinformatica_2026-1/Practica_6/
cd /home/robertolede/bioinformatica_2026-1/Practica_6/
abyss-pe
sudo apt update
sudo apt install abyss -y
cd /home/robertolede/bioinformatica_2026-1/Practica_6
ls
ls sub_*.fastq
for k in 55 65 75 85 95; do   echo "=== Ejecutando ABySS con k=$k ===";      mkdir -p abyss-k$k;      abyss-pe -C abyss-k$k   B=2G   k=$k   name=abyss   in='sub_1.fastq sub_2.fastq';      echo "=== Terminado k=$k ==="; done
in='../sub_1.fastq ../sub_2.fastq'
in='/home/robertolede/bioinformatica_2026-1/Practica_6/sub_1.fastq /home/robertolede/bioinformatica_2026-1/Practica_6/sub_2.fastq'
#!/bin/bash
for k in 55 65 75 85 95; do   echo "=== Ejecutando ABySS con k=$k ===";    mkdir -p abyss-k$k;    abyss-pe -C abyss-k$k     B=2G     k=$k     name=abyss     in='../sub_1.fastq ../sub_2.fastq';    echo "=== Terminado k=$k ==="; done
conda install -c bioconda velvetoptimiser velvet
#!/bin/bash
mkdir -p velvet_results
cd velvet_results
VelvetOptimiser.pl  -s 55 -e 95 -f "-fastq -shortPaired ../sub_1.fastq ../sub_2.fastq"  -t 2  -d output_velvet
velveth output 55 -fastq -shortPaired ../sub_1.fastq ../sub_2.fastq
velvetg output -cov_cutoff auto
VelvetOptimiser.pl  -s 21 -e 31  -f "-fastq -shortPaired ../sub_1.fastq ../sub_2.fastq"  -d velvetopt  -p velvetopt  -t 2
conda install -c bioconda velvetoptimiser velvet
VelvetOptimiser.pl --help
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install velvetoptimiser
VelvetOptimiser.pl --help
which VelvetOptimiser.pl
conda install -c bioconda velvetoptimiser velvet
VelvetOptimiser.pl --help
conda search velvetoptimiser
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels defaults
conda config --set channel_priority strict
conda install velvet velvetoptimiser -c bioconda
conda install velvetoptimiser -c bioconda -v
#!/bin/bash
mkdir -p velvet_results
for k in 21 23 25 27 29 31; do   echo "====================================";   echo " Ejecutando Velvet con k = $k";   echo "====================================";    velveth velvet_results/k$k $k -fastq -shortPaired ../sub_1.fastq ../sub_2.fastq;    velvetg velvet_results/k$k -cov_cutoff auto;  done
grep "N50" velvet_results/*/Logfile.txt
quast.py  velvet_results/k21/contigs.fa  velvet_results/k23/contigs.fa  velvet_results/k25/contigs.fa  velvet_results/k27/contigs.fa  velvet_results/k29/contigs.fa  velvet_results/k31/contigs.fa  -o quast_velvet
conda install -c bioconda quast
quast.py --version
cd ~/bioinformatica_2026-1/Practica_6
quast.py --version
cd ~/bioinformatica_2026-1/Practica_6
conda install -c bioconda quast
quast.py --version
conda create -n quast_env python=3.10 -y
conda activate quast_env
conda install -c bioconda quast -y
cd ~/bioinformatica_2026-1/Practica_6
mkdir -p assemblies/abyss
cp ../abyss-k55/abyss-6.fa assemblies/abyss/k55.fa
cp ../abyss-k65/abyss-6.fa assemblies/abyss/k65.fa
cp ../abyss-k75/abyss-6.fa assemblies/abyss/k75.fa
cp ../abyss-k85/abyss-6.fa assemblies/abyss/k85.fa
cp ../abyss-k95/abyss-6.fa assemblies/abyss/k95.fa
cd ~/bioinformatica_2026-1/Practica_6
mkdir -p assemblies/spades
cp ../spades/contigs.fasta assemblies/spades/spades.fa
ls
mkdir -p assemblies/velvet
cp velvet_results/k21/contigs.fa assemblies/velvet/k21.fa
cp velvet_results/k23/contigs.fa assemblies/velvet/k23.fa
cp velvet_results/k25/contigs.fa assemblies/velvet/k25.fa
cp velvet_results/k27/contigs.fa assemblies/velvet/k27.fa
cp velvet_results/k29/contigs.fa assemblies/velvet/k29.fa
cp velvet_results/k31/contigs.fa assemblies/velvet/k31.fa
mkdir -p assemblies/abyss
mkdir -p assemblies/spades
mkdir -p assemblies/velvet
mkdir -p resultados_quast
cp ~/bioinformatica_2026-1/Practica_6/abyss-k55/abyss-6.fa assemblies/abyss/k55.fa
cp ~/bioinformatica_2026-1/Practica_6/abyss-k65/abyss-6.fa assemblies/abyss/k65.fa
cp ~/bioinformatica_2026-1/Practica_6/abyss-k75/abyss-6.fa assemblies/abyss/k75.fa
cp ~/bioinformatica_2026-1/Practica_6/abyss-k85/abyss-6.fa assemblies/abyss/k85.fa
cp ~/bioinformatica_2026-1/Practica_6/abyss-k95/abyss-6.fa assemblies/abyss/k95.fa
cp ~/bioinformatica_2026-1/Practica_6/spades/contigs.fasta assemblies/spades/spades.fa
find ~/bioinformatica_2026-1/Practica_6 -name "contigs.fa"
cp $(find ~/bioinformatica_2026-1/Practica_6 -path "*k21*/contigs.fa") assemblies/velvet/k21.fa
cp $(find ~/bioinformatica_2026-1/Practica_6 -path "*k23*/contigs.fa") assemblies/velvet/k23.fa
cp $(find ~/bioinformatica_2026-1/Practica_6 -path "*k25*/contigs.fa") assemblies/velvet/k25.fa
cp $(find ~/bioinformatica_2026-1/Practica_6 -path "*k27*/contigs.fa") assemblies/velvet/k27.fa
cp $(find ~/bioinformatica_2026-1/Practica_6 -path "*k29*/contigs.fa") assemblies/velvet/k29.fa
cp $(find ~/bioinformatica_2026-1/Practica_6 -path "*k31*/contigs.fa") assemblies/velvet/k31.fa
cd ~/bioinformatica_2026-1/Practica_6
quast.py  assemblies/abyss/*.fa  assemblies/spades/spades.fa  assemblies/velvet/*.fa  -o resultados_quast
conda create -n quast_env python=3.10 -y
conda activate quast_env
conda install -c bioconda quast -y
quast.py --version
which quast.py
conda info --envs
conda list | grep quast
conda activate quast_env
conda install -c bioconda -c conda-forge quast -y
conda list | grep quast
conda list | grep -i quast
conda install -c bioconda -c conda-forge quast -y
cd ~/bioinformatica_2026-1/Practica_6
python -m quast  assemblies/abyss/*.fa  assemblies/spades/spades.fa  assemblies/velvet/*.fa  -o resultados_quast
conda activate quast_env
conda install -c bioconda -c conda-forge quast -y
conda list | grep quast
conda install -c bioconda quast=5.2.0 -y
cd ~/bioinformatica_2026-1/Practica_6
python -m quast  assemblies/abyss/*.fa  assemblies/spades/spades.fa  assemblies/velvet/*.fa  -o resultados_quast
conda activate quast_env
conda install -c bioconda quast -y
python -m quast --version
conda activate quast_env
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install quast=5.2.0 -y
quast.py --version
cd ~/bioinformatica_2026-1/Practica_6
quast.py  assemblies/abyss/*.fa  assemblies/spades/spades.fa  assemblies/velvet/*.fa  -o resultados_quast
xdg-open resultados_quast/report.html
sudo apt install xdg-utils
xdg-open resultados_quast/report.html
cd ~/bioinformatica_2026-1/Practica_6/resultados_quast
ls
google-chrome report.html
cd ~/bioinformatica_2026-1/Practica_6/resultados_quast
cat report.tsv
firefox ~/bioinformatica_2026-1/Practica_6/resultados_quast/report.html
xdg-open ~/bioinformatica_2026-1/Practica_6/resultados_quast/report.html
firefox ~/bioinformatica_2026-1/Practica_6/resultados_quast/report.html
cp ~/bioinformatica_2026-1/Practica_6/resultados_quast/report.html .
cd ~/bioinformatica_2026-1/Practica_6/resultados_quast
column -t report.tsv | less
cat ~/bioinformatica_2026-1/Practica_6/resultados_quast/report.tsv
quast.py abyss_k95.fasta spades_contigs.fasta velvet_k31.fasta -o quast_results
awk '/^>/ {if (seqlen){print seqlen}; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' archivo.fasta | sort -n
ls
ls -lh ~/bioinformatica_2026-1/Practica_6/spades
find ~/bioinformatica_2026-1/Practica_6/spades -name "contigs.fasta"
spades/contigs.fasta
pwd
ls -lh
cd ~/bioinformatica_2026-1/Practica_6
ls -R
find . -maxdepth 2
ls -lh contigs.fasta
cd ~/bioinformatica_2026-1/Practica_6/spades
ls -lh contigs.fasta
find ~/bioinformatica_2026-1/Practica_6/spades -name "contigs.fasta"
cd ~/bioinformatica_2026-1/Practica_6
mkdir -p assemblies/abyss
mkdir -p assemblies/velvet
mkdir -p assemblies/spades
mkdir -p resultados_quast
cp ~/bioinformatica_2026-1/Practica_6/abyss-k55/abyss-6.fa assemblies/abyss/k55.fa
cp ~/bioinformatica_2026-1/Practica_6/abyss-k65/abyss-6.fa assemblies/abyss/k65.fa
cp ~/bioinformatica_2026-1/Practica_6/abyss-k75/abyss-6.fa assemblies/abyss/k75.fa
cp ~/bioinformatica_2026-1/Practica_6/abyss-k85/abyss-6.fa assemblies/abyss/k85.fa
cp ~/bioinformatica_2026-1/Practica_6/abyss-k95/abyss-6.fa assemblies/abyss/k95.fa
cp ~/bioinformatica_2026-1/Practica_6/velvet_results/k21/contigs.fa assemblies/velvet/k21.fa
cp ~/bioinformatica_6/Practica_6/velvet_results/k23/contigs.fa assemblies/velvet/k23.fa
cp ~/bioinformatica_6/Practica_6/velvet_results/k25/contigs.fa assemblies/velvet/k25.fa
cp ~/bioinformatica_6/Practica_6/velvet_results/k27/contigs.fa assemblies/velvet/k27.fa
cp ~/bioinformatica_6/Practica_6/velvet_results/k29/contigs.fa assemblies/velvet/k29.fa
cp ~/bioinformatica_6/Practica_6/velvet_results/k31/contigs.fa assemblies/velvet/k31.fa
cp ~/bioinformatica_2026-1/Practica_6/spades/contigs.fasta assemblies/spades/spades.fa
cd ~/bioinformatica_2026-1/Practica_6
quast.py  assemblies/abyss/*.fa  assemblies/velvet/*.fa  assemblies/spades/spades.fa  -o resultados_quast
cd resultados_quast
column -t report.tsv | less
cd ~/bioinformatica_2026-1/Practica_6/resultados_quast
cp ~/bioinformatica_2026-1/Practica_6/resultados_quast/report.html .
code report.html
column -t report.tsv | less
cp ~/bioinformatica_2026-1/Practica_6/resultados_quast/report.html .
cd ~/bioinformatica_2026-1/Practica_6
cp resultados_quast/report.html .
column -t resultados_quast/report.tsv | less
column -t report.txt
grep -E "Assembly|# contigs|Total length|GC|N50|L50" report.txt | sed 's/  */\t/g' | column -t
cd ~/bioinformatica_2026-1/Practica_6
mkdir -p prokka_spades
prokka  --outdir prokka_spades  --force  --prefix LMX  --genus Erwinia  --species amylovora  --kingdom Bacteria  --usegenus  --evalue 1e-12  assemblies/spades/spades.fa
ls prokka_spades
prokka --setupdb
prokka  --outdir prokka_spades  --force  --prefix LMX  --genus Erwinia  --species amylovora  --kingdom Bacteria  --usegenus  --evalue 1e-12  assemblies/spades/spades.fa
cd ~/bioinformatica_2026-1/Practica_6/prokka_spades
ls -lh
ls LMX.*
grep -c "CDS" LMX.gff
prokka --setupdb
prokka --version
cd ~/bioinformatica_2026-1/Practica_6
prokka  --outdir prokka_spades  --force  --prefix LMX  --genus Erwinia  --species amylovora  --kingdom Bacteria  --usegenus  --evalue 1e-12  assemblies/spades/spades.fa
ls prokka_spades
ps aux | grep prokka
ls ~/miniconda3/db
ls ~/yes/db
cd ~/bioinformatica_2026-1/Practica_6
pwd
conda activate quast_env
cd ~/bioinformatica_2026-1/Practica_6
prokka --setupdb
prokka  --outdir prokka_spades  --force  --prefix LMX  --genus Erwinia  --species amylovora  --kingdom Bacteria  --usegenus  --evalue 1e-12  assemblies/spades/spades.fa
ls prokka_spades
pwd
ls -lh
find ~/bioinformatica_2026-1/Practica_6 -type d
find ~/bioinformatica_2026-1/Practica_6 -name "contigs.fasta"
pwd && ls -lh
cd ~/bioinformatica_2026-1/Practica_6
ls spades | head
ls spades/contigs.fasta
conda activate quast_env
prokka --setupdb
prokka  --outdir prokka_spades  --force  --prefix Erw2  --genus Erwinia  --species amylovora  --kingdom Bacteria  --usegenus  --evalue 1e-12  spades/contigs.fasta
ls prokka_spades
prokka --version
cd ~/bioinformatica_2026-1/Practica_6
prokka  --outdir prokka_spades  --force  --prefix Erw2  --genus Erwinia  --species amylovora  --kingdom Bacteria  --usegenus  --evalue 1e-12  spades/contigs.fasta
prokka --setupdb
cd ~/bioinformatica_2026-1/Practica_6
prokka  --outdir prokka_spades  --force  --prefix Erw2  --genus Erwinia  --species amylovora  --kingdom Bacteria  --usegenus  --evalue 1e-12  spades/contigs.fasta
ls prokka_spades
~/.prokka
rm -rf ~/.prokka
rm -rf ~/.cache/prokka
which prokka
conda install -c bioconda prokka --force-reinstall
cd ~/bioinformatica_2026-1/Practica_6
prokka  --outdir prokka_spades  --force  --prefix Erw2  --genus Erwinia  --species amylovora  --kingdom Bacteria  spades/contigs.fasta
conda list | grep prokka
which prokka
prokka --version
echo $CONDA_DEFAULT_ENV
conda create -n prokka_env -c bioconda prokka -y
conda activate prokka_env
which prokka
prokka --version
cd ~/bioinformatica_2026-1/Practica_6
prokka  --outdir prokka_spades  --force  --prefix Erw2  --genus Erwinia  --species amylovora  --kingdom Bacteria  spades/contigs.fasta
blastp -version
cd ~/bioinformatica_2026-1/Practica_6
prokka  --outdir prokka_spades  --force  --prefix Erw2  --genus Erwinia  --species amylovora  --kingdom Bacteria  spades/contigs.fasta
blastp -version
conda deactivate
conda remove -n prokka_env --all -y
conda create -n prokka_env -c bioconda -c conda-forge prokka blast hmmer -y
conda activate prokka_env
prokka --version
blastp -version
cd ~/bioinformatica_2026-1/Practica_6
prokka  --outdir prokka_spades  --force  --prefix Erw2  --genus Erwinia  --species amylovora  --kingdom Bacteria  spades/contigs.fasta
which prokka
which blastp
blastp -version
conda deactivate
conda remove -n prokka_env --all -y
conda clean --all -y
conda create -n prokka_env -c bioconda -c conda-forge prokka blast hmmer aragorn barrnap -y
conda activate prokka_env
which prokka
prokka --version
blastp -version
cd ~/bioinformatica_2026-1/Practica_6
prokka  --outdir prokka_spades  --force  --prefix Erw2  --genus Erwinia  --species amylovora  --kingdom Bacteria  spades/contigs.fasta
which -a prokka
which -a blastp
conda activate prokka_env
which prokka
cd ~/bioinformatica_2026-1/Practica_6
prokka  --outdir prokka_spades  --force  --prefix Erw2  --genus Erwinia  --species amylovora  --kingdom Bacteria  spades/contigs.fasta
which prokka
conda deactivate
conda remove -n prokka_env --all -y
conda create -n prokka_env -c conda-forge -c bioconda prokka=1.14 blast hmmer aragorn barrnap -y
conda activate prokka_env
prokka --version
blastp -version
cd ~/bioinformatica_2026-1/Practica_6
prokka  --outdir prokka_spades  --force  --prefix Erw2  --genus Erwinia  --species amylovora  --kingdom Bacteria  spades/contigs.fasta
cd ~/bioinformatica_2026-1/Practica_6/prokka_spades
ls -lh
head Erw2.gff
grep -c "CDS" Erw2.gff
head Erw2.faa
column -t Erw2.tsv | less
head Erw2.tsv
ls *.fna | wc -l
cd /home/robertolede/Base_de_datos
cd /home/robertolede/bioinformatica_2026-1/Base_de_datos
ls *.fna | wc -l
cd /home/robertolede/bioinformatica_2026-1/Base_de_datos/ncbi_dataset/datos/
cd /home/robertolede/bioinformatica_2026-1/Base_de_datos/ncbi_dataset/
cd datos
cd data
ls *.fna
cd /home/robertolede/bioinformatica_2026-1/Base_de_datos/solo_20
ls *.fna | wc -l
ls *.fna
ls *.fna | head -n 20 > lista_20.txt
ls *.fna | wc -l
cat *.fna > Erwinia_20_genomas.fasta
ls
grep ">" Erwinia_20_genomas.fasta | wc -l
mkdir /home/robertolede/bioinformatica_2026-1/Base_de_datos/solo_20_filtrado
for file in /home/robertolede/bioinformatica_2026-1/Base_de_datos/solo_20/*.fna; do   awk '/^>/{if(seen++) exit} {print}' "$file" > /home/robertolede/bioinformatica_2026-1/Base_de_datos/solo_20_filtrado/$(basename "$file"); done
cd /home/robertolede/bioinformatica_2026-1/Base_de_datos/solo_20_filtrado
cat *.fna > Erwinia_20_genomas.fasta
grep ">" Erwinia_20_genomas.fasta | wc -l
conda activate quast_env
cd /home/robertolede/bioinformatica_2026-1/Practica_6
quast.py /home/robertolede/bioinformatica_2026-1/Practica_6/spades/contigs.fasta /home/robertolede/bioinformatica_2026-1/Practica_6/abyss-k85/abyss-contigs.fa /home/robertolede/bioinformatica_2026-1/Practica_6/velvetopt/contigs.fa --labels SPAdes,ABySS_k85,Velvet_k31 -o /home/robertolede/bioinformatica_2026-1/Practica_6/quast_comparacion_final
ls /home/robertolede/bioinformatica_2026-1/Practica_6/velvetopt
find /home/robertolede/bioinformatica_2026-1/Practica_6/velvetopt -name "*.fa"
ls
find /home/robertolede/bioinformatica_2026-1/Practica_6/velvet_results -name "contigs.fa"
conda activate quast_env
quast.py /home/robertolede/bioinformatica_2026-1/Practica_6/spades/contigs.fasta /home/robertolede/bioinformatica_2026-1/Practica_6/abyss-k85/abyss-contigs.fa /home/robertolede/bioinformatica_2026-1/Practica_6/velvet_results/auto_data_31/contigs.fa --labels SPAdes,ABySS_k85,Velvet_k31 -o /home/robertolede/bioinformatica_2026-1/Practica_6/quast_comparacion_final
find /home/robertolede/bioinformatica_2026-1/Practica_6 -name "abyss-6.fa"
conda activate quast_env
quast.py /home/robertolede/bioinformatica_2026-1/Practica_6/spades/contigs.fasta /home/robertolede/bioinformatica_2026-1/Practica_6/abyss-k85/abyss-6.fa /home/robertolede/bioinformatica_2026-1/Practica_6/velvet_results/auto_data_31/contigs.fa --labels SPAdes,ABySS_k85,Velvet_k31 -o /home/robertolede/bioinformatica_2026-1/Practica_6/quast_comparacion_final
conda activate quast_env
quast.py /home/robertolede/bioinformatica_2026-1/Practica_6/spades/contigs.fasta /home/robertolede/bioinformatica_2026-1/Practica_6/abyss-k85/abyss-6.fa /home/robertolede/bioinformatica_2026-1/Practica_6/velvet_results/auto_data_31/contigs.fa --labels SPAdes,ABySS_k85,Velvet_k31 -o /home/robertolede/bioinformatica_2026-1/Practica_6/quast_comparacion_final
WARNING: Python locale settings can't be changed
ERROR! File not found (contigs): /home/robertolede/bioinformatica_2026-1/Practica_6/velvet_results/auto_data_31/contigs.fa

In case you have troubles running QUAST, you can write to quast.support@cab.spbu.ru
or report an issue on our GitHub repository https://github.com/ablab/quast/issues
Please provide us with quast.log file from the output directory.

find /home/robertolede/bioinformatica_2026-1/Practica_6/velvet_results -name "*.fa"
conda activate quast_env
quast.py /home/robertolede/bioinformatica_2026-1/Practica_6/spades/contigs.fasta /home/robertolede/bioinformatica_2026-1/Practica_6/abyss-k85/abyss-6.fa /home/robertolede/bioinformatica_2026-1/Practica_6/velvet_results/auto_data_31/contigs.fa --labels SPAdes,ABySS_k85,Velvet_k31 -o /home/robertolede/bioinformatica_2026-1/Practica_6/quast_comparacion_final
cd /home/robertolede/bioinformatica_2026-1/Practica_6
conda activate quast_env
ls spades/contigs.fasta
ls abyss-k85/abyss-6.fa
ls velvet_results/velvet_results/k31/contigs.fa
quast.py spades/contigs.fasta abyss-k85/abyss-6.fa velvet_results/velvet_results/k31/contigs.fa --labels SPAdes,ABySS_k85,Velvet_k31 -o quast_comparacion_final
quast.py spades/contigs.fasta abyss-k85/abyss-6.fa velvet_results/velvet_results/k31/contigs.fa --labels SPAdes,ABySS_k85,Velvet_k31 -o quast_comparacion_final
unzip GB_Final_U8972884G0
cd
cd Descargas
cd /Descargas
cd ~/Descargas
C:\Users\rober\Downloads\GB_Final_U8972884G0.zip
unzip GB_Final_U8972884G0.zip
cd ~/GB_Final_U8972884G0
cd /home/robertolede/GB_Final_U8972884G0
cd /home/robertolede/GB_Final_U8972884G0/
pdw
pwd
ls
unzip GB_Final_U8972884G0.zip 
ls /home/robertolede
unzip /home/robertolede/GB_Final_U8972884G0.zip
unzip -o /home/robertolede/GB_Final_U8972884G0.zip
ls /home/robertolede
ls
unzip GB_Final_U8972884G0.zip
ls
cd /home/robertolede/bioinformatica_2026-1/Practica_7
unzip Genomas_ref-20260511T225219Z-3-001.zip
pip install pandas bottleneck tqdm numpy biopython scikit-learn gensim
git clone https://github.com/alfonsomt/vamphyre.git
cd vamphyre/
chmod +x Install_VAMPhyRE.sh
./Install_VAMPhyRE.sh
clear
ls
cd ..
ls
clear
ls
cd VAMPhyRE/
ls
cd ..
cd vamphyre/
ls
cd ..
ls
cd VAMPhyRE/
ls
clear
ls
cd bin/
ls
./VH5cmdl 
clear
ls
cd ..
ls
cd Genomas_ref/
ls
cd ..
ñs
ls
mv spades.fa Genomas_ref
cd Genomas_ref/
ls
grep ">" -c *.fasta
ls
cd ..
ls
nano ContigsConcatenation.sh
clear
nano ContigsConcatenation.sh
cd Genomas_ref/
nano ContigsConcatenation.sh
clear
ls
chmod +x ContigsConcatenation.sh
./ContigsConcatenation.sh 
ls
mv spades.fa spades.fasta
./ContigsConcatenation.sh 
ls
rm *.join.join.fasta
ls
grep ">" -c *.join.fasta
ls
clear
ls
mv VAMPhyRE/ bioinformatica_2026-1/
cd bioinformatica_2026-1/
ls
mv VAMPhyRE/ Practica_7/
ls
cd bioinformatica_2026-1/
ls
cd Practica_7/
ls
clear
ls
cd Genomas_ref/
ls
mkdir originales
mkdir join
mv *.join.fasta 
mv *.join.fasta  join/
ls
mv *.fasta originales/
ls
cd join/
ls
mv * ../
cd ..
ls
rm -r join/
clear
ls
ls -lh
clear
ls
ls *.fasta > list.txt
ld
ls
nano list.txt 
clear
ls
cd VAMPhyRE/
ls
cd bin/
ls
pwd
cd ..
ls
/home/robertolede/bioinformatica_2026-1/Practica_7/Genomas_ref/VAMPhyRE/bin/VH5cmdl 
/home/robertolede/bioinformatica_2026-1/Practica_7/Genomas_ref/VAMPhyRE/bin/VH5cmdl -PROBEFILE vps13.txt -TARGETLIST list.txt -OUTFILE vh_salida.txt -MISMATCHES 1 -STRAND both
/home/robertolede/bioinformatica_2026-1/Practica_7/Genomas_ref/VAMPhyRE/bin/VFAT -VHFILE vh_salida.txt -TARGETLIST list.txt -OUTFILE dist_salida -LEFTEXT 4 -RIGHTEXT 4 -THRESHOLD 19 -MODE DISTANCE
history > historial.txt
cd Practica_6/
pwd
cd /home/robertolede/bioinformatica_2026-1/Practica_6/assemblies/spades/
./ContigsConcatenation.sh
spades./ContigsConcatenation.sh
chmod +spades
chmod +x
/home/robertolede/bioinformatica_2026-1/Practica_6/assemblies/spades
/home/robertolede/bioinformatica_2026-1/Practica_6/assemblies/spades/
cd /home/robertolede/bioinformatica_2026-1/Practica_6/assemblies/spades
conda activate vamphyre
conda env list
conda info --envs
conda env list
conda activate Vamphyre
find ~/ -iname "*vamphyre*" 2>/dev/null
cd VAMPhyRE
cd vamphyre
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto
ls
cd VAMPhyRE
ls -1
chmod +x VH5cmdl
chmod +x VFAT
bin/
cd bin/
ls
chmod +x VH5cmdl
chmod +x VFAT
chmod +x *
./VH5cmdl -h
./VH5cmdl --help
./VH5cmdl + x spades
../
cd /bioinformatica_2026-1/Practica_7_proyecto/
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto
ls
cd VAMPhyRE
ls
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/VAMPhyRE/bin
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/VAMPhyRE/bin
ls
chmod +x VH5cmdl
chmod +x VFAT
./VH5cmdl +h
./VH5cmdl
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto
mkdir analisis_vamphyre
cd analisis_vamphyre
cp /home/robertolede/bioinformatica_2026-1/Practica_6/assemblies/spades/contigs.fasta .
cp /home/robertolede/bioinformatica_2026-1/Practica_6/assemblies/spades/spades
cp /home/robertolede/bioinformatica_2026-1/Practica_6/assemblies/spades/
cd /home/robertolede/bioinformatica_2026-1/Practica_6/assemblies/spades
ls
cd /home/robertolede/bioinformatica_2026-1/Practica_6/assemblies/spades
mv spades.fa spades.fasta
ls
cp spades.fasta /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
ls
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
ls
cp /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genoma20_Erw2_grupo_externo/* .
cp /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genoma20_Erw2_grupo_externo/
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/VAMPhyRE/bin
./VH5cmdl +h
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
ls
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/VAMPhyRE/bin/VH5cmdl +x spades.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
ls
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/VAMPhyRE/bin/VH5cmdl +x spades.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/VAMPhyRE/bin/VH5cmdl +x spades.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/VAMPhyRE/bin/VH5cmdl spades.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/VAMPhyRE/bin
./VH5cmdl
help
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/VAMPhyRE/bin
ls -l VH5cmdl
./VH5cmdl
pwd
./VH5cmdl spades.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
./VH5cmdl spades.fasta
cd /home/robertolede/bioinformatica_2026-1
ls
cd Practica_7_proyecto
ls
cd VAMPhyRE
ls
cd bin
ls
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
ls
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/VAMPhyRE/bin/VH5cmdl spades.fasta
ls
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
ls
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/VAMPhyRE/bin/VH5cmdl Genoma20_Erw2_grupo_externo.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
ls
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
rm "Genoma20_Erw2_grupo_externo.fasta:Zone.Identifier"
head Genoma20_Erw2_grupo_externo.fasta
dos2unix Genoma20_Erw2_grupo_externo.fasta
sudo apt install dos2unix
dos2unix Genoma20_Erw2_grupo_externo.fasta
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/VAMPhyRE/bin/VH5cmdl Genoma20_Erw2_grupo_externo.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
mkdir separados
awk '/^>/{f="separados/genoma"++i".fasta"} {print > f}' Genoma20_Erw2_grupo_externo.fasta
ls separados
cd separados
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/VAMPhyRE/bin/VH5cmdl genoma1.fasta
head -5 genoma1.fasta
grep ">" genoma1.fasta
wc -l genoma1.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
rm -r separados
mkdir separados
csplit -z -f separados/genoma_ Genoma20_Erw2_grupo_externo.fasta '/^>/' '{*}'
ls separados
cd separados
head genoma_00
mv genoma_00 genoma1.fasta
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/VAMPhyRE/bin/VH5cmdl genoma1.fasta
head -20 genoma1.fasta
file genoma1.fasta
head -20 genoma1.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/VAMPhyRE/bin
./VH5cmdl
genoma1.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/VAMPhyRE/bin
strings VH5cmdl | grep fasta
strings VH5cmdl | grep ".vhm"
strings VH5cmdl | grep usage
./VH5cmdl -i /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/separados/genoma1.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
pwd
conda create -n Vamphyre perl -y
conda activate Vamphyre
echo $CONDA_DEFAULT_ENV
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
ls
perl /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/ContigsConcatenation.pl spades
ls
perl /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/ContigsConcatenation.pl Genoma20_Erw2_grupo_externo
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
ls
perl /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/ContigsConcatenation.pl spades.fasta
ls
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
ls
perl /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/ContigsConcatenation.pl spades.fasta
perl /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/ContigsConcatenation.pl Genoma20_Erw2_grupo_externo.fasta
ls
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
ls
chmod +x ContigsConcatenation.pl
chmod +x ContigsConcatenation.sh
cd originales
perl ../ContigsConcatenation.pl spades.fasta > spades.join
perl ../ContigsConcatenation.pl Genoma20_Erw2_grupo_externo.fasta > Genoma20_Erw2_grupo_externo.join
ls
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin/VH5cmdl spades.join
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
ls
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin/VH5cmdl spades.join.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin
chmod +x VH5cmdl
ls -l VH5cmdl
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin/VH5cmdl spades.join.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
ls
rm -f *.join *.vhm *.vf *.txt *.tmp *.log
ls
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
ls
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
ls -lh
perl ../ContigsConcatenation.pl spades.fasta
perl ../ContigsConcatenation.pl Genoma20_Erw2_grupo_externo.fasta
ls
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
perl ContigsConcatenation.pl originales/spades.fasta > originales/spades.join
perl ContigsConcatenation.pl originales/Genoma20_Erw2_grupo_externo.fasta > originales/Genoma20_Erw2_grupo_externo.join
ls originales
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin/VH5cmdl originales/spades.join
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin/VH5cmdl spades.fasta
ls *.fasta > lista.txt
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin/VH5cmdl lista.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin
./VH5cmdl
./VH5cmdl --help
./VH5cmdl -h
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
bash ContigsConcatenation.sh
find . -name "*.vh*" -o -name "*.vhm" -o -name "*.out"
bash /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/ContigsConcatenation.sh
ls
cd originales
ls
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin/VH5cmdl spades.join
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
ls *.fasta
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin/VH5cmdl spades.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
bash ContigsConcatenation.sh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
rm -f *.join *.vhm *.txt *.out
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
bash ContigsConcatenation.sh
find . -type f | grep -E "vhm|vf|out|join"
bash ContigsConcatenation.sh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin/VH5cmdl spades.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin
./VH5cmdl
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
spades.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
ls
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
ls
nano targets.txt
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin/VH5cmdl -PROBEFILE spades.fasta -TARGETLIST targets.txt -OUTFILE resultado_vamphyre.txt -MISMATCHES 0 -STRAND direct
head -n 20000 spades.fasta > spades_small.fasta
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin/VH5cmdl -PROBEFILE spades_small.fasta -TARGETLIST targets.txt -OUTFILE resultado_vamphyre.txt -MISMATCHES 0 -STRAND direct
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin/VH5cmdl -PROBEFILE spades.fasta -TARGETLIST targets.txt -OUTFILE resultado_vamphyre.txt -MISMATCHES 1 -STRAND direct
head -n 20000 spades.fasta > spades_small.fasta
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin/VH5cmdl -PROBEFILE spades_small.fasta -TARGETLIST targets.txt -OUTFILE resultado_vamphyre.txt -MISMATCHES 1 -STRAND direct -RESULTSFORMAT reduced
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
cp spades.fasta spades.join
cp Genoma20_Erw2_grupo_externo.fasta Genoma20_Erw2_grupo_externo.join
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin/VH5cmdl -PROBEFILE spades.join -TARGETLIST targets.txt -OUTFILE resultado_vamphyre.txt -MISMATCHES 1 -STRAND direct -RESULTSFORMAT reduced
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
awk '/^>/ {print ">seq_" ++i; next} {print}' spades.fasta > spades_clean.fasta
awk '/^>/ {print ">ref_" ++i; next} {print}' Genoma20_Erw2_grupo_externo.fasta > ref_clean.fasta
echo "/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/ref_clean.fasta" > targets.txt
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin/VH5cmdl -PROBEFILE spades_clean.fasta -TARGETLIST targets.txt -OUTFILE resultado_vamphyre.txt -MISMATCHES 1 -STRAND direct -RESULTSFORMAT reduced
head -n 5000 spades_clean.fasta > spades_small.fasta
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin/VH5cmdl -PROBEFILE spades_small.fasta -TARGETLIST targets.txt -OUTFILE resultado_vamphyre.txt -MISMATCHES 1 -STRAND direct -RESULTSFORMAT reduced
echo -e ">test\nATGC" > test.fasta
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin/VH5cmdl -PROBEFILE test.fasta -TARGETLIST targets.txt -OUTFILE test_out.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
seqkit head -n 40 spades.fasta > spades_20contigs.fasta
grep -c ">" spades_20contigs.fasta
seqkit head -n 20 spades.fasta > spades_20contigs.fasta
grep -c ">" spades_20contigs.fasta
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin/VH5cmdl -PROBEFILE spades_20contigs.fasta -TARGETLIST targets.txt -OUTFILE resultado_vamphyre.txt -MISMATCHES 1 -STRAND direct -RESULTSFORMAT reduced
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
cd originales
ls
rm -f spades_small* spades_20* resultado_vamphyre.txt targets.txt ref_clean.fasta spades_clean.fasta
awk '/^>/ {print ">seq_" ++i; next} {print}' spades.fasta > spades_clean.fasta
awk '/^>/ {print ">ref_" ++i; next} {print}' Genoma20_Erw2_grupo_externo.fasta > ref_clean.fasta
echo "/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/ref_clean.fasta" > targets.txt
grep -c ">" spades_clean.fasta
grep -c ">" ref_clean.fasta
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin/VH5cmdl -PROBEFILE spades_clean.fasta -TARGETLIST targets.txt -OUTFILE resultado_vamphyre.txt -MISMATCHES 1 -STRAND direct -RESULTSFORMAT reduced
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
ls -lh
cd originales
ls -lh
grep -c ">" spades.fasta
awk '/^>/ {if (c==20) exit; c++; print; next} {print}' spades.fasta > spades_20.fasta
grep -c ">" spades_20.fasta
echo "/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/Genoma20_Erw2_grupo_externo.fasta" > lista.txt
ls -lh spades_20.fasta Genoma20_Erw2_grupo_externo.fasta lista.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin
cd ../VPS
ls -lh
cd ../bin
./VH5cmdl -PROBEFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/spades_20.fasta -TARGETLIST /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/lista.txt -OUTFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/resultado_vamphyre.txt -MISMATCHES 1 -STRAND both -RESULTSFORMAT reduced
head -n 10 spades_20.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
ls -lh
echo "Ambiente VAMPhyRE activo en: $(pwd)"
ls -lh VAMPhyRE/bin/VH5cmdl
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
ls -lh
grep -c ">" spades.fasta
awk '/^>/ {if (c==20) exit; c++; print; next} {print}' spades.fasta > spades_20.fasta
grep -c ">" spades_20.fasta
echo "/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/Genoma20_Erw2_grupo_externo.fasta" > lista.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin
./VH5cmdl -PROBEFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/spades_20.fasta -TARGETLIST /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/lista.txt -OUTFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/resultado_vamphyre.txt -MISMATCHES 1 -STRAND both -RESULTSFORMAT reduced
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
cd originales
ls -lh
grep -c ">" spades.fasta
awk '
BEGIN{c=0}
$0 ~ /^>/ {c++; if(c>25) exit}
{print}
' spades.fasta > spades_25.fasta
grep -c ">" spades_25.fasta
echo "/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/Genoma20_Erw2_grupo_externo.fasta" > lista.txt
awk '/^>/ {print ">seq_" ++i; next} {print}' spades_25.fasta > spades_vh5.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin
./VH5cmdl -PROBEFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/spades_vh5.fasta -TARGETLIST /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/lista.txt -OUTFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/resultado_vamphyre.txt -MISMATCHES 1 -STRAND both -RESULTSFORMAT reduced
# 1. Subsample de SPAdes aún más agresivo
awk '/^>/ {if (c==10) exit; c++; print; next} {print}' spades.fasta > spades_min.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
cd originales
ls -lh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
ls -lh
cd ~/bioinformatica_2026-1/Practica_7
pwd
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
ls -lh
ls -lh originales
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
rm -f spades_20.fasta spades_25.fasta spades_vh5.fasta resultado_vamphyre.txt
ls -lh
grep -c ">" Genoma20_Erw2_grupo_externo.fasta
grep ">" Genoma20_Erw2_grupo_externo.fasta
ls -lh Genoma20_Erw2_grupo_externo.fasta
head -n 20 Genoma20_Erw2_grupo_externo.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
rm -f originales/spades_20.fasta originales/spades_25.fasta originales/spades_vh5.fasta originales/resultado_vamphyre.txt
ls -lh originales
grep -c ">" originales/Genoma20_Erw2_grupo_externo.fasta
awk '
BEGIN{c=0}
$0 ~ /^>/ {c++; if(c>15) exit; print; next}
{print}
' originales/spades.fasta > originales/spades_vh5.fasta
grep -c ">" originales/spades_vh5.fasta
awk '/^>/ {print ">seq_" ++i; next} {print}' originales/spades_vh5.fasta > originales/spades_vh5_clean.fasta
echo "/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/Genoma20_Erw2_grupo_externo.fasta" > originales/lista.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin
./VH5cmdl -PROBEFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/spades_vh5_clean.fasta -TARGETLIST /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/lista.txt -OUTFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/resultado_vamphyre.txt -MISMATCHES 1 -STRAND direct -RESULTSFORMAT reduced
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin
chmod +x VH5cmdl
ls -lh VH5cmdl
./VH5cmdl
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
awk '
BEGIN{c=0}
$0 ~ /^>/ {c++; if(c>5) exit; print; next}
{print}
' spades.fasta > spades_vh5_safe.fasta
awk '/^>/ {print ">seq_" ++i; next} {print}' spades_vh5_safe.fasta > spades_vh5_final.fasta
grep -c ">" spades_vh5_final.fasta
echo "/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/Genoma20_Erw2_grupo_externo.fasta" > lista.txt
./VH5cmdl -PROBEFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/spades_vh5_final.fasta -TARGETLIST /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/lista.txt -OUTFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/resultado_vamphyre.txt -MISMATCHES 1 -STRAND direct -RESULTSFORMAT reduced
./VH5cmdl -PROBEFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/spades_vh5_final.fasta -TARGETLIST /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/lista.txt -OUTFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/resultado_vamphyre.txt -MISMATCHES 1 -STRAND direct -RESULTSFORMAT reduced
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin
ls -lh
chmod +x VH5cmdl
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin
pwd
ls -lh VH5cmdl
chmod +x VH5cmdl
./VH5cmdl
./VH5cmdl -PROBEFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/spades_vh5_final.fasta -TARGETLIST /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/lista.txt -OUTFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/resultado_vamphyre.txt -MISMATCHES 1 -STRAND direct -RESULTSFORMAT reduced
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
pwd
cd originales
ls -lh
rm -f spades_vh5.fasta spades_vh5_clean.fasta spades_vh5_final.fasta spades_vh5_safe.fasta resultado_vamphyre.txt
ls -lh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
mash sketch -o spades spades.fasta
sudo apt update
sudo apt install mash -y
mash --help
conda create -n mash_env -c bioconda mash -y
conda activate mash_env
mash --help
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
pwd
ls -lh
conda deactivate
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
pwd
ls -lh
grep -c ">" Genoma20_Erw2_grupo_externo.fasta
grep -c ">" spades.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
awk '
BEGIN{c=0}
$0 ~ /^>/ {c++; if(c>20) exit; print; next}
{print}
' spades.fasta > spades_vh5.fasta
awk '/^>/ {print ">seq_" ++i; next} {print}' spades_vh5.fasta > spades_vh5_clean.fasta
grep -c ">" spades_vh5_clean.fasta
echo "/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/Genoma20_Erw2_grupo_externo.fasta" > lista.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin
./VH5cmdl -PROBEFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/spades_vh5_clean.fasta -TARGETLIST /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/lista.txt -OUTFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/resultado_vamphyre.txt -MISMATCHES 1 -STRAND both -RESULTSFORMAT reduced
./VH5cmdl -PROBEFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/spades_vh5_clean.fasta -TARGETLIST /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/lista.txt -OUTFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales/resultado_vamphyre.txt -MISMATCHES 1 -STRAND both -RESULTSFORMAT reduced
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
pwd
ls -lh
awk '/^>/ {print ">Query_Spades"; next} {print}' spades.fasta > spades_mega.fasta
awk '/^>/ {print ">Ref_" ++i; next} {print}' Genoma20_Erw2_grupo_externo.fasta > ref_mega.fasta
grep -c ">" spades_mega.fasta
grep -c ">" ref_mega.fasta
cat spades_mega.fasta ref_mega.fasta > dataset_mega.fasta
grep -c ">" dataset_mega.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
awk '
BEGIN{seq=""}
/^>/ {
if(seq!="") print seq
print ">Query_Spades"
seq=""
next
}
{seq=seq $0}
END{print seq}
' spades.fasta > spades_single.fasta
grep -c ">" spades_single.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
awk '
BEGIN{seq=""}
$0 ~ /^>/ {
if(seq!="") print seq
seq=""
next
}
{seq=seq $0}
END{print seq}
' spades.fasta > spades_consensus.fasta
grep -c ">" spades_consensus.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
grep -c ">" spades.fasta
grep -c ">" Genoma20_Erw2_grupo_externo.fasta
awk '
BEGIN{seq=""}
$0 ~ /^>/ {next}
{seq=seq $0}
END{
print ">Query_Spades"
print seq
}
' spades.fasta > spades_single.fasta
grep -c ">" spades_single.fasta
awk '/^>/ {print ">Ref_" ++i; next} {print}' Genoma20_Erw2_grupo_externo.fasta > ref_clean.fasta
grep -c ">" ref_clean.fasta
cat spades_single.fasta ref_clean.fasta > dataset_final.fasta
grep -c ">" dataset_final.fasta
FastTree -nt dataset_final.fasta > arbol.nwk
sudo apt update
sudo apt install fasttree -y
which FastTree
FastTree -help
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
ls -lh
FastTree -nt dataset_final.fasta > arbol.nwk
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
cp spades.fasta spades_clean.fasta
awk '/^>/ {print ">Ref_" ++i; next} {print}' Genoma20_Erw2_grupo_externo.fasta > ref_clean.fasta
cat spades_clean.fasta ref_clean.fasta > dataset_final.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
rm -f spades_vh5*.fasta
rm -f spades_single.fasta
rm -f spades_consensus.fasta
rm -f spades_mega.fasta
rm -f dataset_mega.fasta
rm -f dataset_final.fasta
rm -f ref_clean.fasta
rm -f aligned.fasta
rm -f arbol.nwk
rm -f resultado_vamphyre.txt
rm -f dist_salida*
rm -f vh_salida.txt
ls -lh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
rm -f spades_clean.fasta
rm -f ref_mega.fasta
rm -f spades_probe.fasta
rm -f vh_salida.txt
rm -f dist_salida*
rm -f arbol.nwk
grep -c ">" spades.fasta
chmod +x ContigsConcatenation.sh
./ContigsConcatenation.sh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
find . -name "ContigsConcatenation*"
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
chmod +x ContigsConcatenation.sh
cd originales
../ContigsConcatenation.sh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
ls -lh ContigsConcatenation.*
chmod +x ContigsConcatenation.sh
chmod +x ContigsConcatenation.pl
./ContigsConcatenation.sh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
chmod +x ContigsConcatenation.sh ContigsConcatenation.pl
./ContigsConcatenation.sh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
ls -lh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
ls -lh
find . -type f -mmin -10
ls -lh spades*
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
ls -lh
find . -type f -mmin -10
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
awk '
BEGIN{c=0}
/^>/ {c++; if(c>25) exit}
{print}
' spades.fasta > spades_probe.fasta
grep -c ">" spades_probe.fasta
../VAMPhyRE/bin/VH5cmdl -PROBEFILE spades_probe.fasta -TARGETLIST lista.txt -OUTFILE vh_salida.txt -MISMATCHES 1 -STRAND both -RESULTSFORMAT reduced
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
awk '
/^>/ {n++; file="genoma_" n ".fasta"}
{print > file}
' Genoma20_Erw2_grupo_externo.fasta
ls genoma_*.fasta | wc -l
ls genoma_*.fasta > lista.txt
../VAMPhyRE/bin/VH5cmdl -PROBEFILE spades_probe.fasta -TARGETLIST lista.txt -OUTFILE vh_salida.txt -MISMATCHES 1 -STRAND both -RESULTSFORMAT reduced
head -n 1 lista.txt > test_list.txt
../VAMPhyRE/bin/VH5cmdl -PROBEFILE spades_probe.fasta -TARGETLIST test_list.txt -OUTFILE test_vh.txt -MISMATCHES 1 -STRAND both -RESULTSFORMAT reduced
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/originales
rm -f lista.txt
rm -f vh_salida.txt
rm -f test_vh.txt
rm -f dist_salida*
rm -f arbol.nwk
rm -f spades_probe.fasta
rm -f genoma_*.fasta
rm -f spades_clean.fasta
rm -f ref_mega.fasta
rm -f dataset*.fasta
rm -f *.meg
rm -f *.msh
ls -lh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
conda activate vamphyre
mkdir -p 01_input 02_vh5 03_vfat 04_mega
mv originales/spades.fasta 01_input/
mv originales/Genoma20_Erw2_grupo_externo.fasta 01_input/
ls -lh 01_input
awk '
BEGIN{c=0}
/^>/ {c++; if(c>25) exit}
{print}
' 01_input/spades.fasta > 02_vh5/spades_probe.fasta
grep -c ">" 02_vh5/spades_probe.fasta
echo "/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/01_input/Genoma20_Erw2_grupo_externo.fasta" > 02_vh5/lista.txt
cat 02_vh5/lista.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin
chmod +x VH5cmdl
./VH5cmdl -PROBEFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/02_vh5/spades_probe.fasta -TARGETLIST /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/02_vh5/lista.txt -OUTFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/02_vh5/vh_salida.txt -MISMATCHES 1 -STRAND both -RESULTSFORMAT reduced
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/01_input
head -n 50 Genoma20_Erw2_grupo_externo.fasta > mini_target.fasta
echo "/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/01_input/mini_target.fasta" > ../02_vh5/lista_min.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/VAMPhyRE/bin
./VH5cmdl -PROBEFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/02_vh5/spades_probe.fasta -TARGETLIST /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/02_vh5/lista_min.txt -OUTFILE test_min.txt -MISMATCHES 1 -STRAND both -RESULTSFORMAT reduced
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
mkdir -p 03_vfat
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/01_input/*.fasta > 03_vfat/lista.txt
mkdir -p 03_vfat
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre/01_input/*.fasta > 03_vfat/lista.txt
cat 03_vfat/lista.txt
VFAT -VHFILE vh_salida.txt -TARGETLIST 03_vfat/lista.txt -OUTFILE 03_vfat/dist_salida -LEFTEXT 4 -RIGHTEXT 4 -THRESHOLD 19 -MODE DISTANCE
find /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre -name "VFAT*"
./ruta/VFAT ...
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
rm -rf 03_vfat 04_mega 02_vh5
rm -f *.txt
rm -f *.meg
rm -f *.nwk
rm -f *.out
rm -f *.log
rm -f spades_probe*
rm -f genoma_*.fasta
mkdir -p 01_input 02_vh5 03_vfat 04_mega
mv originales/spades.fasta 01_input/
mv originales/Genoma20_Erw2_grupo_externo.fasta 01_input/
find . -maxdepth 3 -type f
ls -lh originales
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/analisis_vamphyre
find . -type f -name "*.fasta"
ps -ef | grep VH5cmdl
top
ls -lh /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/vh_salida.txt
ps -ef | grep VH5cmdl
top
while true; do echo "===== $(date) ====="; ps -ef | grep VH5cmdl | grep -v grep; echo ""; ls -lh /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/vh_salida.txt; echo ""; sleep 60; done
while true; do echo "===== $(date) =====";  if ps -ef | grep VH5cmdl | grep -v grep > /dev/null; then     echo "VH5: SIGUE CORRIENDO"; else     echo "VH5: TERMINADO";     echo "LISTO VH5" | tee -a vh_monitor_log.txt;     break; fi;  ls -lh /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/vh_salida.txt;  echo ""; sleep 60; done
while true; do echo "===== $(date) =====";  if ps -ef | grep VH5cmdl | grep -v grep > /dev/null; then     echo "VH5: SIGUE CORRIENDO"; else     echo "VH5: TERMINADO";     echo "LISTO VH5" | tee -a vh_monitor_log.txt;     break; fi;  ls -lh /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/vh_salida.txt;  echo ""; sleep 60; done
top
ls -lh vh_salida.txt
ls -lh /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/vh_salida.txt
head -n 20 /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/vh_salida.txt
cat /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/lista.txt
grep -c ">" /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/Genoma20.join.fasta
grep -c ">" /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/probes.fasta
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
conda activate vamphyre
source activate vamphyre
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
ls -lh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
chmod +x ContigsConcatenation.sh ContigsConcatenation.pl
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales
ls -lh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
./ContigsConcatenation.sh
cd originales
ls -lh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales
ls -lh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
sed -n '1,200p' ContigsConcatenation.sh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales
mkdir -p backup_originales
cp *.fasta backup_originales/
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
nano ContigsConcatenation.sh
chmod +x ContigsConcatenation.sh
./ContigsConcatenation.sh
ls -lh ../originales
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto
find . -type d -name "originales"
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto
find . -name "*.join.fasta"
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales
ls -lh *.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
perl ContigsConcatenation.pl
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales
perl ../ContigsConcatenation.pl spades.fasta spades.join.fasta spades
perl ../ContigsConcatenation.pl Genoma20_Erw2_grupo_externo.fasta Genoma20.join.fasta Genoma20
ls -lh *.join.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
ls -lh VH5cmdl
chmod +x VH5cmdl
cd ../originales
nano lista.txt
cd ../originales
cp spades.join.fasta probes.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales
ls -lh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
chmod +x VH5cmdl
ls -lh VH5cmdl
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales
cp spades.join.fasta probes.fasta
nano lista.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales
rm -f lista.txt
nano lista.txt
cat lista.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
chmod +x VH5cmdl
ls -lh VH5cmdl
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/probes.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales
cp spades.join.fasta probes.fasta
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/probes.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
./VH5cmdl -PROBEFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/probes.fasta -TARGETLIST /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/lista.txt -OUTFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/vh_salida.txt -MISMATCHES 1 -STRAND both -RESULTSFORMAT reduced
mkdir -p /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/split_genomes
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales
awk '/^>/{f="genoma_"++i".fasta"}{print >> f}' Genoma20_Erw2_grupo_externo.fasta
ls genoma_*.fasta > lista.txt
ls -lh genoma_*.fasta
wc -l genoma_*.fasta
awk '
/^>/{close(out); out="genoma_"++i".fasta"}
{print >> out}
' Genoma20_Erw2_grupo_externo.fasta
ls -lh genoma_*.fasta
pwd
head lista.txt
wc -l lista.txt
ls genoma_*.fasta | sort -V > lista.txt
wc -l lista.txt
head lista.txt
tail lista.txt
realpath genoma_*.fasta | sort -V > lista.txt
wc -l lista.txt
tail lista.txt
head -n 1 lista.txt
ls -lh $(head -n 1 lista.txt)
seqkit --version
mkdir genomas_split
seqkit split2   -i Genoma20_Erw2_grupo_externo.fasta   -O genomas_split
sudo apt install seqkit
mkdir genomas_split
seqkit split2   -i Genoma20_Erw2_grupo_externo.fasta   -O genomas_split
sudo apt install seqkit
mkdir genomas_split
seqkit split2   -i Genoma20_Erw2_grupo_externo.fasta   -O genomas_split
mkdir genomas_split
seqkit split2   -i Genoma20_Erw2_grupo_externo.fasta   -O genomas_split
cd originales/
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales
unzip ncbi_dataset
unzip ncbi_dataset.zip
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto
ls originales
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/originales
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto
mkdir genomas_concat
chmod +x ContigsConcatenation.sh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto
mkdir genomas_concat
chmod +x ContigsConcatenation.sh
chmod +x /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/ContigsConcatenation.sh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
chmod +x ContigsConcatenation.sh
originales/
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
chmod +x /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/ContigsConcatenation.sh
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/ContigsConcatenation.sh /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/*.fasta > /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt
cat /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
chmod +x VH5cmdl
./VH5cmdl -PROBEFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/vps13/vps13.txt -TARGETLIST /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt -OUTFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
mkdir -p /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join
for f in /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/*.fasta; do     base=$(basename "$f" .fasta);      awk '
    /^>/ {header=$0; next}
    {
        gsub(/N+/, "") 
        if(length($0) > 200) print header"\n"$0
    }' "$f"     > /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join/${base}.join.fasta; done
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join/*.join.fasta > /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt
cat /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
rm -rf originales_join
rm -f list.txt
rm -f vh_salida.txt
rm -f dist_salida*
mkdir /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join
for f in /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/*.fasta; do     base=$(basename "$f" .fasta);      awk '
    /^>/ {h=$0; next}
    {
        gsub(/N+/, "")
        if(length($0) > 200) print h"\n"$0
    }' "$f"     > /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join/${base}.join.fasta; done
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join/*.join.fasta > /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt
wc -l /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join
ls -lh /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
rm -rf originales_join
rm -f list.txt
mkdir originales_join
for f in /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/*.{fna,fasta}; do     [ -e "$f" ] || continue;      base=$(basename "$f");     base="${base%.*}";      awk '
    /^>/ {h=$0; next}
    {
        gsub(/N+/, "")
        if(length($0) > 200) print h"\n"$0
    }' "$f"     > /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join/${base}.join.fasta; done
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join/*.join.fasta > /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt
wc -l /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt
sort /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt | uniq -d
while read f; do ls "$f"; done < /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt
wc -l list.txt
while read f; do ls "$f"; done < /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt
./VH5cmdl ...
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join | grep -i out
wc -l /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt
wc -l list.txt
→ 22
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
./VH5cmdl -PROBEFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/vps13/vps13.txt -TARGETLIST /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt -OUTFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
chmod +x ContigsConcatenation.sh
./ContigsConcatenation.sh /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join/*.join.fasta > /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
./VH5cmdl -PROBEFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/vps13/vps13.txt -TARGETLIST /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt -OUTFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/vh_salida.txt -MISMATCHES 1 -STRAND both
while read f; do     if [ ! -f "$f" ]; then         echo "❌ NO EXISTE: $f";     else         echo "OK: $f";     fi; done < /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt
ls -lh archivo.join.fasta
while read f; do     if [ ! -f "$f" ]; then         echo "❌ NO EXISTE: $f";     else         echo "OK: $f";     fi; done < /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt
ls -lh /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join/
ls -lh /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join/NOMBRE_DEL_ARCHIVO.join.fasta
for f in /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales/*.{fna,fasta}; do     [ -e "$f" ] || continue;      base=$(basename "$f");     base="${base%.*}";      awk '
    /^>/ {h=$0; print h; next}
    {
        seq=$0
        gsub(/N+/, "", seq)
        if(length(seq) > 200) print seq
    }' "$f"     > /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join/${base}.join.fasta; done
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join/*.join.fasta > /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
./VH5cmdl -PROBEFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/vps13/vps13.txt -TARGETLIST /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt -OUTFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
cat list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join
ls -lh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
./ContigsConcatenation.sh originales
cd originales_join
ls -lh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
./ContigsConcatenation.sh originales
cd originales_join
ls -lh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
rm -rf originales_join
mkdir originales_join
for f in originales/*.{fna,fasta}; do     [ -e "$f" ] || continue;      base=$(basename "$f");     base="${base%.*}";      awk '
    /^>/ {print $0; next}
    {
        gsub(/N+/, "", $0)
        if(length($0) > 200) print $0
    }' "$f"     > originales_join/${base}.join.fasta; done
ls originales_join/*.join.fasta > list.txt
ls -lh originales_join
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
rm -rf originales_join
mkdir originales_join
for f in originales/*.{fna,fasta}; do     [ -e "$f" ] || continue;      base=$(basename "$f");     base="${base%.*}";      awk '
    /^>/ {print $0; next}
    {
        seq=$0
        gsub(/N+/, "", seq)
        seq = seq seq
        if(length(seq) > 200) print seq
    }' "$f"     > originales_join/${base}.join.fasta; done
ls -lh originales_join
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
rm -rf originales_join
mkdir originales_join
cd originales
for f in *.fna; do     base=$(basename "$f" .fna);     cp "$f" "${base}.fasta"; done
cd ..
for f in originales/*.fasta; do     base=$(basename "$f" .fasta);      awk '
    /^>/ {h=$0; next}
    {
        seq=$0
        gsub(/N+/, "", seq)
        if(length(seq) > 200) print seq
    }' "$f"     > originales_join/${base}.join.fasta; done
ls originales_join/*.join.fasta > list.txt
cd originales_join
ls -lh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
rm -rf originales_join
mkdir originales_join
for f in originales/*.fasta; do     base=$(basename "$f" .fasta);      awk '
    /^>/ {print $0; next}
    {
        gsub(/[ \t]/, "")
        seq = seq $0
    }
    END {
        if(length(seq) > 0)
            print seq
    }' "$f"     > originales_join/${base}.join.fasta; done
ls originales_join/*.join.fasta > list.txt
ls -lh originales_join
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
ls originales_join/*.join.fasta > list.txt
cd VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../../vps13/vps13.txt -TARGETLIST ../../list.txt -OUTFILE ../../vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
cat -A list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
rm -f list.txt
ls -1 originales_join/*.join.fasta > list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
while read f; do     echo "Checking: $f";     ls -lh "$f"; done < list.txt
cd VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../../vps13/vps13.txt -TARGETLIST ../../list.txt -OUTFILE ../../vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
cat -A list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
ls -1 originales_join/*.join.fasta | sed 's|^|../|' > list.txt
ls -1 originales_join/*.join.fasta > list.txt
while read f; do     echo "Checking: $f";     test -f "$f" && echo "OK" || echo "FAIL"; done < list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
./VAMPhyRE/bin/VH5cmdl -PROBEFILE vps13/vps13.txt -TARGETLIST list.txt -OUTFILE vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join
for f in *.join.fasta; do     echo "Checking $f";     awk 'NR==1{ok=1} END{if(ok) print "OK"}' "$f" || echo "PROBLEMA: $f"; done
grep -L "^>" *.join.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join/*.join.fasta > /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt
./VH5cmdl -PROBEFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/vps13/vps13.txt -TARGETLIST /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt -OUTFILE /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
cp /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join/*.join.fasta .
ls *.join.fasta > list.txt
cp /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/vps13/vps13.txt .
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
find . -name "vps13.txt"
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
cp ../../vps13.txt .
ls *.join.fasta > list.txt
./VH5cmdl -PROBEFILE vps13.txt -TARGETLIST list.txt -OUTFILE vh_salida.txt -MISMATCHES 1 -STRAND both
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/vh_salida.txt
./VFAT -VHFILE vh_salida.txt -TARGETLIST list.txt -OUTFILE dist_salida -LEFTEXT 4 -RIGHTEXT 4 -THRESHOLD 20 -MODE DISTANCE
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
ls -l VFAT
chmod +x VFAT
ls -l VFAT
./VFAT -VHFILE vh_salida.txt -TARGETLIST list.txt -OUTFILE dist_salida -LEFTEXT 4 -RIGHTEXT 4 -THRESHOLD 20 -MODE DISTANCE
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
ls -lh dist_salida*
/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin/dist_salida.meg
file dist_salida
head dist_salida
ls -lh dist_salida*
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
head -n 30 dist_salida.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
nano dist_salida.meg
/dist_salida.meg
cd /dist_salida.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
head -n 40 dist_salida.meg
sed -n '25,60p' dist_salida.meg
-MODE DISTANCE -OUTFORMAT PHYLIP
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
ls -lh dist_salida.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
sed -E 's/^\[[0-9]+\]\s+#/>/' dist_salida.meg > dist_salida_FIXED.meg
head dist_salida_FIXED.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
awk '
BEGIN {print "22 22"}
/^\[/ {
  gsub(/^\[[0-9]+\] #/, "", $0)
  split($0,a," ")
  name=a[1]
  $1=""; $2=""
  print name, $0
}' dist_salida.meg > dist_salida.phy
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
sed -E 's/^\[[0-9]+\]\s+#/>/' dist_salida.meg > dist_salida_FIXED.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
cp dist_salida.meg dist_salida_backup.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
awk '
BEGIN{flag=0}
{
    if($0 ~ /^>/){
        gsub(/^>/,"[ ")
        print
        next
    }
    print
}' dist_salida.meg > dist_salida_fixed.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
awk '
/^\[/{
    if($0 !~ /#/){
        sub(/\] /, "] #")
    }
}
{print}
' dist_salida_fixed.meg > dist_salida_clean.meg
grep -c "^\[" dist_salida_clean.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
wc -l list.txt
cat list.txt
22 target sequences read...
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
rm -f dist_salida_fixed.meg dist_salida_clean.meg
rm -f dist_salida.meg
./VFAT -VHFILE vh_salida.txt -TARGETLIST /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt -OUTFILE dist_salida -LEFTEXT 4 -RIGHTEXT 4 -THRESHOLD 19 -MODE DISTANCE
grep -c "\[" dist_salida.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
grep -c "^\[" dist_salida.meg
wc -l /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt
sort list.txt | uniq -c | sort -nr | head
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
rm -f list.txt
ls originales_join/*.join.fasta | sort > list.txt
wc -l list.txt
cd VAMPhyRE/bin
rm -f vh_salida.txt dist_salida.*
./VH5cmdl -PROBEFILE vps13.txt -TARGETLIST ../list.txt -OUTFILE vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
while read f; do     if [ ! -f "$f" ]; then         echo "❌ FALTA: $f";     else         echo "OK: $f";     fi; done < list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join
for f in *.join.fasta; do     echo "Revisando $f";     head -n 1 "$f"; done
for f in *.join.fasta; do     if ! grep -q "^>" "$f"; then         echo "❌ SIN HEADER FASTA: $f";     fi; done
cd originales_join
dos2unix *.join.fasta 2>/dev/null
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../../list.txt -OUTFILE vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
ls -1 originales_join/*.join.fasta > list.txt
realpath originales_join/*.join.fasta > list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
dos2unix list.txt
dos2unix originales_join/*.join.fasta
while read f; do ls -lh ../"$f"; done < ../list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
ls -lh list.txt
ls -1 originales_join/*.join.fasta > list.txt
while read f; do ls -lh "$f"; done < list.txt
while read f; do     if [ -f "$f" ]; then         echo "OK: $f";     else         echo "FALTA: $f";     fi; done < list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
while read f; do     echo "---- $f ----";     head -n 2 "$f"; done < list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join
mkdir solo_chromosome
for f in *.join.fasta; do     awk '/^>/ {keep=($0 ~ /chromosome/)} keep' "$f" > solo_chromosome/"$f"; done
grep -c ">" solo_chromosome/*.fasta
grep -L "^>" solo_chromosome/*.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
rm solo_chromosome/spades.join.fasta
realpath solo_chromosome/*.fasta > list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
ls -lh
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
ls originales_join/*.join.fasta > list.txt
cat list.txt
cd VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
ls -lh list.txt
cat list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
realpath originales_join/*.join.fasta > list.txt
head list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
while read f; do     ls -lh "$f" || echo "❌ NO EXISTE: $f"; done < ../list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
while read f; do     ls -lh "/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/$f"     || echo "❌ NO EXISTE: $f"; done < /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
realpath originales_join/*.join.fasta > list.txt
head list.txt
cd VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
while read f; do     echo "CHECK -> $f";     if [ ! -e "$f" ]; then         echo "❌ NO EXISTE";         break;     fi; done < list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
find $(pwd)/originales_join -name "*.join.fasta" | sort > list.txt
while read f; do     test -s "$f" || echo "❌ PROBLEMA: $f"; done < list.txt
ls -lh /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join/solo_chromosome/
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
rm list.txt
find originales_join -type f -name "*.join.fasta" -size +0c | sort > list.txt
while read f; do     test -s "$f" || echo "❌ PROBLEMA: $f"; done < list.txt
cd VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
cd VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
while read f; do     echo "ABRIENDO -> $f";     head -n 1 "$f" >/dev/null || {         echo "💥 FALLA AQUÍ: $f";         break;     }; done < list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
rm -f list.txt
find originales_join -type f -name "*.join.fasta" -size +1k | sort > list.txt
while read f; do     echo "OK -> $f";     test -s "$f" || echo "💥 MAL: $f"; done < list.txt
cd VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
grep -v "NODE_" list.txt > list_clean.txt
chmod -R 755 originales_join/
dos2unix originales_join/*.fasta 2>/dev/null
chmod -R 755 /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join/
dos2unix /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join/*.fasta 2>/dev/null
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST list_clean.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
while read f; do     echo "OK: $f";     head -n 1 "$f" || echo "💥 FALLA: $f"; done < list_clean.txt
grep -v "spades" list_clean.txt > list_no_spades.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST list_no_spades.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
realpath originales_join/*.join.fasta > list.txt
head list.txt
realpath originales_join/spades.join.fasta >> list.txt
sort list.txt -o list.txt
cd VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
mkdir -p solo_chromosome
grep -L "plasmid\|NODE" originales_join/*.join.fasta
for f in originales_join/*.join.fasta; do     if grep -q "chromosome" "$f"; then         cp "$f" solo_chromosome/;     fi; done
realpath solo_chromosome/*.join.fasta > list.txt
realpath originales_join/spades.join.fasta >> list.txt
cat list.txt
while read f; do     if [ ! -s "$f" ]; then         echo "❌ PROBLEMA: $f";     else         echo "OK: $f";     fi; done < list.txt
cd VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
while read f; do     echo "PROBANDO: $f";     head -n 2 "$f" > /dev/null || echo "❌ FALLA: $f"; done < list.txt
while read f; do     if [ ! -s "$f" ]; then         echo "❌ VACÍO O CORRUPTO: $f";     elif ! head -n 1 "$f" | grep -q "^>"; then         echo "❌ NO FASTA: $f";     else         echo "OK: $f";     fi; done < list.txt
mkdir -p solo_chromosome_clean
for f in solo_chromosome/*.fasta; do     awk '/^>/{if(n++) exit} {print}' "$f" > solo_chromosome_clean/$(basename "$f"); done
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/solo_chromosome_clean/*.fasta > list.txt
cd VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
head -n 5 list.txt
ls GCF_000026185.1_ASM2618v1_genomic.join.fasta
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
realpath originales_join/*.join.fasta > list.txt
realpath originales_join/spades.join.fasta >> list.txt
head list.txt
cd VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
cd Genomas_ref
mkdir -p vh_input
cp originales_join/*.join.fasta vh_input/
ls vh_input/*.fasta > list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
mkdir -p vh_input
cp originales_join/*.join.fasta vh_input/
ls vh_input | head
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/vh_input/*.fasta > list.txt
head list.txt
while read f; do     if [ ! -s "$f" ]; then         echo "ERROR: $f";     else         echo "OK: $f";     fi; done < list.txt
cd VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
pwd
ls vh_input | head
head list.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
ls
mkdir -p vh_input
cp originales_join/*.join.fasta vh_input/
ls /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/vh_input/*.fasta > list.txt
head list.txt
while read f; do     test -s "$f" && echo "OK: $f" || echo "ERROR: $f"; done < list.txt
cd VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
#!/bin/bash
set -e
BASE="/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref"
echo "📁 Entrando a proyecto..."
cd "$BASE"
# -----------------------------
# 1. Crear carpeta de trabajo limpia
# -----------------------------
echo "🧹 Creando vh_input limpio..."
rm -rf vh_input
mkdir -p vh_input
# -----------------------------
# 2. Copiar genomas
# -----------------------------
echo "📦 Copiando genomas..."
cp originales_join/*.join.fasta vh_input/
# -----------------------------
# 3. Crear lista absoluta
# -----------------------------
echo "🧾 Generando list.txt..."
ls "$BASE"/vh_input/*.fasta > list.txt
# -----------------------------
# 4. Validación básica
# -----------------------------
echo "🔍 Validando archivos..."
while read f; do     if [ ! -s "$f" ]; then         echo "❌ ERROR archivo vacío: $f";         exit 1;     fi; done < list.txt
echo "✅ Todos los archivos OK"
# -----------------------------
# 5. Ir a VAMPhyRE
# -----------------------------
echo "🧬 Ejecutando VH5cmdl..."
cd "$BASE/VAMPhyRE/bin"
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
#!/bin/bash
set -e
BASE="/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref"
echo "📁 Entrando a proyecto..."
cd "$BASE"
# -----------------------------
# 1. Crear carpeta de trabajo limpia
# -----------------------------
echo "🧹 Creando vh_input limpio..."
rm -rf vh_input
mkdir -p vh_input
# -----------------------------
# 2. Copiar genomas
# -----------------------------
echo "📦 Copiando genomas..."
cp originales_join/*.join.fasta vh_input/
# -----------------------------
# 3. Crear lista absoluta
# -----------------------------
echo "🧾 Generando list.txt..."
ls "$BASE"/vh_input/*.fasta > list.txt
# -----------------------------
# 4. Validación básica
# -----------------------------
echo "🔍 Validando archivos..."
while read f; do     if [ ! -s "$f" ]; then         echo "❌ ERROR archivo vacío: $f";         exit 1;     fi; done < list.txt
echo "✅ Todos los archivos OK"
# -----------------------------
# 5. Ir a VAMPhyRE
# -----------------------------
echo "🧬 Ejecutando VH5cmdl..."
cd "$BASE/VAMPhyRE/bin"
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
#!/bin/bash
BASE="/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref"
echo "📁 Entrando a proyecto..."
cd "$BASE" || { echo "❌ No se pudo entrar a $BASE"; exit 1; }
echo "🧹 Limpiando input..."
rm -rf vh_input
mkdir -p vh_input
echo "📦 Copiando genomas..."
cp originales_join/*.join.fasta vh_input/ 2>/dev/null
echo "🧾 Generando list.txt..."
ls "$BASE"/vh_input/*.fasta > list.txt 2>/dev/null
echo "🔍 Validando archivos..."
while read f; do     if [ ! -s "$f" ]; then         echo "⚠️ Archivo vacío o inválido: $f";     else         echo "OK: $f";     fi; done < list.txt
echo "🧬 Ejecutando VH5cmdl..."
cd "$BASE/VAMPhyRE/bin" || exit 1
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
VH_STATUS=$?
echo "VH5cmdl exit code: $VH_STATUS"
if [ $VH_STATUS -ne 0 ]; then     echo "❌ VH5cmdl falló, pero el script continúa"; fi
chmod +x VFAT
PROBE_SIZE=13
LEFT=4
RIGHT=4
MM=1
THRESHOLD=$((PROBE_SIZE + LEFT + RIGHT - MM))
echo "📊 THRESHOLD = $THRESHOLD"
./VFAT -VHFILE ../vh_salida.txt -TARGETLIST ../list.txt -OUTFILE ../dist_salida -LEFTEXT $LEFT -RIGHTEXT $RIGHT -THRESHOLD $THRESHOLD -MODE DISTANCE
VFAT_STATUS=$?
echo "VFAT exit code: $VFAT_STATUS"
echo "🧬 Ajustando MEGA..."
cd "$BASE"
if [ -f dist_salida.meg ]; then     sed -i 's/,.*$//' dist_salida.meg;     dos2unix dist_salida.meg 2>/dev/null;     echo "✅ MEGA listo"; else     echo "⚠️ No se generó dist_salida.meg"; fi
echo "🎉 FIN DEL PIPELINE (sin cierre de terminal)"
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
mkdir -p vh_input_clean
for f in originales_join/*.join.fasta; do     awk '/^>/ {print; next} {gsub(/[^ATCGN]/,""); print}' "$f" > vh_input_clean/$(basename "$f"); done
cd ~/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join
awk '
/^>/ {next}
{printf "%s", $0}
END {print ""}
' spades.join.fasta > spades.one_line.fasta
echo ">spades_assembly" > spades_clean.fasta
cat spades.one_line.fasta >> spades_clean.fasta
awk '
/^>/ {if(seq){print seq}; seq=""; next}
{seq=seq $0}
END {print seq}
' spades.join.fasta > tmp.seq
echo ">spades_assembly" > spades_clean.fasta
cat tmp.seq >> spades_clean.fasta
cd ~/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
mkdir -p vh_input_final
cp originales_join/*.join.fasta vh_input_final/
cp originales_join/spades_clean.fasta vh_input_final/
ls vh_input_final/*.fasta > list_final.txt
grep ">" vh_input_final/spades_clean.fasta
cd VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list_final.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
cat ../list_final.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
ls vh_input_final
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
realpath originales_join/*.join.fasta > list_final.txt
realpath originales_join/spades_clean.fasta >> list_final.txt
cat list_final.txt
sort list_final.txt | uniq > list_final.tmp && mv list_final.tmp list_final.txt
cd VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list_final.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
while read f; do     if [ -f "$f" ]; then         echo "OK: $f";     else         echo "❌ NO EXISTE: $f";     fi; done < list_final.txt
while read f; do     if [ -f "$f" ]; then         echo "OK: $f";     else         echo "❌ NO EXISTE: $f";     fi; done < ../list_final.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
ls -lh list_final.txt
cd Genomas_ref
ls -lh list_final.txt
while read f; do     if [ ! -f "$f" ]; then         echo "❌ NO EXISTE: $f";     fi; done < list_final.txt
while read f; do     if [ ! -f "$f" ]; then         echo "❌ NO EXISTE: $f";     fi; done < list_final.txt
ls -lh ../vps13/vps13.txt
head -n 5 ../vps13/vps13.txt
dos2unix ../vps13/vps13.txt
dos2unix list_final.txt
head -n 3 list_final.txt
tail -n 3 list_final.txt
grep "GCF_" list_final.txt > list_vh5_clean.txt
wc -l list_vh5_clean.txt
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list_vh5_clean.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../../vps13/vps13.txt -TARGETLIST ../../list_vh5_clean.txt -OUTFILE ../../vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
ls -lh VH5cmdl
touch test_write.txt
strace -o vh_debug.txt ./VH5cmdl -PROBEFILE ../../vps13/vps13.txt -TARGETLIST ../../list_vh5_clean.txt -OUTFILE ../../vh_salida.txt -MISMATCHES 1 -STRAND both
sudo apt update
sudo apt install strace
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/VAMPhyRE/bin
strace -o vh_debug.txt ./VH5cmdl -PROBEFILE ../../vps13/vps13.txt -TARGETLIST ../../list_vh5_clean.txt -OUTFILE ../../vh_salida.txt -MISMATCHES 1 -STRAND both
grep "ENOENT" vh_debug.txt | tail
./VH5cmdl -PROBEFILE ../../vps13.txt -TARGETLIST ../../list_vh5_clean.txt -OUTFILE ../../vh_salida.txt -MISMATCHES 1 -STRAND both
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
ls -lh originales_join/spades_clean.fasta
grep "GCF_" list_final.txt > list_tmp.txt
echo "/home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref/originales_join/spades_clean.fasta" >> list_tmp.txt
while read f; do     if [ ! -f "$f" ]; then         echo "❌ NO EXISTE: $f";     else         echo "OK: $f";     fi; done < list_tmp.txt
mv list_tmp.txt list_final.txt
wc -l list_final.txt
cd VAMPhyRE/bin
./VH5cmdl -PROBEFILE ../../vps13.txt -TARGETLIST ../../list_final.txt -OUTFILE ../../vh_salida.txt -MISMATCHES 1 -STRAND both
./VFAT -VHFILE ../../vh_salida.txt -TARGETLIST ../../list_final.txt -OUTFILE ../../dist_salida -LEFTEXT 4 -RIGHTEXT 4 -THRESHOLD 20 -MODE DISTANCE
head -n 30 dist_salida.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
ls -lh dist_salida*
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
ls -lh dist_salida*
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
head -n 20 dist_salida.meg
head -n 30 dist_salida.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
sed -i 's/#spades_assembly/SPADES/' dist_salida.meg
tail -n 5 dist_salida.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
head -n 40 dist_salida.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
sed -i 's/#spades_assembly/SPADES/g' dist_salida.meg
sed -i 's/#//g' dist_salida.meg
head -n 30 dist_salida.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
sed -i '1d' dist_salida.meg
head -n 5 dist_salida.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
awk 'BEGIN{print "#MEGA"} {print}' dist_salida.meg > dist_salida_meg_fix.meg
grep "\[" dist_salida_meg_fix.meg | head
grep -c "\[" dist_salida_meg_fix.meg
grep "\[" dist_salida.meg | head -n 5
[ 1] NC_010694.1 Erwinia tasmaniensis Et1/99, complete sequence
[ 2] NZ_CP024970.1 Erwinia amylovora strain E-2 chromosome, complete genome
[ 3] NZ_CP050258.1 Erwinia amylovora strain FB-86 chromosome, complete genome
[ 4] NZ_CP050244.1 Erwinia amylovora strain TS3238 chromosome, complete genome
[ 5] NZ_CP050263.1 Erwinia amylovora strain FB-207 chromosome, complete genome
tail -n 5 dist_salida.meg
tail -n 20 dist_salida.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
sed '/^!/d' dist_salida.meg > mega_clean.meg
echo "#MEGA" | cat - mega_clean.meg > final_mega.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
awk '
BEGIN {print "22"} 
/^\\[/ { 
    printf "TAX%d ", NR
    for(i=2;i<=NF;i++) printf "%s ", $i
    printf "\n"
}
' dist_salida.meg > dist_salida.phy
head dist_salida.phy
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
awk '
BEGIN {n=0}
# detecta filas de matriz
/^\[[0-9]+\]/ {
    n++
    label[n]=$2
    row=""
    for(i=3;i<=NF;i++) row=row" "$i
    data[n]=row
}
END {
    print n
    for(i=1;i<=n;i++){
        printf "TAX%d %s\n", i, data[i]
    }
}
' dist_salida.meg > dist_salida.phy
head dist_salida.phy
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
grep -E "^\[[0-9]+\]" dist_salida.meg > matriz_raw.txt
wc -l matriz_raw.txt
grep "^\[" dist_salida.meg
grep -c "^\[" dist_salida.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
grep -E "^\[[0-9]+\]" dist_salida.meg | head -n 22 > matrix_only.txt
wc -l matrix_only.txt
awk '
{
    row=""
    for(i=3;i<=NF;i++) row=row" "$i
    print "Taxa_"NR row
}
' matrix_only.txt > dist_clean.phy
sed -i '1i 22' dist_clean.phy
head dist_clean.phy
wc -l matrix_only.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
awk '
{
    printf "Taxa_%02d ", NR
    for(i=1;i<=NF;i++) printf "%s ", $i
    printf "\n"
}
' matrix_only.txt > dist_clean.phy
sed -i '1i 22' dist_clean.phy
head dist_clean.phy
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
awk '
BEGIN {print "22"}
{
    name = sprintf("T%02d", NR)
    printf "%-10s", name
    for(i=1;i<=NF;i++) printf " %s", $i
    printf "\n"
}
' matrix_only.txt > dist_clean.phy
head dist_clean.phy
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
awk '
BEGIN {print "22"}
{
    printf "T%02d", NR
    for(i=1;i<=NF;i++) {
        if($i ~ /^[0-9.]+$/ || $i ~ /^-?[0-9.]+$/)
            printf " %s", $i
    }
    printf "\n"
}
' matrix_only.txt > dist_clean.phy
head dist_clean.phy
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
head -n 5 dist_salida.meg
grep -E "^\[[0-9]+\]" dist_salida.meg > raw_matrix.txt
awk '
BEGIN {print "22"}
{
    printf "T%02d", NR
    for(i=2;i<=NF;i++) {
        if($i ~ /^[0-9.]+$/ || $i ~ /^-?[0-9.]+$/)
            printf " %s", $i
    }
    printf "\n"
}
' raw_matrix.txt > dist_clean.phy
head dist_clean.phy
head -n 30 raw_matrix.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
grep -E "^\[[0-9]+\]" dist_salida.meg | grep -v "SPADES" > matrix_only.txt
head matrix_only.txt
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
grep -E "^\[[0-9]+\][[:space:]]+[0-9]" dist_salida.meg > matrix_only.txt
head matrix_only.txt
awk '
BEGIN {print "22"}
{
    printf "T%02d", NR
    for(i=2;i<=NF;i++)
        printf " %s", $i
    printf "\n"
}
' matrix_only.txt > dist_clean.phy
head dist_clean.phy
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
cp dist_salida.meg dist_final.meg
nano dist_final.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
ls -lh dist_final.meg
cd /home/robertolede/bioinformatica_2026-1/Practica_7_proyecto/Genomas_ref
nano dist_final.meg
history | grep -E "VH5cmdl|mega|nano|cp|ls|cat|grep|awk|sed|cd" > pipeline_bioinfo.txt
cat << 'EOF' > pipeline_final.sh
# 1. Preparación de genomas
cp originales_join/*.fasta vh_input/

# 2. Crear lista
ls vh_input/*.fasta > list_final.txt

# 3. Virtual Hybridization (VAMPhyRE)
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list_final.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both

# 4. (salida MEGA ya generada)
EOF

tar -czvf backup_bioinfo.tar.gz dist_final.meg list_final.txt vh_salida.txt vh_input/
