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

Para Árbol filogenómico: # Entrar al directorio de trabajo
cd ~/bioinformatica_2026-1/Practica_7_proyecto

# Descomprimir datos de NCBI
unzip ncbi_dataset.zip

# Localizar archivos fasta
find ncbi_dataset -name "*.fna"

# Renombrar o copiar el multifasta
cp archivo_original.fna Genoma20_Erw2_grupo_externo.fasta

# Contar secuencias del multifasta
grep -c "^>" Genoma20_Erw2_grupo_externo.fasta

# Crear carpeta para genomas individuales
mkdir originales

# Separar multifasta en archivos independientes
awk '
/^>/{
 if(out) close(out);
 c++;
 out=sprintf("genoma_%02d.fasta",c)
}
{
 print > out
}' Genoma20_Erw2_grupo_externo.fasta

# Mover archivos generados
mv genoma_*.fasta originales/

# Verificar cantidad
ls originales | wc -l

# Crear lista de genomas
find $(pwd)/originales -name "*.fasta" | sort > lista.txt

# Verificar lista
wc -l lista.txt

# Revisar contenido
cat lista.txt

# Dar permisos a VH
chmod +x VH

# Ejecutar VH
./VH lista.txt probes.txt salida.vh

# Revisar salida
less salida.vh

# Ejecutar VFAT
java -jar VFAT.jar salida.vh

# Revisar archivos generados
ls

# Convertir para MEGA
# (depende de la versión de VFAT)

# Abrir archivo en MEGA
megacc -a infer_NJ.mao archivo.meg
