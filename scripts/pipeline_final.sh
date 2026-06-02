# 1. Preparación de genomas
cp originales_join/*.fasta vh_input/

# 2. Crear lista
ls vh_input/*.fasta > list_final.txt

# 3. Virtual Hybridization (VAMPhyRE)
./VH5cmdl -PROBEFILE ../vps13/vps13.txt -TARGETLIST ../list_final.txt -OUTFILE ../vh_salida.txt -MISMATCHES 1 -STRAND both

# 4. (salida MEGA ya generada)
