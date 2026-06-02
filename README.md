# Confirmación taxonómica y análisis genómico de Erwinia amylovora

## Descripción

Este repositorio contiene el análisis bioinformático realizado para la identificación taxonómica y caracterización genómica de un aislamiento bacteriano mediante secuenciación de nueva generación (NGS).

El flujo de trabajo incluyó evaluación de calidad de lecturas, clasificación taxonómica, ensamblaje de novo, evaluación de ensamblajes, anotación funcional y análisis filogenómico.

## Objetivo general

Confirmar la identidad taxonómica del microorganismo analizado e identificar genes asociados con virulencia, adaptación y resistencia.

## Herramientas utilizadas

- FastQC
- One Codex
- SPAdes
- QUAST
- Prokka
- VAMPhyRE
- MEGA
- iToL

## Resultados principales

### Calidad de lecturas

- 9,202,195 lecturas paired-end
- Longitud de lectura: 151 pb
- Contenido GC: 53 %
- Calidad Phred alta
- Sin secuencias de baja calidad

### Clasificación taxonómica

- Erwinia amylovora: 97.02 %

### Ensamblaje seleccionado (SPAdes)

- Contigs: 41
- N50: 378,097 pb
- Longitud total: 3,783,407 pb
- GC: 53.54 %

### Genes identificados mediante Prokka

- hrpB
- hrpN
- hrcN
- amsI
- rpoS
- rpsL
- strB

## Estructura del proyecto

- data/
- results/
- figures/
- scripts/
- docs/
- report/

## Autor

Roberto Ledezma Ruiz
Universidad Autónoma de Chihuahua
