# scEPCOT

## Usage

```
wget -O - https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz| gunzip -c > black_list.bed
bamCoverage --bam %s.bam -o ../normalize_atac/tmp_%s_atac.bigWig --outFileFormat bigwig --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --Offset 1 --binSize 1 --numberOfProcessors 12 --blackListFileName /scratch/drjieliu_root/drjieliu/zhenhaoz/ATAC-seq/bam/black_list.bed


```
