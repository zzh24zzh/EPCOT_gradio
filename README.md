# scEPCOT

## Usage

### pre-processing input data
Prepare the input ATAC-seq files
```
wget -O - https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz| gunzip -c > black_list.bed

bamCoverage --bam GM12878.bam -o GM12878_atac.bigWig --outFileFormat bigwig --normalizeUsing RPGC --effectiveGenomeSize 2913022398 
--Offset 1 --binSize 1 --numberOfProcessors 12 --blackListFileName /scratch/drjieliu_root/drjieliu/zhenhaoz/ATAC-seq/bam/black_list.bed

#usage: python process_atac.py -i <bigWig_file> -p <prefix> -o <output_directory>
python process_atac.py -i GM12878_atac.bigWig -p GM12878 -o atacseq/
```


Download and extract the processed input reference sequence
```
!mkdir refSeq
!gdown 1iqOKKNFwjl9hMZovxxhG-t_Y1MZw45R0 --output refSeq/hg38.tar.gz
!tar -xvf refSeq/hg38.tar.gz -C refSeq/
```

Generate the Excel file and other supported files for [nuclesome browser](https://github.com/nucleome/nucleserver) from the prediction result files
```
#usage: python makesheet,py -i <input prediction files for browser> -o <out file name.xlsx>
python makesheet.py -i tmps/browser_data/tmp_11-10550000-10950000.zip -o tmps/tmp.xlsx
```
