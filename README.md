# scEPCOT

## Usage

### Download source code

```
git clone https://github.com/zzh24zzh/scEPCOT.git
```

### Pre-process input data
Prepare the input ATAC-seq files
```
wget -O - https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz| gunzip -c > black_list.bed

bamCoverage --bam GM12878.bam -o GM12878_atac.bigWig --outFileFormat bigwig --normalizeUsing RPGC --effectiveGenomeSize 2913022398 
--Offset 1 --binSize 1 --numberOfProcessors 12 --blackListFileName /scratch/drjieliu_root/drjieliu/zhenhaoz/ATAC-seq/bam/black_list.bed

#usage: python process_atac.py -i <bigWig_file> -p <prefix> -o <output_directory>
python process_atac.py -i GM12878_atac.bigWig -p GM12878 -o atacseq/
```

Download models, input reference genome sequence and example data
```
### the files are downloaded under 'models', 'refSeq', and 'examples' folder
python download.py
```


### Run demo

```
python gradio_epcot.py
```

### Process prediction files
The files of prediction results will appear under the folder 'tmps' after running the demo.

Generate the Excel file and other supported files for [nuclesome browser](https://github.com/nucleome/nucleserver) from the prediction result files
```
#usage: python makesheet,py -i <input prediction files for browser> -o <out file name.xlsx>
python makesheet.py -i tmps/browser_data/tmp_11-10550000-10950000.zip -o tmps/tmp.xlsx
```
