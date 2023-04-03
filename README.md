# EPCOT Gradio

## Dependencies

deepTools-3.5.1

samtools-1.16.1

## Usage

### Step 1: Download source code

```
git clone https://github.com/zzh24zzh/EPCOT_gradio.git
```

### Step 2: Download trained models and reference genome data

```
python download.py
```

### Step 3: Prepare input ATAC-seq data 
Required packages: deepTools, samtools
```
python process_atac.py -b <ATAC-seq bam file> -p <number of processors>
```

### Download EPCOT models and input reference genome sequence data
```
###the files are downloaded under 'models' and 'refSeq' folders. Add '-d' if not to download the reference sequnece 
python download.py
```


### Run the demo

```
python gradio_epcot.py
```

### Process prediction files
The files of prediction results will appear under the folder 'tmps' after running the demo.

Generate the Excel file and other supported files for [nuclesome browser](https://github.com/nucleome/nucleserver) from the prediction result files
```
###usage: python makesheet,py -i <inputFilename> -o <outputFilename.xlsx>
python makesheet.py -i tmps/browser_data/tmp_11-10550000-10950000.zip -o tmps/nb.xlsx
```
