# EPCOT Gradio demo

## Dependencies

* deepTools-3.5.1

* samtools-1.16.1

## Usage

### Step 1: Clone the repository

```
git clone https://github.com/zzh24zzh/EPCOT_gradio.git
```

### Step 2: Download trained models and reference genome data (download.py)

The pre-trained models will be downloaded to the "models/" directory, and the reference sequence data will be downloaded to the "refSeq/" directory. 
```
python download.py
```


### Step 3: Prepare input ATAC-seq data (process_atac.py)
#### Input
* an ATAC-seq profile in <span style="color:red"> .bam </span> format
* the number of processors used in bamConverage function

#### Output
* processed ATAC-seq file in .npz format


Required packages: deepTools, samtools
```
python process_atac.py -b <ATAC-seq bam file> -p <number of processors>
```



### Step 4: Run the demo

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
