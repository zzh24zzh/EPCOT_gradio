# EPCOT Gradio demo


Mapping ATAC-seq to epigenome, transcriptome, and high-resolution chromatin contact maps with an easy-to-use Gradio interface.

## Prerequisites

* Python 3.9.12


## Usage

Follow these steps to run the demo locally:

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
* An ATAC-seq profile in **.bam** format
* The number of processors to use in deepTools's bamCoverage

#### Output
* A processed ATAC-seq file in **.npz** format, which can be uploaded to the demo to excute models


#### Required packages: 
* deepTools-3.5.1
* samtools-1.16.1
* pyBigWig-0.3.17
```
python process_atac.py -b <ATAC-seq bam file> -p <number of processors>
```



### Step 4: Run Gradio demo

```
python gradio_epcot.py
```

