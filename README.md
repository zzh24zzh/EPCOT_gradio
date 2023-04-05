# EPCOT Gradio demo

Mapping ATAC-seq to epigenome, transcriptome, and high-resolution chromatin contact maps with an easy-to-use Gradio interface.

The demo includes two interfaces: (1) to run the model, and (2) to visualize prediction results. Users only need to provide a processed ATAC-seq file, which can be generated from a .bam file using the provided script [atac_process.py](https://github.com/zzh24zzh/EPCOT_gradio/blob/main/atac_process.py) (refer to step 3 in Usage).

## Prerequisites

* Python 3.9.12


## Usage

The demo can be deployed on Google Colab notebook where the GPUs can speed up running the models. Here is an example

Follow these steps to run the demo locally:

####  1. Clone the repository

```
git clone https://github.com/zzh24zzh/EPCOT_gradio.git
```

####  2. Download trained models and reference genome data (download.py)

The pre-trained models will be downloaded to the "models/" directory, and the reference sequence data will be downloaded to the "refSeq/" directory. 
```
python download.py
```


####  3. Prepare input ATAC-seq data (process_atac.py)
##### Input
* An ATAC-seq profile in **.bam** format
* The number of processors to use in deepTools's bamCoverage

##### Output
* A processed ATAC-seq file in **.npz** format, which can be uploaded to the demo to excute models


##### Required packages: 
* deepTools-3.5.1
* samtools-1.16.1
* pyBigWig-0.3.17
```
python process_atac.py -b <ATAC-seq bam file> -p <number of processors>
```



####  4. Launch the Gradio demo

```
python gradio_epcot.py
```

After completing these steps, you will be able to use the Gradio interface to run the EPCOT model on your ATAC-seq data and visualize the results.
