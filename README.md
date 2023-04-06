# EPCOT Gradio demo

This repository contains a Gradio demo for mapping ATAC-seq data to epigenomes, transcriptomes, and high-resolution chromatin contact maps.

The demo features two user-friendly interfaces: (1) to run the model, and (2) to visualize prediction results. Users only need to provide a processed ATAC-seq file, which can be generated from a .bam file using the provided script [atac_process.py](https://github.com/zzh24zzh/EPCOT_gradio/blob/main/atac_process.py) (refer to step 3 in Usage below).

![](https://github.com/zzh24zzh/EPCOT_gradio/blob/main/data/interface1.gif)


![](https://github.com/zzh24zzh/EPCOT_gradio/blob/main/data/interface2.gif)



## Prerequisites

* Python 3.9.12


## Usage

The demo can be deployed on a Google Colab notebook where GPUs can accelerate the model execution. Click [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/zzh24zzh/EPCOT_gradio/blob/main/gradio.ipynb) to open the example.

To run the demo locally, follow these steps:

####  1. Clone the repository

```
git clone https://github.com/zzh24zzh/EPCOT_gradio.git
```

####  2. Download trained models and reference genome data in hg38 version (download.py)

The pre-trained models will be downloaded to the "models/" directory, and the reference sequence data will be downloaded to the "refSeq/" directory. 
```
python download.py
```


####  3. Prepare input ATAC-seq data (atac_process.py)

 Process an ATAC-seq bam file (hg38 version) to ouput a processed ATAC-seq file in **.pickle** format, which can be uploaded to the demo to excute models.
 It is recommended to generate the .bam file using the [ENCODE ATAC-seq processing pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline).


**Required packages**: 
* deepTools-3.5.1
* samtools-1.16.1
* pyBigWig-0.3.17
```
python atac_process.py -b <ATAC-seq bam file> -p <number of processors>
```



####  4. Launch the Gradio demo

The demo has two interfaces: **(1) Run Model**, and **(2) Visualize Prediction Results**.

In the first interface, users can enter a genomic region and execute models to generate predictions, including

* a file named **"prediction_xxxx.npz"**, which can be uploaded to the second interface for visualization,
* a file named **"formatted_xxxx.zip"**, which contains ChIP-seq and CAGE-seq data in .bigWig format, and contact maps in .bedpe format.

The two files can also be found under the "results/" directory.


```
python gradio_epcot.py
```

After completing these steps, you will be able to use the Gradio interface to run the EPCOT model on your ATAC-seq data and visualize the results.


The predicted modalities include:

| Data       | Standard |
| ----------- | ----------- |
| ChIP-seq (245 epigenomic features)      | arcsinh-transformed signal p-values     |
| CAGE-seq   | log2(x+1)-transformed signal |
| Micro-C   | Observed over expectation count ratio|
| Hi-C| Observed over expectation count ratio|
| ChIP-PET| log2(x+1)-transformed observed over expectation count ratio|
