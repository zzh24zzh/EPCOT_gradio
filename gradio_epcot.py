import gradio as gr
from util import predict
import os
# from math import sqrt
# import matplotlib

# matplotlib.use("Agg")

# import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import load_npz
# import plotly.express as px
# import pandas as pd



with open('epigenomes.txt','r') as f:
    epis=f.read().splitlines()
inputs = [
    gr.Dropdown([str(i) for i in range(1,23)],label='Choose a chromosome'),
    # gr.Text(label='Region start'),
    gr.Text(label='Genomic region in 1 Mb'),
    gr.Dropdown(epis[:-9],label='Select TFs'),
    gr.CheckboxGroup(epis[-9:],label='Select histone marks'),
    # gr.Dropdown(['None','CTCF','RAD21'],label='Choose the epigenomic feature to be predicted'),
    gr.Files(label="Upload the ATAC-seq and reference genome files of the chosen chromosome"),
]
outputs = gr.Text()

demo = gr.Interface(
    fn=predict,
    inputs=inputs,
    outputs=outputs,
    examples=[
        ["11","7750000-8750000", "CTCF",["H3K4me2","H3K4me3"],["examples/chr11.npz","examples/GM12878_atac_chr11.npz"]],
    ],
    # cache_examples=True,
)


if __name__ == "__main__":
    demo.launch()
