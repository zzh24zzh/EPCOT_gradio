import gradio as gr
from func_gradio import predict_func,make_plots
import os


with open(os.path.abspath('epigenomes.txt'),'r') as f:
    epis=f.read().splitlines()
inputs = [
    gr.Dropdown([str(i) for i in range(1,23)],label='Choose a chromosome'),
    gr.Dropdown(['Micro-C (enter a 500 kb region)', 'Hi-C (enter a 1 Mb region)']
                ,label='Choose the type of chromatin contact maps to be predicted'),
    gr.Text(label='Enter a region of interest'),

    gr.Files(label="Upload the ATAC-seq and reference genome files of the chosen chromosome"),
]
outputs = [
    gr.File(label='Download the prediction files for visualization in the app'),
    gr.File(label='Download the prediction files for browser')
]

app1 = gr.Interface(
    fn=predict_func,
    inputs=inputs,
    outputs=outputs,
    examples=[
        ["11","Micro-C (enter a 500 kb region)","10500000-11000000",
         [os.path.join(os.path.dirname(__file__),"examples/chr11.npz"),
          os.path.join(os.path.dirname(__file__),"examples/GM12878_atac_chr11.npz")]
         ],
        ["11","Hi-C (enter a 1 Mb region)","7750000-8750000",
         [os.path.join(os.path.dirname(__file__),"examples/chr11.npz"),
          os.path.join(os.path.dirname(__file__),"examples/GM12878_atac_chr11.npz")]
         ]
    ],
)

inputs = [
    gr.File(label="Upload the prediction files under 'tmps' folder"),
    gr.Slider(maximum=16,label="Range of contact map values",value=6),
    gr.Slider(maximum=10,label="Range of epigenomic feature track values",value=6),
    gr.CheckboxGroup(epis,label='Choose epigenomes to be visualized'),

]
outputs = gr.Plot(label='Prediction visualization')

app2 = gr.Interface(
    fn=make_plots,
    inputs=inputs,
    outputs=outputs,
    examples=[
        [os.path.join(os.path.dirname(__file__),"examples/GM12878_prediction_11-10550000-10950000.npz"),"6","6",["CTCF","H3K4me3"]],
         [os.path.join(os.path.dirname(__file__),"examples/prediction_11-7850000-8650000.npz"),"3","6",["CTCF","POLR2A","H3K4me3"]],
    ],
)
demo = gr.TabbedInterface([app1, app2], ["Modality prediction", "Prediction visualization"])

if __name__ == "__main__":
    demo.launch()
