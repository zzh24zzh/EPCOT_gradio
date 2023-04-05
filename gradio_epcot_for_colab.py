import gradio as gr
import os
from func_gradio import predict_func_for_colab,make_plots

inputs = [
    gr.Dropdown([str(i) for i in range(1,23)],label='Chromosome',default='1'),
    gr.Dropdown(['Micro-C', 'Hi-C (ChIA-PET)']
                ,label='Chromatin contact map', info='One type of contact map is predicted for each time'),
    gr.Number(label='Region of interest (500kb for Micro-C and 1Mb for Hi-C)',info='From'),
    gr.Number(info='To',show_label=False),
    gr.Textbox(
            label="ATAC-seq file",
            info="Path to the processed ATAC-seq file",
            lines=1,
            # value="examples/atac_GM12878.pickle",
        ),
]
outputs = [
    gr.Files(label='Download the results'),
]
app1 = gr.Interface(
    fn=predict_func_for_colab,
    inputs=inputs,
    outputs=outputs,
    title='A computational tool to use ATAC-seq to impute epigenome, transcriptome, and high-resolution chromatin contact maps',
    allow_flagging='manual',
    description='<a href="https://github.com/zzh24zzh/EPCOT_gradio" class="built-with svelte-1lyswbr" target="_blank" '
                'style="font-size: 15px; font-color: black; font-weight:bold" rel="noreferrer">'
                'View Documentation </a>',
    # examples=[["11","Micro-C","10500000","11000000","examples/atac_GM12878.pickle"]]
)


with open(os.path.abspath('data/epigenomes.txt'), 'r') as f:
    epis=f.read().splitlines()
inputs1 = [
    gr.File(label="Prediction file (in .npz format))"),
    gr.Markdown(value='### Visualization options'),
    gr.Dropdown(epis,label='Epigenome features',multiselect=True,max_choices=8,value=['CTCF','H3K4me3']),
    gr.Radio(choices=['Signal p-values (archsinh)','Binding probability'], label='Choose the type of epigenomic feature data'
             , value='Signal p-values (archsinh)'),
    gr.Slider(maximum=16,label="Range of values displayed on the contact map plot",value=6),
    gr.Slider(maximum=12,label="Range of values displayed on the epigenomic feature plot",value=4),
    gr.Slider(maximum=12,label="Range of values displayed on the CAGE-seq plot",value=8),
]
outputs1 = gr.Plot(label='Plots')
app2 = gr.Interface(
    fn=make_plots,
    inputs=inputs1,
    outputs=outputs1,
    live=True,
    allow_flagging='manual'
)

demo = gr.TabbedInterface([app1, app2], ["Run model", "Visualize prediction results"],
                          theme=gr.themes.Soft())
if __name__=='__main__':
    demo.launch(debug=True)