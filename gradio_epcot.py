import gradio as gr
import os
from func_gradio import predict_func,make_plots

inputs = [
    gr.Dropdown([str(i) for i in range(1,23)],label='Choose a chromosome'),
    gr.Dropdown(['Micro-C', 'Hi-C (ChIA-PET)']
                ,label='Choose a type of chromatin contact maps to be predicted', info='One contact map type is predicted for each time'),
    gr.Number(label='Enter a 500kb genomic region if Hi-C is selected, or a 1Mb genomic region if Micro-C is selected.',info='From'),
    gr.Number(info='To',show_label=False),
    gr.File(label='Upload or Drop the processed ATAC-seq file (in **.pickle** format)'),
    # gr.Examples(["11","Micro-C","10500000","11000000",
    #       os.path.join(os.path.dirname(__file__),"atac_36edf5f4-8824-4e1d-b031-451c5df505f0.pickle")])

]

outputs = [
    gr.Files(label='Download the prediction results once model running is completed.'),
]



app1 = gr.Interface(
    fn=predict_func,
    inputs=inputs,
    outputs=outputs,
    title='A computational tool to use ATAC-seq to impute epigenome, transcriptome, and high-resolution chromatin contact maps',
    allow_flagging='manual',
    description='<a href="https://github.com/zzh24zzh/EPCOT_gradio" class="built-with svelte-1lyswbr" target="_blank" '
                'style="font-size: 15px; font-color: black; font-weight:bold" rel="noreferrer">'
                'View Documentation </a>',
    interpretation='explain',
)


with open(os.path.abspath('data/epigenomes.txt'), 'r') as f:
    epis=f.read().splitlines()
inputs1 = [
    gr.File(label="Upload the prediction file (in .npz format))"),
    gr.Markdown(value='### Visualization options'),
    gr.Dropdown(epis,label='Select epigenome features to be plotted',multiselect=True,max_choices=8,value=['CTCF','H3K4me3']),
    gr.Radio(choices=['Signal (arcsinh signal p-values)','Binding probability'], label='Choose the type of epigenomic feature data'
             , value='Signal (arcsinh signal p-values)'),
    gr.Slider(maximum=16,label="Set the range of values displayed on the contact map plot",value=6),
    gr.Slider(maximum=12,label="Set the range of values displayed on the epigenomic feature plot",value=4),
    gr.Slider(maximum=12,label="Set the range of values displayed on the CAGE-seq plot",value=8),
]
outputs1 = gr.Plot(label='Plots of predicted modalities')

app2 = gr.Interface(
    fn=make_plots,
    inputs=inputs1,
    outputs=outputs1,
    live=True,
    allow_flagging='manual'
)

demo = gr.TabbedInterface([app1, app2], ["Run model", "Visualize prediction results"],
                          css=".gradio-container-3-24-1 {font-weight:bold}"
                              "h1 {font-family: verdana; color:grey} div {border-color: black}"
                                ".selected.svelte-1g805jl {border-color: white;border-width:2.5px; font-weight:bold;border-bottom-width:1px} "
                              ".tab-nav.svelte-1g805jl{border-bottom-color: white;border-bottom-width:2.5px}"
                              "div.tabitem {border-color: white;border-width:2.5px}" 
                              "button.svelte-1g805jl{font-weight:bold}"
                              "body{background-color: white;}"
                ".gradio-container{background-color: #e5e5f7;"
                          "background-image:  repeating-radial-gradient( circle at 0 0, transparent 0, #e5e5f7 10px ),"
                          " repeating-linear-gradient(#fffcdc,#d9a7c7) ",
                          theme=gr.themes.Soft())


if __name__ == "__main__":
    demo.launch(share=True,debug=True)
