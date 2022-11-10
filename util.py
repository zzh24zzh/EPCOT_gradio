import argparse
import numpy as np
import torch
from pretrain.track.model import build_track_model
from cage.model import build_cage_model
from cop.micro_model import build_microc_model
from scipy.sparse import load_npz
# from COP.microc.model import build_pretrain_model_microc
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def parser_args():
    """
    Hyperparameters for the pre-training model
    """
    # add_help = False
    parser = argparse.ArgumentParser(add_help = False)
    parser.add_argument('--num_class', default=245, type=int,help='the number of epigenomic features to be predicted')
    parser.add_argument('--seq_length', default=1600, type=int,help='the length of input sequences')
    parser.add_argument('--nheads', default=4, type=int)
    parser.add_argument('--hidden_dim', default=512, type=int)
    parser.add_argument('--dim_feedforward', default=1024, type=int)
    parser.add_argument('--enc_layers', default=1, type=int)
    parser.add_argument('--dec_layers', default=2, type=int)
    parser.add_argument('--dropout', default=0.2, type=float)
    args, unknown = parser.parse_known_args()
    return args,parser
def get_args():
    args,_ = parser_args()
    return args,_

def parser_args_epi(parent_parser):
    """
    Hyperparameters for the downstream model to predict 1kb-resolution CAGE-seq
    """
    parser=argparse.ArgumentParser(parents=[parent_parser],add_help = False)
    parser.add_argument('--bins', type=int, default=500)
    parser.add_argument('--crop', type=int, default=50)
    parser.add_argument('--embed_dim', default=768, type=int)
    parser.add_argument('--return_embed', default=False, action='store_true')
    args, unknown = parser.parse_known_args()
    return args

def parser_args_cage(parent_parser):
    """
    Hyperparameters for the downstream model to predict 1kb-resolution CAGE-seq
    """
    parser=argparse.ArgumentParser(parents=[parent_parser],add_help = False)
    parser.add_argument('--bins', type=int, default=250)
    parser.add_argument('--crop', type=int, default=25)
    parser.add_argument('--embed_dim', default=360, type=int)
    args, unknown = parser.parse_known_args()
    return args

def parser_args_hic(parent_parser):
    """
    Hyperparameters for the downstream model to predict 5kb-resolution Hi-C and ChIA-PET
    """
    parser=argparse.ArgumentParser(parents=[parent_parser],add_help = False)
    parser.add_argument('--bins', type=int, default=200)
    parser.add_argument('--crop', type=int, default=4)
    parser.add_argument('--pretrain_path', type=str, default='none')
    parser.add_argument('--embed_dim', default=256, type=int)
    parser.add_argument('--trunk',  type=str, default='transformer')
    parser.add_argument('--fine_tune', default=False, action='store_true')
    args, unknown = parser.parse_known_args()
    return args

def parser_args_microc(parent_parser):
    """
    Hyperparameters for the downstream model to predict 1kb-resolution Micro-C
    """
    parser=argparse.ArgumentParser(parents=[parent_parser],add_help = False)
    parser.add_argument('--bins', type=int, default=500)
    parser.add_argument('--crop', type=int, default=10)
    parser.add_argument('--embed_dim', default=768, type=int)
    parser.add_argument('--return_embed', default=True, action='store_false')
    args, unknown = parser.parse_known_args()
    return args


# cage_args=parser_args_cage(parser)
# hic_args=parser_args_hic(parser)
# microc_args=parser_args_microc(parser)



def check_region(chrom,region,ref_genome):
    region=region.replace(',','')
    start,end=list(map(int, region.split('-')))
    if end-start != 1000000:
        raise ValueError("Please enter a 1Mb region")
    if start<300 or end > ref_genome.shape[1]-300:
        raise ValueError("The start of input region should be greater than 300 and "
                         "the end of the region should be less than %s"%(ref_genome.shape[1]-300))
    return int(chrom),start,end

def generate_input(start,end,ref_genome,atac_seq):
    inputs=[]
    for loc in range(start,end,1000):
        tmp_seq = ref_genome[:, loc - 300:loc + 1300]
        tmp_atac = atac_seq[:, loc - 300:loc + 1300]
        inputs.append(np.vstack([tmp_seq,tmp_atac]))
    return np.stack(inputs)

def search_tf(tf):
    with open('epigenomes.txt', 'r') as f:
        epigenomes = f.read().splitlines()
    tf_idx= epigenomes.index(tf)
    return tf_idx


def predict(input_chrom,input_region,input_tf,input_hm,input_file):
    if input_chrom=='' or input_region=='':
        return ValueError("Please specify a genomic region")
    if len(input_file) !=2:
        return ValueError("Please upload two files only, one for the reference genome and one for ATAC-seq")
    with open('epigenomes.txt', 'r') as f:
        epigenomes = f.read().splitlines()
    filea=load_npz(input_file[0]).toarray()
    fileb=load_npz(input_file[1]).toarray()
    ref_genome=filea if filea.shape[0]==4 else fileb
    atac_seq=filea if filea.shape[0]==1 else fileb
    chrom,start,end=check_region(input_chrom,input_region,ref_genome)
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

    epis=[input_tf]+input_hm
    epi_idx=np.array([epigenomes.index(epi) for epi in epis])

    out_epi = predict_epis('models/epi_track.pt', [start, end], ref_genome, atac_seq, device)[:,:,epi_idx]
    print(out_epi.shape)
    # out_cage = predict_cage(model_path, [start, end], ref_genome, atac_seq, device)
    # out_microc = predict_microc(model_path, [start, end], ref_genome, atac_seq, device)

    return 1






def predict_epis(
        model_path,
        region, ref_genome,atac_seq,
        device
):
    args, parser = get_args()
    epi_args = parser_args_epi(parser)
    pretrain_model = build_track_model(epi_args)
    pretrain_model.load_state_dict(torch.load(model_path,map_location=torch.device(device)))
    pretrain_model.eval()
    pretrain_model.to(device)
    inputs=[]
    start,end=region
    for loc in range(start+100000,end-100000,400000):
        inputs.append(generate_input(loc-50000,loc+450000,ref_genome,atac_seq))
    inputs=np.stack(inputs)
    inputs=torch.tensor(inputs).float().to(device)
    pred_epi=[]
    with torch.no_grad():
        for i in range(inputs.shape[0]):
            pred_epi.append(pretrain_model(inputs[i:i+1]).detach().cpu().numpy())
    return np.vstack(pred_epi)

def predict_cage(
        model_path,
        region, ref_genome, atac_seq,
        device
):
    args, parser = get_args()
    cage_args = parser_args_cage(parser)
    cage_model=build_cage_model(cage_args)
    cage_model.load_state_dict(torch.load(model_path,map_location=torch.device(device)))
    cage_model.eval()
    cage_model.to(device)
    inputs = []
    start, end = region
    for loc in range(start + 100000, end - 100000, 200000):
        inputs.append(generate_input(loc - 25000, loc + 225000, ref_genome, atac_seq))
    inputs = np.stack(inputs)
    inputs = torch.tensor(inputs).float().to(device)
    pred_cage=[]
    with torch.no_grad():
        for i in range(inputs.shape[0]):
            pred_cage.append(cage_model(inputs[i:i + 1]).detach().cpu().numpy().squeeze())
    return np.concatenate(pred_cage)

def arraytouptri(arrays,args):
    effective_lens=args.bins-2*args.crop
    triu_tup = np.triu_indices(effective_lens)
    temp=np.zeros((effective_lens,effective_lens))
    temp[triu_tup]=arrays
    return temp
def complete_mat(mat):
    temp = mat.copy()
    np.fill_diagonal(temp,0)
    mat= mat+temp.T
    return mat


# def predict_hic(model,chrom,start,end,dnase,fasta_extractor,cross_cell_type=False):
#     if (end-start)!=1000000:
#         raise ValueError('Please input a 1Mb region')
#     input_dnase = read_dnase_pickle(dnase, chrom[3:])
#     inputs=generate_input(fasta_extractor,chrom,start,end, input_dnase)
#     device = next(model.parameters()).device
#     inputs=inputs.unsqueeze(0).float().to(device)
#     if cross_cell_type:
#         for m in model.modules():
#             if m.__class__.__name__.startswith('BatchNorm'):
#                 m.train()
#     with torch.no_grad():
#         pred_hic = model(inputs).detach().cpu().numpy().squeeze()
#     pred_hic=complete_mat(arraytouptri(pred_hic.squeeze(), hic_args))
#     return pred_hic


def predict_microc(
        model_path,
        region, ref_genome,atac_seq,
        device
):
    args, parser = get_args()
    microc_args = parser_args_microc(parser)
    microc_model = build_microc_model(microc_args)
    microc_model.load_state_dict(torch.load(model_path,map_location=torch.device(device)))
    microc_model.eval()
    microc_model.to(device)
    inputs=[]
    start,end=region
    for loc in range(start+20000,end-20000,480000):
        inputs.append(generate_input(loc-10000,loc+490000,ref_genome,atac_seq))
    inputs=np.stack(inputs)
    inputs=torch.tensor(inputs).float().to(device)
    pred_epi=[]
    with torch.no_grad():
        for i in range(inputs.shape[0]):
            temp=microc_model(inputs[i:i+1]).detach().cpu().numpy().squeeze()
            pred_epi.append(complete_mat(arraytouptri(temp, microc_args)))
    return np.vstack(pred_epi)


def make_plots(out_epis,out_cages,out_microc,start,end,epis):
    num_mod=out_epis.shape[1]+1+1

    fig = plt.figure(figsize=(4, 5.5))
    gs = GridSpec(6, 4, wspace=0, hspace=0.2)



def plot_atac(ax,val,color='#17becf'):
    ax.fill_between(np.arange(val.shape[0]), 0, val, color=color)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.margins(x=0)
    ax.set_xticks([])

def plot_bindings(ax, val, chr, start, end, color='#17becf'):
    ax.fill_between(np.arange(val.shape[0]), 0, val, color=color)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks(np.arange(0, val.shape[0], val.shape[0] // 5))
    ax.set_ylim(0, 1)
    ax.set_xticklabels(np.arange(start, end, (end - start) // 5))
    ax.margins(x=0)

def plot_cage(ax,val,chr,start,end,color='#17becf'):
    ax.fill_between(np.arange(val.shape[0]), 0, val, color=color)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xlabel('%s:%s-%s'%(chr,start,end))
    ax.set_xticks(np.arange(0,val.shape[0],val.shape[0]//5))
    ax.set_xticklabels(np.arange(start,end,(end-start)//5))
    ax.margins(x=0)
def plot_hic(ax, mat,cmap='RdBu_r', vmin=0, vmax=5):
    ax.imshow(mat,cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks([])
    ax.set_yticks([])