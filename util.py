import argparse
import numpy as np
import torch
from pretrain.track.model import build_track_model
from cage.model import build_cage_model
from cop.micro_model import build_microc_model
from cop.hic_model import build_hic_model






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
    parser.add_argument('--embed_dim', default=256, type=int)
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




def check_region(chrom,region,ref_genome,region_len):
    region=region.replace(',','')
    start,end=list(map(int, region.split('-')))
    if end-start != region_len:
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






def predict_epis(
        model_path,
        region, ref_genome,atac_seq,
        device,
        cop_type
):
    args, parser = get_args()
    epi_args = parser_args_epi(parser)
    pretrain_model = build_track_model(epi_args)
    pretrain_model.load_state_dict(torch.load(model_path,map_location=torch.device(device)))
    pretrain_model.eval()
    pretrain_model.to(device)
    inputs=[]
    start,end=region
    if cop_type == 'Micro-C (enter a 500 kb region)':
        inputs.append(generate_input(start,end,ref_genome,atac_seq))
    else:
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
        device,
        cop_type
):
    args, parser = get_args()
    cage_args = parser_args_cage(parser)
    cage_model=build_cage_model(cage_args)
    cage_model.load_state_dict(torch.load(model_path,map_location=torch.device(device)))
    cage_model.eval()
    cage_model.to(device)
    inputs = []
    start, end = region
    if cop_type == 'Micro-C (enter a 500 kb region)':
        for loc in range(start + 50000, end - 50000, 200000):
            inputs.append(generate_input(loc - 25000, loc + 225000, ref_genome, atac_seq))
    else:
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


def predict_hic(
        model_path,
        region, ref_genome,atac_seq,
        device
):
    args, parser = get_args()
    hic_args = parser_args_hic(parser)
    hic_model = build_hic_model(hic_args)
    hic_model.load_state_dict(torch.load(model_path,map_location=torch.device(device)))
    hic_model.eval()
    hic_model.to(device)
    start,end=region
    inputs=np.stack([generate_input(start,end,ref_genome,atac_seq)])
    inputs=torch.tensor(inputs).float().to(device)
    with torch.no_grad():
        temp=hic_model(inputs).detach().cpu().numpy().squeeze()
    return np.stack([complete_mat(arraytouptri(temp[:,i], hic_args)) for i in range(temp.shape[-1])])


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
    start,end=region
    inputs=np.stack([generate_input(start,end,ref_genome,atac_seq)])
    inputs=torch.tensor(inputs).float().to(device)
    with torch.no_grad():
        temp=microc_model(inputs).detach().cpu().numpy().squeeze()
    return complete_mat(arraytouptri(temp, microc_args))


def filetobrowser(out_epis,out_cages,out_cop,chrom,start,end):
    import pyBigWig,os
    from zipfile import ZipFile
    with open('epigenomes.txt', 'r') as f:
        epigenomes = f.read().splitlines()

    out_zipfile=ZipFile("tmps/browser_data/tmp_%s-%s-%s.zip"%(chrom,start,end), "w")
    hdr=[]
    with open('examples/chrom_size_hg38.txt','r') as f:
        for line in f:
            tmp=line.strip().split('\t')
            hdr.append((tmp[0],int(tmp[1])))
        # bwOutput.addHeader([('chr18', 80373285)])

    for i in range(out_epis.shape[1]):
        bwfile = pyBigWig.open("tmps/browser_data/%s.bigWig"%epigenomes[i], 'w')
        bwfile.addHeader(hdr)
        # for idx in range((end-start)//1000):
        bwfile.addEntries(['chr' + str(chrom)]*out_epis.shape[0],[loc for loc in range(start,end,1000)],
                          ends=[loc+1000 for loc in range(start,end,1000)],values=out_epis[:,i].tolist())
        bwfile.close()
        out_zipfile.write("tmps/browser_data/%s.bigWig"%epigenomes[i])
        os.remove("tmps/browser_data/%s.bigWig"%epigenomes[i])
    bwfile = pyBigWig.open("tmps/browser_data/cage.bigWig",'w')
    bwfile.addHeader(hdr)

    bwfile.addEntries(['chr' + str(chrom)] * out_cages.shape[0], [loc for loc in range(start, end, 1000)],
                      ends=[loc + 1000 for loc in range(start, end, 1000)], values=out_cages.tolist())
    bwfile.close()
    out_zipfile.write("tmps/browser_data/cage.bigWig")
    os.remove("tmps/browser_data/cage.bigWig")
    cop_lines=[]

    interval=1000 if out_cop.shape[-1]==400 else 5000
    if out_cop.shape[-1]==400:
        for bin1 in range(out_cop.shape[-1]):
            for bin2 in range(bin1,out_cop.shape[-1],1):
                tmp=['chr' + str(chrom),str(start+bin1*interval),str(start+(bin1+1)*interval),'chr' + str(chrom),
                                  str(start + bin2 * interval), str(start + (bin2 + 1) * interval),'.',str(np.around(out_cop[bin1,bin2],2)),'.','.'
                     ]
                cop_lines.append('\t'.join(tmp)+'\n')
        with open("tmps/browser_data/microc.bedpe",'w') as f:
            f.writelines(cop_lines)
        out_zipfile.write("tmps/browser_data/microc.bedpe")
        os.remove("tmps/browser_data/microc.bedpe")
    else:
        types=['CTCF ChIA-PET','POLR2 ChIA-PET','Hi-C']
        for i in range(3):
            for bin1 in range(out_cop.shape[-1]):
                for bin2 in range(bin1, out_cop.shape[-1], 1):
                    tmp = ['chr' + str(chrom), str(start + bin1 * interval), str(start + (bin1 + 1) * interval),
                           'chr' + str(chrom),
                           str(start + bin2 * interval), str(start + (bin2 + 1) * interval), '.',
                           str(np.around(out_cop[i,bin1, bin2], 2)), '.', '.'
                           ]
                    cop_lines.append('\t'.join(tmp) + '\n')
            with open("tmps/browser_data/%s.bedpe"%types[i], 'w') as f:
                f.writelines(cop_lines)
            out_zipfile.write("tmps/browser_data/%s.bedpe"%types[i])
            os.remove("tmps/browser_data/%s.bedpe"%types[i])

    out_zipfile.close()
    return "tmps/browser_data/tmp_%s-%s-%s.zip"%(chrom,start,end)











