import argparse
import numpy as np
import torch
from pretrain.model import build_epd_model
from pretrain.track.model import build_track_model
from cage.model import build_cage_model
from cop.micro_model import build_microc_model
from cop.hic_model import build_hic_model
from einops import rearrange





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
    parser.add_argument('--crop', type=int, default=10)
    parser.add_argument('--embed_dim', default=768, type=int)
    parser.add_argument('--return_embed', default=False, action='store_true')
    args, unknown = parser.parse_known_args()
    return args

def parser_args_cage(parent_parser):
    """
    Hyperparameters for the downstream model to predict 1kb-resolution CAGE-seq
    """
    parser=argparse.ArgumentParser(parents=[parent_parser],add_help = False)
    parser.add_argument('--bins', type=int, default=500)
    parser.add_argument('--crop', type=int, default=10)
    parser.add_argument('--embed_dim', default=768, type=int)
    parser.add_argument('--return_embed', default=True, action='store_false')
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




def check_region(chrom,start,end,ref_genome,region_len):
    start,end=int(start),int(end)
    if end-start != region_len:
        raise ValueError("Please enter a region with the correct length")
    if start<300 or end > ref_genome.shape[1]-300:
        raise ValueError("The start of input region should be greater than 300 and "
                         "the end of the region should be less than %s"%(ref_genome.shape[1]-300))
    return int(chrom),start,end

def generate_input(start,end,ref_genome,atac_seq):
    # inputs=[]
    pad_left=np.expand_dims(np.vstack((ref_genome[:,start-300:start],atac_seq[:,start-300:start])),0)
    pad_right=np.expand_dims(np.vstack((ref_genome[:,end:end+300],atac_seq[:,end:end+300])),0)
    center=np.vstack((ref_genome[:,start:end],atac_seq[:,start:end]))
    center=rearrange(center,'n (b l)-> b n l',l=1000)
    dmatrix = np.concatenate((pad_left, center[:, :, -300:]), axis=0)[:-1, :, :]
    umatrix = np.concatenate((center[:, :, :300], pad_right), axis=0)[1:, :, :]
    return np.concatenate((dmatrix, center, umatrix), axis=2)

    # for loc in range(start,end,1000):
    #     tmp_seq = ref_genome[:, loc - 300:loc + 1300]
    #     tmp_atac = atac_seq[:, loc - 300:loc + 1300]
    #     inputs.append(np.vstack([tmp_seq,tmp_atac]))
    # return np.stack(inputs)

def search_tf(tf):
    with open('data/epigenomes.txt', 'r') as f:
        epigenomes = f.read().splitlines()
    tf_idx= epigenomes.index(tf)
    return tf_idx



def predict_epb(
        model_path,
        region, ref_genome,atac_seq,
        device,
        cop_type
):
    args, parser = get_args()

    pretrain_model = build_epd_model(args)
    pretrain_model.load_state_dict(torch.load(model_path,map_location=torch.device(device)))
    pretrain_model.eval()
    pretrain_model.to(device)
    start,end=region
    inputs=generate_input(start,end,ref_genome,atac_seq)
    inputs=torch.tensor(inputs).float().to(device)
    with torch.no_grad():
        pred_epi=torch.sigmoid(pretrain_model(inputs)).detach().cpu().numpy()
    if cop_type == 'Micro-C (enter a 500 kb region)':
        return pred_epi[10:-10,:]
    else:
        return pred_epi[20:-20,:]


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
        for loc in range(start+20000,end-20000,480000):
            inputs.append(generate_input(loc-10000,loc+490000,ref_genome,atac_seq))
    inputs=np.stack(inputs)
    inputs=torch.tensor(inputs).float().to(device)
    pred_epi=[]
    with torch.no_grad():
        for i in range(inputs.shape[0]):
            pred_epi.append(pretrain_model(inputs[i:i+1]).detach().cpu().numpy())

    out_epi = rearrange(np.vstack(pred_epi), 'i j k -> (i j) k')
    return out_epi

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
        inputs.append(generate_input(start, end, ref_genome, atac_seq))
    else:
        for loc in range(start + 20000, end - 20000, 480000):
            inputs.append(generate_input(loc - 10000, loc + 490000, ref_genome, atac_seq))
    inputs = np.stack(inputs)
    inputs = torch.tensor(inputs).float().to(device)
    pred_cage = []
    with torch.no_grad():
        for i in range(inputs.shape[0]):
            pred_cage.append(cage_model(inputs[i:i + 1]).detach().cpu().numpy().squeeze())
    return np.concatenate(pred_cage)

    # inputs = []
    # start, end = region
    # if cop_type == 'Micro-C (enter a 500 kb region)':
    #     for loc in range(start + 50000, end - 50000, 200000):
    #         inputs.append(generate_input(loc - 25000, loc + 225000, ref_genome, atac_seq))
    # else:
    #     for loc in range(start + 100000, end - 100000, 200000):
    #         inputs.append(generate_input(loc - 25000, loc + 225000, ref_genome, atac_seq))
    # inputs = np.stack(inputs)
    # inputs = torch.tensor(inputs).float().to(device)
    # pred_cage=[]
    # with torch.no_grad():
    #     for i in range(inputs.shape[0]):
    #         pred_cage.append(cage_model(inputs[i:i + 1]).detach().cpu().numpy().squeeze())
    # return np.concatenate(pred_cage)

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


def filetobrowser(out_epis,out_cages,out_cop,chrom,start,end,file_id):
    import pyBigWig,os
    from zipfile import ZipFile
    import zipfile
    import shutil
    import uuid
    with open('data/epigenomes.txt', 'r') as f:
        epigenomes = f.read().splitlines()

    # filename = str(uuid.uuid4())
    files_to_zip = file_id
    if os.path.exists(files_to_zip):
        shutil.rmtree(files_to_zip)
    os.mkdir(files_to_zip)

    hdr=[]
    with open('data/chrom_size_hg38.txt', 'r') as f:
        for line in f:
            tmp=line.strip().split('\t')
            hdr.append((tmp[0],int(tmp[1])))


    for i in range(out_epis.shape[1]):
        bwfile = pyBigWig.open(os.path.join(files_to_zip,"%s.bigWig"%epigenomes[i]), 'w')
        bwfile.addHeader(hdr)
        bwfile.addEntries(['chr' + str(chrom)]*out_epis.shape[0],[loc for loc in range(start,end,1000)],
                          ends=[loc+1000 for loc in range(start,end,1000)],values=out_epis[:,i].tolist())
        bwfile.close()
    bwfile = pyBigWig.open(os.path.join(files_to_zip,"cage.bigWig"),'w')
    bwfile.addHeader(hdr)

    bwfile.addEntries(['chr' + str(chrom)] * out_cages.shape[0], [loc for loc in range(start, end, 1000)],
                      ends=[loc + 1000 for loc in range(start, end, 1000)], values=out_cages.tolist())
    bwfile.close()
    cop_lines=[]

    interval=1000 if out_cop.shape[-1]==480 else 5000
    if out_cop.shape[-1]==480:
        for bin1 in range(out_cop.shape[-1]):
            for bin2 in range(bin1,out_cop.shape[-1],1):
                # tmp=['chr' + str(chrom),str(start+bin1*interval),str(start+(bin1+1)*interval),'chr' + str(chrom),
                #                   str(start + bin2 * interval), str(start + (bin2 + 1) * interval),'.',str(np.around(out_cop[bin1,bin2],2)),'.','.'
                #      ]
                tmp = ['0', 'chr' + str(chrom), str(start + bin1 * interval), '0', '0', 'chr' + str(chrom),
                       str(start + bin2 * interval), '1', str(np.around(out_cop[bin1, bin2], 2))]
                cop_lines.append('\t'.join(tmp)+'\n')
        with open(os.path.join(files_to_zip,"microc.bedpe"),'w') as f:
            f.writelines(cop_lines)
    else:
        types=['CTCF_ChIA-PET','POLR2_ChIA-PET','Hi-C']
        for i in range(len(types)):
            for bin1 in range(out_cop.shape[-1]):
                for bin2 in range(bin1, out_cop.shape[-1], 1):
                    tmp=['0','chr' + str(chrom), str(start + bin1 * interval),'0','0','chr' +str(chrom),str(start + bin2 * interval),'1',str(np.around(out_cop[i,bin1, bin2], 2))]

                    # tmp = ['chr' + str(chrom), str(start + bin1 * interval), str(start + (bin1 + 1) * interval),
                    #        'chr' + str(chrom),
                    #        str(start + bin2 * interval), str(start + (bin2 + 1) * interval), '.',
                    #        str(np.around(out_cop[i,bin1, bin2], 2)), '.', '.'
                    #        ]
                    cop_lines.append('\t'.join(tmp) + '\n')
            with open(os.path.join(files_to_zip,"%s.bedpe"%types[i]), 'w') as f:
                f.writelines(cop_lines)

    out_zipfile = ZipFile("results/formatted_%s.zip" % file_id, "w", zipfile.ZIP_DEFLATED)
    for file_to_zip in os.listdir(files_to_zip):
        file_to_zip_full_path = os.path.join(files_to_zip, file_to_zip)
        out_zipfile.write(filename=file_to_zip_full_path, arcname=file_to_zip)

    out_zipfile.close()
    shutil.rmtree(files_to_zip)
    return "results/formatted_%s.zip"%file_id









