import os,sys
from pretrain.track.layers import AttentionPool,CNN
from pretrain.track.transformers import Transformer
from einops.layers.torch import Rearrange
from einops import rearrange
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np

import torchvision.transforms as T

class Convblock(nn.Module):
    def __init__(self,in_channel,kernel_size,dilate_size,dropout=0.1):
        super().__init__()
        self.conv=nn.Sequential(
            nn.Conv2d(
                in_channel, in_channel,
                kernel_size, padding=self.pad(kernel_size, dilate_size),
                dilation=dilate_size),
            nn.GroupNorm(16, in_channel),
            nn.Dropout(dropout)
        )
    def pad(self,kernelsize, dialte_size):
        return (kernelsize - 1) * dialte_size // 2
    def symmetric(self,x):
        return (x + x.permute(0,1,3,2)) / 2
    def forward(self,x):
        identity=x
        out=self.conv(x)
        x=out+identity
        x=self.symmetric(x)
        return F.relu(x)

class dilated_tower(nn.Module):
    def __init__(self,embed_dim,in_channel=48,kernel_size=9,dilate_rate=4):
        super().__init__()
        dilate_convs=[]
        for i in range(dilate_rate+1):
            dilate_convs.append(
                Convblock(in_channel,kernel_size=kernel_size,dilate_size=2**i))

        self.cnn=nn.Sequential(
            Rearrange('b l n d -> b d l n'),
            nn.Conv2d(embed_dim, in_channel, kernel_size=1),
            *dilate_convs,
            nn.Conv2d(in_channel, in_channel, kernel_size=1),
            Rearrange('b d l n -> b l n d'),
        )
    def forward(self,x,crop):
        x=self.cnn(x)
        x=x[:,crop:-crop,crop:-crop,:]
        return x

class discriminator(nn.Module):
    def __init__(self,):
        super(discriminator, self).__init__()
        self.conv=nn.Sequential(
            nn.Conv2d(1,32,kernel_size=7,padding=3),
            nn.ReLU(),
            nn.MaxPool1d(4),
            nn.Conv2d(32,64,kernel_size=7,padding=3),
            nn.ReLU(),
            nn.MaxPool1d(4),
            nn.Conv2d(64, 64, kernel_size=1),
        )
        self.linear=nn.Linear(64,1)

    def complete_mat(self,x,smooth=True):
        tmp=x.copy()
        tmp.fill_diagonal_(0)
        if smooth:
            t = T.GaussianBlur(kernel_size=5, sigma=0.5)
            tmp = t(tmp.T + x)
        else:
            tmp=tmp.T + x
        return tmp
    def forward(self,x1,x2):
        x1,x2=self.complete_mat(x1),self.complete_mat(x2)
        diff_mat=(x1-x2)**2
        diff_mat=self.conv(diff_mat).mean(dim=(-2, -1))
        return self.linear(diff_mat)



class Tranmodel(nn.Module):
    def __init__(self, backbone, transfomer):
        super().__init__()
        self.backbone = backbone
        self.transformer = transfomer
        hidden_dim = transfomer.d_model
        self.input_proj = nn.Conv1d(backbone.num_channels, hidden_dim, kernel_size=1)
    def forward(self, input):
        input=rearrange(input,'b n c l -> (b n) c l')
        src = self.backbone(input)
        src=self.input_proj(src)
        src = self.transformer(src)
        return src

class finetunemodel(nn.Module):
    def __init__(
        self,
        pretrain_model,
        hidden_dim,
        embed_dim,
        device,
        bins=200,
        in_dim=64,
        max_bin=10,
        crop=4,
        output_dim=1
    ):
        super().__init__()
        self.pretrain_model = pretrain_model
        self.bins = bins
        self.max_bin = max_bin
        self.attention_pool = AttentionPool(hidden_dim)
        self.crop = crop
        self.project = nn.Sequential(
            Rearrange('(b n) c -> b c n', n=bins * 5),
            nn.Conv1d(hidden_dim, hidden_dim, kernel_size=15, padding=7, groups=hidden_dim),
            nn.InstanceNorm1d(hidden_dim, affine=True),
            nn.Conv1d(hidden_dim, embed_dim, kernel_size=1),
            nn.ReLU(inplace=True),
            nn.Dropout(0.2)
        )

        self.cnn = nn.Sequential(
            nn.Conv1d(embed_dim, embed_dim, kernel_size=15, padding=7),
            nn.GroupNorm(32, embed_dim),
            nn.MaxPool1d(kernel_size=5, stride=5),
            nn.ReLU(inplace=True),
            nn.Conv1d(embed_dim, embed_dim, kernel_size=1),
            nn.Dropout(0.2),
            Rearrange('b c n -> b n c')
        )
        encoder_layer = nn.TransformerEncoderLayer(d_model=embed_dim, nhead=4, dim_feedforward=2 * embed_dim,
                                                   batch_first=True, norm_first=True)
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=4)

        self.distance_embed = nn.Embedding(max_bin + 1, embed_dim)

        self.dilate_tower = dilated_tower(embed_dim=embed_dim, in_channel=in_dim)
        self.prediction_head1 = nn.Linear(in_dim, output_dim)
        self.dist_dropout = nn.Dropout(0.1)
        self.device = device

    def position_matrix(self, bins, b, maxbin):
        pos1 = np.tile(np.arange(bins), (bins, 1))
        pos2 = pos1.T
        pos = np.abs(pos1 - pos2)
        pos = np.where(pos > maxbin, maxbin, pos)
        pos = np.tile(pos, (b, 1, 1))
        return torch.tensor(pos).long().to(self.device)

    def upper_tri(self, x, bins):
        triu_tup = np.triu_indices(bins)
        d = np.array(list(triu_tup[1] + bins * triu_tup[0]))
        return x[:, d, :]

    def output_head(self, x, dist_embed, bins):
        x1 = torch.tile(x.unsqueeze(1), (1, bins, 1, 1))
        x2 = x1.permute(0, 2, 1, 3)
        mean_out = (x1 + x2) / 2
        dot_out = x1 * x2
        return mean_out + dot_out + dist_embed

    def forward(self, x):
        b = x.shape[0]
        x = self.pretrain_model(x)
        x = self.attention_pool(x)
        x = self.project(x)
        x = self.cnn(x)
        x = self.transformer(x)
        dist_embed = self.dist_dropout(self.distance_embed(self.position_matrix(self.bins, b=b, maxbin=self.max_bin)))
        x = self.output_head(x, dist_embed, self.bins)
        x = self.dilate_tower(x, self.crop)
        x = rearrange(x, 'b l n d -> b (l n) d')
        x = self.upper_tri(x, self.bins - 2 * self.crop)
        x = self.prediction_head1(x)
        return x

def build_backbone():
    model = CNN()
    return model
def build_transformer(args):
    return Transformer(
        d_model=args.hidden_dim,
        dropout=args.dropout,
        nhead=args.nheads,
        dim_feedforward=args.dim_feedforward,
        num_encoder_layers=args.enc_layers,
        num_decoder_layers=args.dec_layers
    )


def build_hic_model(args):
    backbone = build_backbone()
    transformer = build_transformer(args)
    pretrain_model = Tranmodel(
            backbone=backbone,
            transfomer=transformer,
        )

    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    model=finetunemodel(pretrain_model,hidden_dim=args.hidden_dim,embed_dim=args.embed_dim,device=device,bins=args.bins,crop=args.crop,output_dim=3)
    return model