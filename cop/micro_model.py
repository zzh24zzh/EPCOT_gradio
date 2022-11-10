import os,sys
# import inspect
# currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
# parentdir = os.path.dirname(currentdir)
# sys.path.insert(0, parentdir)
import math
from pretrain.track.model import build_track_model
import torch.nn as nn
import torch
import torch.nn.functional as F
from einops.layers.torch import Rearrange
from einops import rearrange
import numpy as np

class Convblock(nn.Module):
    def __init__(self,in_channel,kernel_size,dilate_size,dropout=0.1):
        super().__init__()
        self.conv=nn.Sequential(
            nn.Conv2d(
                in_channel, in_channel,
                kernel_size, padding=self.pad(kernel_size,1)),
            nn.GroupNorm(16, in_channel),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Conv2d(
                in_channel, in_channel,
                kernel_size, padding=self.pad(kernel_size, dilate_size),
                dilation=dilate_size),
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
    def __init__(self,embed_dim,in_channel=64,kernel_size=7,dilate_rate=5):
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

class Downstream_microc_model(nn.Module):
    def __init__(
        self,
        pretrain_model,
        embed_dim,
        hidden_dim=256,
        in_dim=64,
        crop=10,
    ):
        super().__init__()
        self.project = nn.Sequential(
            nn.Linear(embed_dim, 512),
            nn.ReLU(),
            nn.Linear(512, hidden_dim),
        )

        self.pretrain_model=pretrain_model
        self.dilate_tower = dilated_tower(embed_dim=hidden_dim, in_channel=in_dim,dilate_rate=5)
        self.prediction_head = nn.Linear(in_dim, 1)
        self.crop=crop

    def output_head(self, x):
        bins=x.shape[1]
        x1 = torch.tile(x.unsqueeze(1), (1, bins, 1, 1))
        x2 = x1.permute(0, 2, 1, 3)
        mean_out = (x1 + x2) / 2
        dot_out = (x1 * x2)/math.sqrt(x.shape[-1])
        return mean_out + dot_out
    def upper_tri(self, x,bins):
        triu_tup = np.triu_indices(bins)
        d = np.array(list(triu_tup[1] + bins * triu_tup[0]))
        return x[:, d, :]

    def forward(self,x):
        x=self.pretrain_model(x)
        print(x.shape)
        x=self.project(x)
        x = self.output_head(x)
        x = self.dilate_tower(x, self.crop)
        bins = x.shape[1]
        x = rearrange(x, 'b l n d -> b (l n) d')
        x = self.upper_tri(x,bins)
        x = self.prediction_head(x)
        return x



def build_microc_model(args):
    pretrain_model=build_track_model(args)
    model=Downstream_microc_model(
        pretrain_model=pretrain_model,
        embed_dim=args.embed_dim,
        crop=args.crop
    )
    return model

