import os,sys
import math
from pretrain.track.model import build_track_model
import torch.nn as nn


class Downstream_cage_model(nn.Module):
    def __init__(self,pretrain_model,embed_dim,crop):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(embed_dim, 128),
            nn.ReLU(),
            nn.Linear(128,1)
        )
        self.pretrain_model=pretrain_model
        self.crop=crop
    def forward(self,x):
        x=self.pretrain_model(x)
        out=self.mlp(x[:,self.crop:-self.crop,:])
        return out

def build_cage_model(args):
    pretrain_model=build_track_model(args)
    model=Downstream_cage_model(
        pretrain_model=pretrain_model,
        embed_dim=args.embed_dim,
        crop=args.crop
    )
    return model












# import os,sys
# # import inspect
# # currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
# # parentdir = os.path.dirname(currentdir)
# # sys.path.insert(0, parentdir)
# from pretrain.track.layers import AttentionPool,Enformer,CNN
# from pretrain.track.transformers import Transformer
# from einops.layers.torch import Rearrange
# from einops import rearrange
# import torch
# import torch.nn as nn
# import torch.nn.functional as F
#
#
# class Convblock(nn.Module):
#     def __init__(self,in_channel,kernel_size,dilate_size,dropout=0.1):
#         super().__init__()
#         self.conv=nn.Sequential(
#             nn.Conv2d(
#                 in_channel, in_channel,
#                 kernel_size, padding=self.pad(kernel_size, dilate_size),
#                 dilation=dilate_size),
#             nn.GroupNorm(16, in_channel),
#             nn.Dropout(dropout)
#         )
#     def pad(self,kernelsize, dialte_size):
#         return (kernelsize - 1) * dialte_size // 2
#     def symmetric(self,x):
#         return (x + x.permute(0,1,3,2)) / 2
#     def forward(self,x):
#         identity=x
#         out=self.conv(x)
#         x=out+identity
#         x=self.symmetric(x)
#         return F.relu(x)
#
# class dilated_tower(nn.Module):
#     def __init__(self,embed_dim,in_channel=48,kernel_size=9,dilate_rate=4):
#         super().__init__()
#         dilate_convs=[]
#         for i in range(dilate_rate+1):
#             dilate_convs.append(
#                 Convblock(in_channel,kernel_size=kernel_size,dilate_size=2**i))
#
#         self.cnn=nn.Sequential(
#             Rearrange('b l n d -> b d l n'),
#             nn.Conv2d(embed_dim, in_channel, kernel_size=1),
#             *dilate_convs,
#             nn.Conv2d(in_channel, in_channel, kernel_size=1),
#             Rearrange('b d l n -> b l n d'),
#         )
#     def forward(self,x,crop):
#         x=self.cnn(x)
#         x=x[:,crop:-crop,crop:-crop,:]
#         return x
#
#
# class Tranmodel(nn.Module):
#     def __init__(self, backbone, transfomer):
#         super().__init__()
#         self.backbone = backbone
#         self.transformer = transfomer
#         hidden_dim = transfomer.d_model
#         self.input_proj = nn.Conv1d(backbone.num_channels, hidden_dim, kernel_size=1)
#     def forward(self, input):
#         input=rearrange(input,'b n c l -> (b n) c l')
#         src = self.backbone(input)
#         src=self.input_proj(src)
#         src = self.transformer(src)
#         return src
#
# class finetunemodel(nn.Module):
#     def __init__(self, pretrain_model, hidden_dim, embed_dim, bins, crop=25):
#         super().__init__()
#         self.pretrain_model = pretrain_model
#         self.bins = bins
#         self.crop = crop
#         self.attention_pool = AttentionPool(hidden_dim)
#         self.project = nn.Sequential(
#             Rearrange('(b n) c -> b c n', n=bins),
#             nn.Conv1d(hidden_dim, hidden_dim, kernel_size=9, padding=4, groups=hidden_dim),
#             nn.InstanceNorm1d(hidden_dim, affine=True),
#             nn.Conv1d(hidden_dim, embed_dim, kernel_size=1),
#             nn.ReLU(inplace=True),
#             nn.Dropout(0.2)
#         )
#         self.transformer = Enformer(dim=embed_dim, depth=4, heads=6)
#         self.prediction_head = nn.Sequential(
#             nn.Linear(embed_dim, 1)
#         )
#
#
#     def forward(self, x):
#         # x = rearrange(x, 'b n c l -> (b n) c l')
#         x = self.pretrain_model(x)
#         x = self.attention_pool(x)
#         x = self.project(x)
#         x = rearrange(x, 'b c n -> b n c')
#         x = self.transformer(x)
#         x = self.prediction_head(x[:, self.crop:-self.crop, :])
#         return x
#
# def build_backbone():
#     model = CNN()
#     return model
# def build_transformer(args):
#     return Transformer(
#         d_model=args.hidden_dim,
#         dropout=args.dropout,
#         nhead=args.nheads,
#         dim_feedforward=args.dim_feedforward,
#         num_encoder_layers=args.enc_layers,
#         num_decoder_layers=args.dec_layers
#     )
# def build_cage_model(args):
#     backbone = build_backbone()
#     transformer = build_transformer(args)
#     pretrain_model = Tranmodel(
#             backbone=backbone,
#             transfomer=transformer,
#         )
#
#     model=finetunemodel(pretrain_model,hidden_dim=args.hidden_dim,embed_dim=args.embed_dim,
#                         bins=args.bins,crop=args.crop)
#     return model