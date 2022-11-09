import torch
import torch.nn as nn
import numpy as np
import torch.nn.functional as F


class CNN(nn.Module):
    def __init__(self):
        super(CNN, self).__init__()
        conv_kernel_size1 = 10
        conv_kernel_size2 = 8
        pool_kernel_size1 = 5
        pool_kernel_size2 = 4
        self.conv_net = nn.Sequential(
            nn.Conv1d(5, 256, kernel_size=conv_kernel_size1),
            nn.ReLU(inplace=True),
            nn.Dropout(p=0.1),
            nn.Conv1d(256, 256, kernel_size=conv_kernel_size1),
            # nn.GroupNorm(16, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(kernel_size=pool_kernel_size1, stride=pool_kernel_size1),
            nn.Dropout(p=0.1),
            nn.Conv1d(256, 360, kernel_size=conv_kernel_size2),
            nn.ReLU(inplace=True),
            nn.Dropout(p=0.1),
            nn.Conv1d(360, 360, kernel_size=conv_kernel_size2),
            nn.BatchNorm1d(360),
            # nn.GroupNorm(36, 360),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(kernel_size=pool_kernel_size2, stride=pool_kernel_size2),
            nn.Dropout(p=0.1),
            nn.Conv1d(360, 512, kernel_size=conv_kernel_size2),
            nn.ReLU(inplace=True),
            nn.Dropout(p=0.2),
            nn.Conv1d(512, 512, kernel_size=conv_kernel_size2),
            nn.BatchNorm1d(512),
            # nn.GroupNorm(32, 512),
            nn.ReLU(inplace=True),
            nn.Dropout(p=0.2))
        self.num_channels = 512
    def forward(self, x):
        out = self.conv_net(x)
        return out


# borrow from https://github.com/Alibaba-MIIL/ASL/blob/main/src/loss_functions/losses.py
class Balanced_AsymmetricLoss(nn.Module):
    def __init__(self, gamma_neg=4, gamma_pos=1, clip=0.05, alpha=None, eps=1e-8, disable_torch_grad_focal_loss=True):
        super(Balanced_AsymmetricLoss, self).__init__()

        self.gamma_neg = gamma_neg
        self.gamma_pos = gamma_pos
        self.clip = clip
        self.disable_torch_grad_focal_loss = disable_torch_grad_focal_loss
        self.eps = eps
        self.alpha = alpha

    def forward(self, x, y, mask):
        # Calculating Probabilities
        assert y.shape == mask.shape
        x_sigmoid = torch.sigmoid(x)
        xs_pos = x_sigmoid
        xs_neg = 1 - x_sigmoid
        # Asymmetric Clipping
        if self.clip is not None and self.clip > 0:
            xs_neg = (xs_neg + self.clip).clamp(max=1)

        # Basic CE calculation
        los_pos = y * torch.log(xs_pos.clamp(min=self.eps))
        los_neg = (1 - y) * torch.log(xs_neg.clamp(min=self.eps))
        if self.alpha is not None:
            los_pos = self.alpha * los_pos
        loss = los_pos + los_neg
        # Asymmetric Focusing
        if self.gamma_neg > 0 or self.gamma_pos > 0:
            if self.disable_torch_grad_focal_loss:
                torch.set_grad_enabled(False)
            pt0 = xs_pos * y
            pt1 = xs_neg * (1 - y)  # pt = p if t > 0 else 1-p
            pt = pt0 + pt1
            one_sided_gamma = self.gamma_pos * y + self.gamma_neg * (1 - y)
            one_sided_w = torch.pow(1 - pt, one_sided_gamma)
            if self.disable_torch_grad_focal_loss:
                torch.set_grad_enabled(True)
            loss *= one_sided_w
        loss *= mask
        return -loss.sum() / (torch.sum(mask) + self.eps)


