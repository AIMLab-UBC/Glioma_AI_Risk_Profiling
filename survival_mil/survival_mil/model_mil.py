import torch
import torch.nn as nn
import torch.nn.functional as F


class Attn_Net_Gated(nn.Module):
    def __init__(self, L = 1024, D = 256, dropout = 0.25, n_classes = 1):
        r"""
        Attention Network with Sigmoid Gating (3 fc layers)

        args:
            L (int): input feature dimension
            D (int): hidden layer dimension
            dropout (bool): whether to apply dropout (p = 0.25)
            n_classes (int): number of classes
        """
        super(Attn_Net_Gated, self).__init__()
        self.attention_a = [
            nn.Linear(L, D),
            nn.Tanh()]

        self.attention_b = [nn.Linear(L, D), nn.Sigmoid()]

        self.attention_a.append(nn.Dropout(dropout))
        self.attention_b.append(nn.Dropout(dropout))

        self.attention_a = nn.Sequential(*self.attention_a)
        self.attention_b = nn.Sequential(*self.attention_b)
        self.attention_c = nn.Linear(D, n_classes)

    def forward(self, x):
        a = self.attention_a(x)
        b = self.attention_b(x)
        A = a.mul(b)
        A = self.attention_c(A)  # N x n_classes
        return A, x

class MILSurv(nn.Module):
    def __init__(self, dropout=0.25, n_classes=4,dim=512):
        r"""
        Attention MIL Implementation

        Args:
            size_arg (str): Size of NN architecture (Choices: small or large)
            dropout (float): Dropout rate
            n_classes (int): Output shape of NN
            dim (int): feature embedding dimension
        """
        super(MILSurv, self).__init__()

        fc = [nn.Linear(dim, 512), nn.ReLU(), nn.Dropout(dropout)]
        attention_net = Attn_Net_Gated(L=512, D=256, dropout=dropout, n_classes=1)
        fc.append(attention_net)
        self.attention_net = nn.Sequential(*fc)
        self.rho = nn.Sequential(*[nn.Linear(512, 256), nn.ReLU(), nn.Dropout(dropout)])

        self.classifier = nn.Linear(256, n_classes)



    def forward(self, data):

        A, h_path = self.attention_net(data.squeeze())
        A = torch.transpose(A, 1, 0)
        A = F.softmax(A, dim=1)
        h_path = torch.mm(A, h_path)
        h = self.rho(h_path).squeeze()


        logits  = self.classifier(h).unsqueeze(0)
        Y_hat = torch.topk(logits, 1, dim = 1)[1]
        hazards = torch.sigmoid(logits)
        S = torch.cumprod(1 - hazards, dim=1)

        return hazards, S, Y_hat, A,logits
