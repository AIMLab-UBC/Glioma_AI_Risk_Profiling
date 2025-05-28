import numpy as np

import torch
import torch.optim as optim
from torch.utils.data import DataLoader, Sampler,RandomSampler, SequentialSampler

device=torch.device("cuda" if torch.cuda.is_available() else "cpu")

def get_optim(model, args):
    if args.opt == "adam":
        optimizer = optim.Adam(filter(lambda p: p.requires_grad, model.parameters()), lr=args.lr, betas=(0.9, 0.999), weight_decay=args.reg)
    elif args.opt == 'sgd':
        optimizer = optim.SGD(filter(lambda p: p.requires_grad, model.parameters()), lr=args.lr, momentum=0.9, weight_decay=args.reg,nesterov=True)
    elif args.opt == 'nadam':
        optimizer = optim.NAdam(filter(lambda p: p.requires_grad, model.parameters()), lr=args.lr, weight_decay=args.reg)
    elif args.opt == 'adamw':
        optimizer = optim.AdamW(filter(lambda p: p.requires_grad, model.parameters()), lr=args.lr,betas=(0.9, 0.999), weight_decay=args.reg)
    else:
        raise NotImplementedError
    return optimizer

def print_network(net):
    params = 0
    params_train = 0

    for param in net.parameters():
        n = param.numel()
        if param.requires_grad:
            params_train += n
        params += n
    print(net)
    print('Total number of parameters: %d' % params)
    print('Total number of trainable parameters: %d' % params_train)

def collate_survival(batch):
    img = torch.cat([item[0] for item in batch], dim = 0)
    label = torch.LongTensor([item[1] for item in batch])
    event_time = np.array([item[2] for item in batch])
    c = torch.FloatTensor([item[3] for item in batch])

    return [img, label, event_time, c]


def get_split_loader(split_dataset,n_workers,training = False, batch_size=1):
    collate = collate_survival
    kwargs = {'num_workers': n_workers} if device.type == "cuda" else {}
    if training:
        loader = DataLoader(split_dataset, batch_size=batch_size, sampler = RandomSampler(split_dataset), collate_fn = collate, **kwargs)
    else:
        loader = DataLoader(split_dataset, batch_size=batch_size, sampler = SequentialSampler(split_dataset), collate_fn = collate, **kwargs)
    return loader
