from argparse import Namespace
from collections import OrderedDict
import os
import pickle

import numpy as np
from sksurv.metrics import concordance_index_censored

import torch
from torch.utils.data import DataLoader

from . import model_mil
from . import utils


def nll_loss(hazards, S, Y, c, alpha=0.4, eps=1e-7):
    batch_size = len(Y)
    Y = Y.view(batch_size, 1)
    c = c.view(batch_size, 1).float()
    if S is None:
        S = torch.cumprod(1 - hazards, dim=1)
    S_padded = torch.cat([torch.ones_like(c), S], 1)
    Y=Y.type(torch.int64)
    uncensored_loss = -(1 - c) * (torch.log(torch.gather(S_padded, 1, Y).clamp(min=eps)) + torch.log(torch.gather(hazards, 1, Y).clamp(min=eps)))
    censored_loss = - c * torch.log(torch.gather(S_padded, 1, Y+1).clamp(min=eps))
    neg_l = censored_loss + uncensored_loss
    loss = (1-alpha) * neg_l + alpha * uncensored_loss
    return loss

class NLLSurvLoss(object):
    def __init__(self, alpha=0.15):
        self.alpha = alpha

    def __call__(self, hazards, S, Y, c, alpha=None):
        if alpha is None:
            return nll_loss(hazards, S, Y, c, alpha=self.alpha)
        else:
            return nll_loss(hazards, S, Y, c, alpha=alpha)

def train_survival(args,datasets, cur):
    device=torch.device("cuda" if torch.cuda.is_available() else "cpu")
    train_split, val_split = datasets

    model_dict = {"dim":args.feat_size,"dropout": args.drop_out, 'n_classes': args.n_classes}

    model=model_mil.MILSurv(**model_dict)
    model = model.to(device)
    utils.print_network(model)

    optimizer = utils.get_optim(model, args)
    if args.cos_scheduler:
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer,T_max=args.epochs)
    else:
        scheduler=None

    train_loader = utils.get_split_loader(train_split,n_workers=args.workers, training=True)
    val_loader = utils.get_split_loader(val_split,n_workers=args.workers)

    loss_fn = NLLSurvLoss(args.alpha_surv)

    for epoch in range(args.epochs):
        train_loop_survival(epoch, model, train_loader, optimizer, args.n_classes, loss_fn, args.gc,scheduler)
        validate_survival(cur, epoch, model, val_loader, args.n_classes, loss_fn, args.results_dir)

    torch.save(model.state_dict(), os.path.join(args.results_dir, "s_{}_checkpoint.pt".format(cur)))
    model.load_state_dict(torch.load(os.path.join(args.results_dir, "s_{}_checkpoint.pt".format(cur))))
    results_val_dict, val_cindex = summary_survival(model, val_loader, args.n_classes)
    print('Val c-Index: {:.4f}'.format(val_cindex))
    return results_val_dict, val_cindex


def train_loop_survival(epoch, model, loader, optimizer, n_classes, loss_fn=None, gc=32,scheduler=None):
    device=torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.train()
    train_loss_surv, train_loss = 0., 0.

    print('\n')
    all_risk_scores = np.zeros((len(loader)))
    all_censorships = np.zeros((len(loader)))
    all_event_times = np.zeros((len(loader)))

    if scheduler:
        scheduler.step()


    print(f"LR: {optimizer.param_groups[0]['lr']}")
    for batch_idx, (data_WSI, label, event_time, c) in enumerate(loader):
        #data_WSI=[d.to(device) for d in data_WSI]
        data_WSI = data_WSI.to(device)
        label = label.to(device)
        c = c.to(device)
        hazards, S, Y_hat, _, _ = model(data_WSI) # return hazards, S, Y_hat, A_raw, results_dict
        #print(hazards,S,Y_hat)
        loss = loss_fn(hazards=hazards, S=S, Y=label, c=c)
        loss_value = loss.item()


        risk = -torch.sum(S, dim=1).detach().cpu().numpy()
        all_risk_scores[batch_idx] = risk
        all_censorships[batch_idx] = c.item()
        all_event_times[batch_idx] = event_time

        train_loss_surv += loss_value
        train_loss += loss_value

        if (batch_idx + 1) % 100 == 0:
            print('batch {}, loss: {:.4f}, label: {}, event_time: {:.4f}, risk: {:.4f}, bag_size: {}, censored: {}'.format(batch_idx, loss_value, label.item(), float(event_time), float(risk), data_WSI.size(0),c))
            #print('batch {}, loss: {:.4f}, label: {}, event_time: {:.4f}, risk: {:.4f}, bag_size: {}, censored: {}'.format(batch_idx, loss_value, label.item(), float(event_time), float(risk), data_WSI[0].size(0),c))
            #print(hazards, S, Y_hat)
        # backward pass
        loss = loss / gc
        loss.backward()

        if (batch_idx + 1) % gc == 0:
            optimizer.step()
            optimizer.zero_grad()

    # calculate loss and error for epoch
    train_loss_surv /= len(loader)
    train_loss /= len(loader)

    c_index = concordance_index_censored((1-all_censorships).astype(bool), all_event_times, all_risk_scores, tied_tol=1e-08)[0]

    print('Epoch: {}, train_loss_surv: {:.4f}, train_loss: {:.4f}, train_c_index: {:.4f}'.format(epoch, train_loss_surv, train_loss, c_index))



def validate_survival(cur, epoch, model, loader, n_classes, loss_fn=None, results_dir=None):
    device=torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(device)
    model.eval()
    val_loss_surv, val_loss = 0., 0.
    all_risk_scores = np.zeros((len(loader)))
    all_censorships = np.zeros((len(loader)))
    all_event_times = np.zeros((len(loader)))

    for batch_idx, (data_WSI, label, event_time, c) in enumerate(loader):

        data_WSI = data_WSI.to(device)
        label = label.to(device)
        c = c.to(device)


        with torch.no_grad():
            hazards, S, Y_hat, _, _ = model(data_WSI) # return hazards, S, Y_hat, A_raw, results_dict

        loss = loss_fn(hazards=hazards, S=S, Y=label, c=c, alpha=0)
        loss_value = loss.item()


        risk = -torch.sum(S, dim=1).cpu().numpy()
        all_risk_scores[batch_idx] = risk
        all_censorships[batch_idx] = c.cpu().numpy()
        all_event_times[batch_idx] = event_time

        val_loss_surv += loss_value
        val_loss += loss_value

    val_loss_surv /= len(loader)
    val_loss /= len(loader)
    c_index = concordance_index_censored((1-all_censorships).astype(bool), all_event_times, all_risk_scores, tied_tol=1e-08)[0]
    print(f"VAL LOSS: {val_loss}")
    print(f"VAL C-INDEX: {c_index}")



def summary_survival(model, loader, n_classes):
    device=torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.eval()
    test_loss = 0.

    all_risk_scores = np.zeros((len(loader)))
    all_censorships = np.zeros((len(loader)))
    all_event_times = np.zeros((len(loader)))

    slide_ids = loader.dataset.patients_df['slide_id']
    patient_results = {}

    for batch_idx, (data_WSI, label, event_time, c) in enumerate(loader):

        data_WSI = data_WSI.to(device)
        label = label.to(device)


        slide_id = slide_ids.iloc[batch_idx]

        with torch.no_grad():
            hazards, survival, Y_hat, _, _ = model(data_WSI)

        risk = -torch.sum(survival, dim=1).detach().cpu().numpy()
        event_time = event_time
        c = c.item()
        all_risk_scores[batch_idx] = risk
        all_censorships[batch_idx] = c
        all_event_times[batch_idx] = event_time
        patient_results.update({slide_id: {'slide_id': np.array(slide_id), 'risk': risk, 'disc_label': label.item(), 'survival': event_time, 'censorship': c}})

    c_index = concordance_index_censored((1-all_censorships).astype(bool), all_event_times, all_risk_scores, tied_tol=1e-08)[0]
    return patient_results, c_index
