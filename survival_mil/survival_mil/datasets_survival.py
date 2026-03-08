from __future__ import print_function, division
import os
import random
import glob

import h5py
import numpy as np
import pandas as pd

import torch
from torch.utils.data import Dataset


class Generic_WSI_Survival_Dataset(Dataset):
    #code adapted from https://github.com/mahmoodlab/MCAT/blob/master/datasets/dataset_survival.py
    def __init__(self,
        csv_path = 'dataset_csv/ccrcc_clean.csv', seed = 7, n_bins = 4, ignore=[], label_col = None, eps=1e-6):

        self.seed = seed
        slide_data = pd.read_csv(csv_path, low_memory=False)
        self.label_col = label_col


        patients_df = slide_data.drop_duplicates(['case_id']).copy()
        uncensored_df = patients_df[patients_df['censorship'] < 1]
        disc_labels, q_bins = pd.qcut(uncensored_df[label_col], q=n_bins, retbins=True, labels=False)
        q_bins[-1] = slide_data[label_col].max() + eps
        q_bins[0] = slide_data[label_col].min() - eps

        disc_labels, q_bins = pd.cut(patients_df[label_col], bins=q_bins, retbins=True, labels=False, right=False, include_lowest=True)
        patients_df.insert(2, 'label', disc_labels.values.astype(int))


        patient_dict = {}
        slide_data = slide_data.set_index('case_id')

        for patient in patients_df['case_id']:
            slide_ids = slide_data.loc[patient, 'slide_id']
            if isinstance(slide_ids, str):
                slide_ids = np.array(slide_ids).reshape(-1)
            else:
                slide_ids = slide_ids.values
            patient_dict.update({patient:slide_ids})

        self.patient_dict = patient_dict

        #slide_data = patients_df
        patients_df.reset_index(drop=True, inplace=True)
        patients_df = patients_df.assign(slide_id=patients_df['case_id'])


        label_dict = {}
        key_count = 0
        for i in range(len(q_bins)-1):
            for c in [0, 1]:
                print('{} : {}'.format((i, c), key_count))
                label_dict.update({(i, c):key_count})
                key_count+=1


        self.label_dict = label_dict
        for i in patients_df.index:
            key = patients_df.loc[i, 'label']
            patients_df.at[i, 'disc_label'] = key
            censorship = patients_df.loc[i, 'censorship']
            key = (key, int(censorship))
            patients_df.at[i, 'label'] = label_dict[key]

        self.bins = q_bins
        self.num_classes=len(self.label_dict)

        new_cols = list(patients_df.columns[-2:]) + list(patients_df.columns[:-2])
        patients_df = patients_df[new_cols]
        self.patients_df = patients_df
        self.cls_ids_prep()

        self.summarize()


    def cls_ids_prep(self):

        self.cls_ids = [[] for i in range(self.num_classes)]
        for i in range(self.num_classes):
            self.cls_ids[i] = np.where(self.patients_df['label'] == i)[0]



    @staticmethod
    def df_prep(data, n_bins, ignore, label_col):
        mask = data[label_col].isin(ignore)
        data = data[~mask]
        data.reset_index(drop=True, inplace=True)
        disc_labels, bins = pd.cut(data[label_col], bins=n_bins)
        return data, bins

    def __len__(self):
        return len(self.patients_df)

    def summarize(self):
        print("label column: {}".format(self.label_col))
        print("label dictionary: {}".format(self.label_dict))
        print("number of classes: {}".format(self.num_classes))
        print("slide-level counts: ", '\n', self.patients_df['label'].value_counts(sort = False))
        for i in range(self.num_classes):
            print('Patient-LVL; Number of samples registered in class %d: %d' % (i, self.cls_ids[i].shape[0]))


    def get_split_from_df(self, all_splits: dict, split_key: str='train', scaler=None):
        split = all_splits[split_key]
        split = split.dropna().reset_index(drop=True)

        if len(split) > 0:
            mask = self.patients_df['slide_id'].isin(split.tolist())
            df_slice = self.patients_df[mask].reset_index(drop=True)
            split = Generic_Split(df_slice,  data_dir=self.data_dir, label_col=self.label_col, patient_dict=self.patient_dict, num_classes=self.num_classes, label_dict=self.label_dict)
        else:
            split = None

        return split


    def return_splits(self, csv_path):
        assert csv_path
        all_splits = pd.read_csv(csv_path)
        train_split = self.get_split_from_df(all_splits=all_splits, split_key='train')
        val_split = self.get_split_from_df(all_splits=all_splits, split_key='val')

        return train_split, val_split


    def get_list(self, ids):
        return self.patients_df['slide_id'][ids]

    def getlabel(self, ids):
        return self.patients_df['label'][ids]

    def __getitem__(self, idx):
        return None

class Generic_MIL_Survival_Dataset(Generic_WSI_Survival_Dataset):
    def __init__(self, data_dir, **kwargs):
        super(Generic_MIL_Survival_Dataset, self).__init__(**kwargs)
        self.data_dir = data_dir

    def __getitem__(self, idx):
        case_id = self.patients_df['case_id'][idx]
        label = self.patients_df['disc_label'][idx]
        event_time = self.patients_df[self.label_col][idx]
        c = self.patients_df['censorship'][idx]
        slide_ids = self.patient_dict[case_id]


        if type(self.data_dir) == dict:
            source = self.patients_df['oncotree_code'][idx]
            data_dir = self.data_dir[source]
        else:
            data_dir = self.data_dir


        if self.data_dir:
            path_features = []

            for i,slide_id in enumerate(slide_ids):
                wsi_path = os.path.join(data_dir, '{}.h5'.format(slide_id.rstrip('.svs')))
                #print(wsi_path)
                with h5py.File(wsi_path,'r') as hdf5_file:
                    features = hdf5_file['features'][:]
                #print(features.shape)
                path_features.append(torch.tensor(features))
            path_features=random.choice(path_features)
            return (path_features, label, event_time, c)



class Generic_Split(Generic_MIL_Survival_Dataset):
    def __init__(self, patients_df, data_dir=None, label_col=None, patient_dict=None, num_classes=2, label_dict={}):
        self.use_h5 = False
        self.patients_df = patients_df
        self.data_dir = data_dir
        self.num_classes = num_classes
        self.label_col = label_col
        self.patient_dict = patient_dict
        self.cls_ids = [[] for i in range(self.num_classes)]
        self.label_dict=label_dict
        for i in range(self.num_classes):
            self.cls_ids[i] = np.where(self.patients_df['label'] == i)[0]

        patients_df = self.patients_df.drop_duplicates(['case_id']).copy()


    def __len__(self):
        return len(self.patients_df)
