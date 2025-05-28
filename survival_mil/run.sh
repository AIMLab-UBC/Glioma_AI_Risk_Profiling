#!/bin/bash
python main.py \
--data_root_dir="/feats/" \
--csv_path="survival.csv" \
--label_col="survival_months" \
--seed=7 \
--opt="adamw" \
--results_dir="results/" \
--split_dir="/splits/" \
--k=5 \
--n_classes=4 \
--drop_out=0.25 \
--lr=0.0001 \
--epochs=1 \
--gc=32 \
--feat_size=1536
