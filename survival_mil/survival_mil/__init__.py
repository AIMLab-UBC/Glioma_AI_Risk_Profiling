import argparse
import os
import random
from timeit import default_timer as timer

import numpy as np
import pandas as pd
import torch

from . import trainer
from .datasets_survival import Generic_MIL_Survival_Dataset


parser = argparse.ArgumentParser(description="MIL Survival Analysis")
### Checkpoint + Misc. Pathing Parameters
parser.add_argument(
    "--data_root_dir", type=str, default="/media/ssd1/pan-cancer/", help="data directory"
)
parser.add_argument(
    "--csv_path", type=str, default="./csv_files/surv.csv", help="csv file w/ survival data"
)
parser.add_argument(
    "--label_col", type=str, default="survival_months", help="name of column w/ survival time"
)
parser.add_argument(
    "--seed", type=int, default=1, help="Random seed for reproducible experiment (default: 1)"
)
parser.add_argument("--k", type=int, default=5, help="Number of folds (default: 5)")
parser.add_argument(
    "--results_dir", type=str, default="./results", help="Results directory (Default: ./results)"
)
parser.add_argument(
    "--split_dir",
    type=str,
    default="tcga_blca_100",
    help="Which cancer type within ./splits/<which_splits> to use for training.",
)

### Model Parameters.
parser.add_argument("--drop_out", type=float, default=0.0, help="Dropping out random patches.")
parser.add_argument("--feat_size", type=int, default=512, help="input feat size")
parser.add_argument("--opt", type=str, choices=["adam", "sgd", "nadam", "adamw"], default="adam")
parser.add_argument("--gc", type=int, default=32, help="Gradient Accumulation Step.")
parser.add_argument("--workers", type=int, default=4, help="num workers for dataloader")
parser.add_argument(
    "--epochs", type=int, default=20, help="number of epochs to train (default: 20)"
)
parser.add_argument(
    "--n_classes", type=int, default=4, help="number of bins/classes"
)
parser.add_argument("--lr", type=float, default=2e-4, help="learning rate")

parser.add_argument(
    "--reg", type=float, default=1e-5, help="L2-regularization weight decay (default: 1e-5)"
)
parser.add_argument(
    "--alpha_surv", type=float, default=0.0, help="How much to weigh censored patients"
)
parser.add_argument(
    "--cos_scheduler", action="store_true", default=False, help="Enable cos lr scheduler"
)


def seed_torch(device, seed=7):
    random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if device.type == "cuda":
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)  # if you are using multi-GPU.
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True


def main_run(args, dataset):

    if not os.path.isdir(args.results_dir):
        os.mkdir(args.results_dir)

    latest_val_cindex = []
    folds = np.arange(0, args.k)

    for i in folds:
        print(f"training fold {i+1}")
        train_dataset, val_dataset = dataset.return_splits(
            csv_path=f"{args.split_dir}/splits_{i}.csv"
        )

        print(f"training: {len(train_dataset)}, validation: {len(val_dataset)}")
        datasets = (train_dataset, val_dataset)

        val_latest, cindex_latest = trainer.train_survival(args,datasets, i)
        latest_val_cindex.append(cindex_latest)

    results_latest_df = pd.DataFrame({"folds": folds, "val_cindex": latest_val_cindex})
    results_latest_df.to_csv(os.path.join(args.results_dir, "summary.csv"))


def run_pipeline():
    args = parser.parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    seed_torch(device,args.seed)
    print(torch.__version__)
    print(np.version.version)

    settings = {
        "num_splits": args.k,
        "epochs": args.epochs,
        "results_dir": args.results_dir,
        "lr": args.lr,
        "seed": args.seed,
        "gc": args.gc,
        "opt": args.opt,
    }

    dataset = Generic_MIL_Survival_Dataset(
        csv_path=args.csv_path,
        data_dir=args.data_root_dir,
        seed=args.seed,
        n_bins=args.n_classes,
        label_col=args.label_col,
        ignore=[],
    )
    args.results_dir = os.path.join(args.results_dir, f"amil_s{args.seed}")
    if not os.path.isdir(args.results_dir):
        os.makedirs(args.results_dir)


    print("split_dir", args.split_dir)
    assert os.path.isdir(args.split_dir)
    settings.update({"split_dir": args.split_dir})

    with open(args.results_dir + "/experiment.txt", "w") as f:
        print(settings, file=f)
    f.close()
    print("################# Settings ###################")
    for key, val in settings.items():
        print(f"{key}:  {val}")

    current_directory = os.getcwd()
    print("Current Directory:", current_directory)
    start = timer()
    main_run(args, dataset)
    end = timer()
    print("finished!")
    print("end script")
    print("Script Time: %f seconds" % (end - start))
