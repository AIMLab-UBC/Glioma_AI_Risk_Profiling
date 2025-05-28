import glob
import torch
import h5py
import csv
import os

from models.amil import *


feat_files=glob.glob("/PATHTOFILES/*.h5")
device=torch.device("cuda" if torch.cuda.is_available() else "cpu")

model_dict = {"dim":1536,"dropout": 0.1, 'n_classes': 4}

model=ABMILSurv(**model_dict)
hazards=[]
for i in range(10):
	print(f"FOLD: ", i)
	ckpt=f"/PATHTOCHECKPOINTS/s_{i}_checkpoint.pt"
	ckpt=torch.load(ckpt,map_location=device)
	model.load_state_dict(ckpt)
	model.eval()
	hazard =[]
	for f in feat_files:
		with h5py.File(f,'r') as hdf5_file:
			features = hdf5_file['embeddings'][:]
		out=model(torch.tensor(features))
		risk = -torch.sum(out[1], dim=1).detach().cpu().numpy()
		hazard.append(risk[0])

	hazards.append(hazard)

with open(f'results_amil.csv', mode='w', newline='') as file:
	writer = csv.writer(file)
	writer.writerow(['Name', 'Hazard1', 'Hazard2', 'Hazard3', 'Hazard4', 'Hazard5'])
	for i,f in enumerate(feat_files):
		name=os.path.basename(f).split('.')[0]
		hazard_values = [hazards[j][i] for j in range(10)]  # Collect hazards for this file across all 5 folds
		writer.writerow([name] + hazard_values)
