#!/usr/bin/env python3
#SBATCH -t 71:10:00
#SBATCH -A naiss2025-22-1484
#SBATCH -n 32

import os
import shutil


"""
Author: Najla Hosseini
Email: najla.hosseini@gmail.com| najla.hosseini@ki.se
Affiliation: Department of Oncology-Pathology, Karolinska Institutet,
SE-171 77 Stockholm, Sweden

Description:
This code implements a computational framework combining molecular dynamics,
the Jarzynski estimator, and symmetry-adapted perturbation theory (SAPT)
to analyze amyloid protofibril stability.

Citation:
If you use this code in your research, please cite my article.

"""


base_dir = "/proj/naiss2025-22-1484/users/x_najho/order/nnqq1/control"
n_new = 20 

existing = [d for d in os.listdir('.') if d.startswith("jarzynski-slow_") and os.path.isdir(d)]

if existing:
    max_existing = max(int(d.split('_')[1]) for d in existing)
else:
    max_existing = -1  

start_idx = max_existing + 1

for i in range(start_idx, start_idx + n_new):
    dir_sim = f"jarzynski-slow_{i}"
    os.makedirs(dir_sim, exist_ok=True)
    print(f"Created {dir_sim}")

    for filename in ["nnqq1-control.tpr", "pulling-slow.sh"]:
        shutil.copy(os.path.join(base_dir, filename), dir_sim)

    os.chdir(dir_sim)
    os.system("sbatch pulling-slow.sh")
    os.chdir("..")

