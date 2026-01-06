#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import glob
import re


def moving_average(data, window_size):
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')


force_pattern = 'slow_pullf*.xvg'
distance_pattern = 'slow_pullx*.xvg'


def get_index(fname):
    m = re.search(r'(\d+)', fname)
    return int(m.group(1)) if m else -1

force_files = sorted(glob.glob(force_pattern), key=get_index)
distance_files = sorted(glob.glob(distance_pattern), key=get_index)


plt.figure(figsize=(8, 6))


cmap = plt.get_cmap('Purples')
colors = [cmap(0.35 + 0.55 * (i / max(1, len(force_files)-1))) for i in range(len(force_files))]


all_forces = []
all_distances = []


for i, (force_file, distance_file) in enumerate(zip(force_files, distance_files)):
    force_data = np.loadtxt(force_file, skiprows=17, usecols=(1), max_rows=79000)
    distance_data = np.loadtxt(distance_file, skiprows=17, usecols=(1), max_rows=79000)


    sorted_indices = np.argsort(distance_data)
    distance_data = distance_data[sorted_indices]
    force_data = force_data[sorted_indices]


    window_size = 3000
    force_data_smoothed = moving_average(force_data, window_size)
    distance_data_smoothed = moving_average(distance_data, window_size)


    all_forces.append(force_data_smoothed)
    all_distances.append(distance_data_smoothed)

    print(f"{i}: force_file = '{force_file}', distance_file = '{distance_file}'")


    plt.plot(distance_data_smoothed, force_data_smoothed,
             color=colors[i],
             lw=2.2,
             alpha=0.9)


min_len = min(map(len, all_forces))
forces_trunc = np.array([f[:min_len] for f in all_forces])
distances_trunc = np.array([d[:min_len] for d in all_distances])

mean_force = np.mean(forces_trunc, axis=0)
std_force = np.std(forces_trunc, axis=0)
mean_distance = np.mean(distances_trunc, axis=0)


plt.fill_between(mean_distance,
                 mean_force - std_force,
                 mean_force + std_force,
                 color='#6A5ACD',  
                 alpha=0.18)
plt.plot(mean_distance, mean_force,
         color='red', lw=3, alpha=0.95)


plt.xlabel(r'Distance (nm)', fontsize=24)
plt.ylabel(r'Force (kJ/mol/nm)', fontsize=24)
plt.xlim(0, 6)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()


plt.savefig('Force_vs_Distance_PurpleGradient.pdf', dpi=600, bbox_inches='tight')
plt.show()

