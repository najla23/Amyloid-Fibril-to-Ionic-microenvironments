#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
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



def trapezoidal_int1(x, f):
    h = x[1] - x[0]
    cumulative_int = np.zeros_like(x)
    
    for i in range(1, len(x)):
        cumulative_int[i] = cumulative_int[i-1] + (h/2) * (f[i-1] + f[i])
    
    return cumulative_int


def trapezoidal_int(x, f):
    h = x[1] - x[0]
    cumulative_int = np.zeros_like(x)
    
    for i in range(1, len(x)):
        cumulative_int[i] =  (h/2) * (f[i-1] + f[i])
    
    return cumulative_int


k=0.008314
T=298.15
all_force_data_smoothed = []
all_force_data_smoothed2 = []
all_distance_data_smoothed = []
all_distance_data_smoothed2 = []
w1=[]
w444=[]
xx=[]
xx=[]
xxx44=[]
xxx4=[]

i=1
plt.figure(figsize=(8, 6))


color_map = cm.Purples(np.linspace(0.4, 0.9, len(force_files))) 


for (force_file, distance_file), color in zip(zip(force_files, distance_files), color_map):

    force_data = np.loadtxt(force_file, skiprows=17, usecols=(1), max_rows=79000)
    distance_data = np.loadtxt(distance_file, skiprows=17, usecols=(1), max_rows=79000)
    sorted_indices = np.argsort(distance_data)


    sorted_distance = distance_data[sorted_indices]
    sorted_force = force_data[sorted_indices]
    

    window_size = 5000
    force_data_smoothed = moving_average(sorted_force, window_size)
    distance_data_smoothed = moving_average(sorted_distance, window_size)

    all_force_data_smoothed.append(sorted_force)
    all_distance_data_smoothed.append(sorted_distance)
    
    k=0.008314
    T=298.15

    w44=trapezoidal_int(sorted_distance,sorted_force)
    print(i)
    print(w44)

    window_size = 100
    w44_a = moving_average(w44, window_size)
    distance_data_smoothed = moving_average(sorted_distance, window_size)
    max_work = np.max(w44_a)

    print(f"Trajectory {i}: max work = {max_work:.3f} (from {force_file})")


    plt.plot(distance_data_smoothed, w44_a, lw=2.5, color=color, alpha=0.9)


    plt.xlabel('Distance (nm)', fontsize=16)
    plt.ylabel('Work (kJ/mol)', fontsize=16)
    plt.savefig('work.pdf', dpi=100, bbox_inches='tight')
    plt.xlim(1.9,6)
    w66=np.exp(-w44/(k*T))
    w444.append(w66)
    xx.append(sorted_distance)
    i+=1

x=np.mean(np.array(xx), axis=0)
xxx44=np.mean(np.array(w444), axis=0)
xxx444 = -np.log(xxx44)

for i in range (len(xxx444)):
     xxx4444=xxx444*(k*T)

xxx44442=xxx4444*k*T*2.48 
window_size = 9000
force_data_smoothed23 = moving_average(xxx44442, window_size)
distance_data_smoothed23 = moving_average(x, window_size)

plt.figure(figsize=(8, 6))


orange_map = cm.Oranges(np.linspace(0.5, 0.9, 1))
plt.plot(distance_data_smoothed23, force_data_smoothed23, 
         color=orange_map[0], lw=2.5, alpha=0.95)


plt.xlabel('Distance (nm)', fontsize=16)
plt.ylabel('Free energy (kJ/mol)', fontsize=16)
plt.xlim(0,6)
plt.savefig('freeenergy.pdf', dpi=100, bbox_inches='tight')
plt.show()

