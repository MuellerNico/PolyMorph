import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
import pandas as pd
from math import sqrt
from scipy import stats

pos_offset = 19 # source position (data is in absolute coordinates)
n = 100 # number of samples
def stderr(stddev):
    return stddev / sqrt(2*(n-1)) # standard error of the standard deviation

viridis = colormaps.get_cmap('viridis')  # Get the viridis colormap
#plt.rcParams['text.usetex'] = True

# --- cv ---
# positional error
data = pd.read_csv('data/positional_error_cv.csv')
data['readout_pos'] = data['readout_pos'] + pos_offset
grad_cv = data['cv'].unique()
thresholds = data['threshold'].unique()
pos_err = data.groupby(['threshold', 'cv'])['readout_pos'].std().reset_index()
colors = viridis(np.linspace(0, 1, len(thresholds)+1))

plt.figure(figsize=(10, 5))
for i, thresh in enumerate(thresholds): 
    subset = pos_err[pos_err['threshold'] == thresh]
    plt.plot(subset['cv'], subset['readout_pos'], marker='o', label=f'Threshold C_theta = {thresh}', color=colors[i])
    plt.errorbar(subset['cv'], subset['readout_pos'], yerr=stderr(subset['readout_pos']), color=colors[i])

plt.title('Positional error vs gradient variability')
plt.xscale('log')
plt.xlabel('Gradient CV (D,k,p)')
plt.ylabel('Positional error')
plt.grid(True)
plt.legend()
plt.savefig("positional_error_cv.pdf")

# readout position
plt.figure(figsize=(10, 5))
for i, thresh in enumerate(thresholds):
    subset = data[data['threshold'] == thresh]
    readout_pos = subset.groupby('cv')['readout_pos'].mean()
    min = subset.groupby('cv')['readout_pos'].min()
    max = subset.groupby('cv')['readout_pos'].max()
    plt.plot(grad_cv, readout_pos, marker='o', label='Threshold C_theta = ' + f'{thresh}', color=colors[i])
    plt.fill_between(grad_cv, min, max, alpha=0.2, color=colors[i]) # min and max values as error estimates

plt.title('Readout position vs gradient variability')
plt.xscale('log')
plt.xlabel('Gradient CV (D,k,p)')
plt.ylabel('Readout Position / avg. cell radius')
plt.grid(True)
plt.legend()
plt.savefig("readout_position.pdf")

    
# --- width ---
# positional error
data = pd.read_csv('data/positional_error_width.csv')
data['readout_pos'] = data['readout_pos'] + pos_offset
widths = data['width'].unique()
thresholds = data['threshold'].unique()
pos_err = data.groupby(['threshold', 'width'])['readout_pos'].std().reset_index()
colors = viridis(np.linspace(0, 1, len(thresholds+2)))

plt.figure(figsize=(10, 5))
for i, thresh in enumerate(thresholds): 
    subset = pos_err[pos_err['threshold'] == thresh]
    plt.plot(subset['width'], subset['readout_pos'], marker='o', label=f'Threshold C_theta = {thresh}', color=colors[i])
    plt.errorbar(subset['width'], subset['readout_pos'], yerr=stderr(subset['readout_pos']), color=colors[i])

# plot reference line 1/sqrt(width)
x = np.linspace(widths.min(), widths.max(), 100)
y = 2 / np.sqrt(x)
plt.plot(x, y, label='2.0/sqrt(width)', color='gray', linestyle='--')

plt.title('Positional error vs domain width')
plt.xlabel('Domain width / avg. cell radius]')
plt.ylabel('Positional error')
plt.grid(True)
plt.legend()
plt.savefig("positional_error_width.pdf")