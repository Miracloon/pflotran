#!/bin/python
#%%
import os
import matplotlib.pyplot as plt
import numpy as np
import h5py

if not os.path.exists('figures'): os.mkdir('figures')

#%%
# Setup and import
h5f = h5py.File('dakota_results.h5')
parameter_names = list(h5f['_scales/models/simulation/NO_MODEL_ID/variables/continuous_descriptors'].asstr()[...])
response_names = list(h5f['_scales/models/simulation/NO_MODEL_ID/responses/function_descriptors'].asstr()[...])

in_samples = h5f['methods/NO_METHOD_ID/sources/NO_MODEL_ID/variables/continuous'][...]
out_samples = h5f['methods/NO_METHOD_ID/sources/NO_MODEL_ID/responses/functions'][...]

ccs = h5f['methods/NO_METHOD_ID/results/execution:1/simple_correlations'][...]
rank_ccs = h5f['methods/NO_METHOD_ID/results/execution:1/simple_rank_correlations'][...]

n_params = len(parameter_names)
n_responses = len(response_names)
n_samples = in_samples.shape[0]

#%%
# Convenience functions

def smart_sci(x, precision=1):
    s = f"{x:.{precision}e}"
    mantissa, exp = s.split("e")
    exp = int(exp)

    if exp == 0:
        return mantissa
    return f"{mantissa}e{exp:+d}"

def scientific_format_ticklabel(axis):
  ticks = axis.get_ticklocs()[1:-1]
  axis.set_ticks(ticks)
  axis.set_ticklabels( [ smart_sci(t) for t in ticks ] )

#%%
# Histograms
for out_sample, out_label in zip(out_samples.T, response_names):
  fig, ax = plt.subplots(figsize=(4,3))
  ax.hist(out_sample, histtype='step', density=True)
  ax.set_xlabel(out_label)
  ax.set_ylabel('Histogram')
  fig.tight_layout()
  fig.savefig(f'figures/{out_label}_histogram.png', dpi=300)

for in_sample, in_label in zip(in_samples.T, parameter_names):
  fig, ax = plt.subplots(figsize=(4,3))
  ax.hist(in_sample, histtype='step', density=True)
  ax.set_xlabel(in_label)
  ax.set_ylabel('Histogram')
  fig.tight_layout()
  fig.savefig(f'figures/{in_label}_histogram.png', dpi=300)

#%%
# Scatterplots

width=10
aspect_ratio = n_responses/n_params*1.3
fig, axs = plt.subplots(n_responses, n_params, figsize=(width,width*aspect_ratio), sharey=True)
axs = np.atleast_2d(axs)
for i, (out_sample, out_label) in enumerate(zip(out_samples.T, response_names)):
  for j, (in_sample, in_label) in enumerate(zip(in_samples.T, parameter_names)):
      axs[i,j].plot( in_sample, out_sample, '.', ms=4)
      if i==n_responses-1:
        axs[i,j].set_xlabel(in_label)
      if j==0:
        axs[i,j].set_ylabel(out_label)
      axs[i,j].set_box_aspect(1)
      if in_sample.max() < 1e-5:
        scientific_format_ticklabel(axs[i,j].xaxis)

      if out_sample.max() < 1e-5:
        scientific_format_ticklabel(axs[i,j].yaxis)
fig.tight_layout()

#%%
# Correlations
def plot_ccs(ccs, title, filename):
  fig, ax = plt.subplots(figsize=(6,6))
  cm = ax.matshow(ccs, vmin=-1,vmax=1, cmap='seismic')
  fig.colorbar(cm, shrink=0.65)

  for i in range(ccs.shape[0]):
      for j in range(ccs.shape[1]):
          val = ccs[i, j]
          ax.text( j, i, f"{val:.2f}",
              ha="center", va="center", color="white" if abs(val) > 0.4 else "black")
  ax.set_xticks(range(n_params+n_responses));
  ax.set_yticks(range(n_params+n_responses));
  ax.set_xticklabels(parameter_names+response_names, rotation=45, ha='left')
  ax.set_yticklabels(parameter_names+response_names)
  ax.set_title(title)
  fig.tight_layout()
  fig.savefig(f'figures/{filename}', dpi=300)

plot_ccs(ccs, 'Simple Correlation Coefficients', 'simple_ccs.png')
plot_ccs(rank_ccs, 'Rank Correlation Coefficients', 'simple_ccs.png')

#%%
# ECDFs

for out_sample, out_label in zip(out_samples.T, response_names):
  fig, ax = plt.subplots(figsize=(4,3))
  ax.ecdf(out_sample)
  ax.set_xlabel(out_label)
  ax.set_ylabel('Empirical CDF')
  fig.tight_layout()
  fig.savefig(f'figures/{out_label}_cdf.png', dpi=300)


#%%
# Plotting time history samples

# Which locations to plot
#plot_location = 'repo'
plot_location = 'aquifer'

# Read all of the PFLOTRAN output:

time = []
I129 = []
tracer1 = []
tracer2 = []
tracer3 = []
tracer4 = []
print('Reading PFLOTRAN output. . . ')
for i in range(n_samples):
  I129.append([])
  tracer1.append([])
  tracer2.append([])
  tracer3.append([])
  tracer4.append([])
  time.append([])
  os.chdir('samples/workdir.' + str(i+1))
  f = open('pflotran-obs-0.pft','r')
  k = 0
  for line in f:
    k = k + 1
    if k > 1:  # There is 1 header line we don't want to read
      if plot_location == 'aquifer':
        time[i].append(float(line.split()[0]))
        I129[i].append(float(line.split()[6]))
        tracer1[i].append(float(line.split()[7]))
        tracer2[i].append(float(line.split()[8]))
        tracer3[i].append(float(line.split()[9]))
        tracer4[i].append(float(line.split()[10]))
      elif plot_location == 'repo':
        time[i].append(float(line.split()[0]))
        I129[i].append(float(line.split()[22]))
        tracer1[i].append(float(line.split()[23]))
        tracer2[i].append(float(line.split()[24]))
        tracer3[i].append(float(line.split()[25]))
        tracer4[i].append(float(line.split()[26]))
  f.close()
  os.chdir('../../')

output_name = 'I-129 Concentration at Observation Point'
fig, ax = plt.subplots(1,1,tight_layout=True,figsize=(4.0, 5.0),squeeze=False)
for i in range(n_samples):
  ax[0,0].loglog(time[i],I129[i])
#ax[0,0].plot((LowerCI_Mean, UpperCI_Mean), (100.0, 100.0), 'k-')
ax[0,0].set_ylabel('I-129 Concentration [M]') 
ax[0,0].set_xlabel(output_name)
ax[0,0].set_title(str(n_samples) + ' PFLOTRAN simulations')
fig.savefig('figures/i129_conc_' + plot_location + '.png', dpi=300)
plt.close(fig)

output_name = 'Tracer 1 Concentration at Observation Point'
fig, ax = plt.subplots(1,1,tight_layout=True,figsize=(4.0, 5.0),squeeze=False)
for i in range(n_samples):
  ax[0,0].loglog(time[i],tracer1[i])
#ax[0,0].plot((LowerCI_Mean, UpperCI_Mean), (100.0, 100.0), 'k-')
ax[0,0].set_ylabel('Tracer1 Concentration [M]')
ax[0,0].set_xlabel(output_name)
ax[0,0].set_title(str(n_samples) + ' PFLOTRAN simulations')
fig.savefig('figures/tracer1_conc_' + plot_location + '.png', dpi=300)
plt.close(fig)

output_name = 'Tracer 2 Concentration at Observation Point'
fig, ax = plt.subplots(1,1,tight_layout=True,figsize=(4.0, 5.0),squeeze=False)
for i in range(n_samples):
  ax[0,0].loglog(time[i],tracer2[i])
#ax[0,0].plot((LowerCI_Mean, UpperCI_Mean), (100.0, 100.0), 'k-')
ax[0,0].set_ylabel('Tracer2 Concentration [M]')
ax[0,0].set_xlabel(output_name)
ax[0,0].set_title(str(n_samples) + ' PFLOTRAN simulations')
fig.savefig('figures/tracer2_conc_' + plot_location + '.png', dpi=300)
plt.close(fig)

output_name = 'Tracer 3 Concentration at Observation Point'
fig, ax = plt.subplots(1,1,tight_layout=True,figsize=(4.0, 5.0),squeeze=False)
for i in range(n_samples):
  ax[0,0].loglog(time[i],tracer3[i])
#ax[0,0].plot((LowerCI_Mean, UpperCI_Mean), (100.0, 100.0), 'k-')
ax[0,0].set_ylabel('Tracer3 Concentration [M]')
ax[0,0].set_xlabel(output_name)
ax[0,0].set_title(str(n_samples) + ' PFLOTRAN simulations')
fig.savefig('figures/tracer3_conc_' + plot_location + '.png', dpi=300)
plt.close(fig)

output_name = 'Tracer 4 Concentration at Observation Point'
fig, ax = plt.subplots(1,1,tight_layout=True,figsize=(4.0, 5.0),squeeze=False)
for i in range(n_samples):
  ax[0,0].loglog(time[i],tracer4[i])
#ax[0,0].plot((LowerCI_Mean, UpperCI_Mean), (100.0, 100.0), 'k-')
ax[0,0].set_ylabel('Tracer4 Concentration [M]')
ax[0,0].set_xlabel(output_name)
ax[0,0].set_title(str(n_samples) + ' PFLOTRAN simulations')
fig.savefig('figures/tracer4_conc_' + plot_location + '.png', dpi=300)
plt.close(fig)