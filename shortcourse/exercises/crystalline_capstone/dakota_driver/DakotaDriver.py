#!/usr/bin/env python3

#############
#  Imports  #
#############

# Common imports
import subprocess
import pandas as pd
import sys
import dakota.interfacing as di

###################
#  Preprocessing  #
###################

# Get the Dakota Parameters and Response objects.
params, results = di.read_parameters_file()

# Substitute the parameter values into the templatized
# input file, outputting a valid input file for the 
# black-box code.
di.dprepro(template='pflotran.template', 
           parameters=params, 
           output='pflotran.in')

##################################
#  Running black-box simulation  #
##################################

# For more complex calls (e.g., PFLOTRAN in parallel) 
# this may instead invoke a run script.
subprocess.run('pflotran')

####################
#  Postprocessing  #
####################
 
# Importing results from PFLOTRAN outputs

# Get the names of columns in the file
with open( 'pflotran-obs-1-agg-1.pft', 'r') as f: 
    cols = f.readlines()[1].strip().replace('"','').split(',')

# Read in the file
output_data = pd.read_csv('pflotran-obs-1-agg-1.pft', skiprows=2, sep=r'\s+', names=cols )

# Check that the final timestep was reached, else error. 
if not float(output_data['Time [y]'].iloc[-1]) == 1e6:
    sys.exit(f'ERROR: PFLOTRAN did not complete for eval id {params.eval_id}.')

# 'Total I129 [M]' reports the max I129 concentration in the aquifer region 
# at each timestep. To get the overall peak we take the max.
peak_i129 = output_data['Total I129 [M]'].max()

# Writing results to Dakota-formatted results file using the Response object.
results['Peak_I129_in_Aquifer_M'].function = peak_i129
results.write()