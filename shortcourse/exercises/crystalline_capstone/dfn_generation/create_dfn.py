#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DFN generating script for PFLOTRAN Advanced Short Course Practice Problem
Summer 2023Created on Thu Sep 10 12:04:03 2020
"""

import os, sys
from pydfnworks import * 

src_path = os.getcwd()

jobname = f"{src_path}/output"   
        
#dfnFlow_file = f"{src_path}/dfn_explicit.in"
DFN = DFNWORKS(jobname)
#               dfnFlow_file=dfnFlow_file)

DFN.params['domainSize']['value'] = ...
DFN.params['h']['value'] = 1.0
DFN.params['seed']['value'] = 20201153 #seed for random generator 0 seeds off clock

DFN.add_fracture_family(shape=...,
                        distribution=...,
                        p32 = ...,
                        beta_distribution = 1,
                        beta = 0,
                        number_of_points = 12,
                        kappa=...,
                        theta=...,
                        phi=...,
                        alpha=...,
                        min_radius=...,
                        max_radius=...,
                        hy_variable='permeability',
                        hy_function=...,
                        hy_params={
                            "alpha": ...,
                            "beta": ...,
                            "sigma": ...
                        })


DFN.add_fracture_family(shape=...,
                        distribution=...,
                        p32 = ...,
                        beta_distribution = 1,
                        beta = 0,
                        number_of_points = 12,
                        kappa=...,
                        theta=...,
                        phi=...,
                        alpha=...,
                        min_radius=...,
                        max_radius=...,
                        hy_variable='permeability',
                        hy_function=...,
                        hy_params={
                            "alpha": ...,
                            "beta": ...,
                            "sigma": ...
                        })

DFN.add_fracture_family(shape=...,
                        distribution=...,
                        p32 = ...,
                        beta_distribution = 1,
                        beta = 0,
                        number_of_points = 12,
                        kappa=...,
                        theta=...,
                        phi=...,
                        alpha=...,
                        min_radius=...,
                        max_radius=...,
                        hy_variable='permeability',
                        hy_function=...,
                        hy_params={
                            "alpha": ...,
                            "beta": ...,
                            "sigma": ...
                        })

DFN.add_user_fract(shape=...,
                    radii=...,
                    translation=...,
                    aspect_ratio=1,
                    beta=0,
                    normal_vector=...,
                    number_of_vertices=5,
                    aperture=...)

DFN.add_user_fract(shape=...,
                    radii=...,
                    translation=...,
                    aspect_ratio=1,
                    beta=0,
                    normal_vector=...,
                    number_of_vertices=5,
                    aperture=...)

DFN.make_working_directory(delete=True)
DFN.print_domain_parameters()
DFN.check_input()
DFN.create_network()
DFN.dump_hydraulic_values()


