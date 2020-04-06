# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 15:24:02 2020

@author: HP
"""
# this is the inputs file. Major inputs will be written in this file and accessed by all modules
##################################
# Thermodynamic Inputs
T_hot_in = 500 #Celsius
T_cold_in = 142 # Celsius
P_hot_in = 90 #Bar
P_cold_in = 210 #Bar
m_channel = 3e-4# channel mass flow rate, Kg/sec
##################################
# Geometric Inputs
t = 2#mm
w = 1#mm
d= 1#Hydraulic diameter of channel - mm (square channel)
L = 0.5# metre, Length
##################################
# Material inputs
hot_fluid = 'co2'
cold_fluid = 'co2'
k_solid = 385 # w/mk-Copper at 350 Celsius

##################################
# Wall Boundary Conditions - adiabatic or wall heat loss, to be specified as a string

##################################
# Convert all inputs to SI
T_hot_in = T_hot_in + 273.15
T_cold_in = T_cold_in + 273.15
P_hot_in =P_hot_in *1e5
P_cold_in =P_cold_in *1e5
d = d/1000
w = w/1000
t = t/1000
##################################
##################################
# Solver Controls
# To solve 2D TRN (i.e to model one column) set c_ch = 1 and select adiabatic wall BC
a = 2 # number of nodes in x direction(direction of flow)
b_ch = 1 # pair of hot and cold channels along height (Y direction) - 1 pair of channels => b_ch =1
c_ch = 1 # no. of channels along width (Z direction) - one channel means c_ch = 1
b = 4*b_ch + 1# number of nodes in y direction
c = 2*c_ch + 1# number of nodes in z direction
delta_x = L/a # length along x of each CV
##################################
