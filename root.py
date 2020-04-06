# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 15:46:14 2020

@author: HP
"""
import numpy as np
import matplotlib.pyplot as plt 
from CoolProp.CoolProp import PropsSI
import inputs as I
import initialize
import find
import update


Tg,Txy,Tyz,Tzx = initialize.temperature()
# print("Tg")
# print(".................")
# print(Tg[:,:,1]- 273.15)

# print(".................")
# print("Tyz")
# print(".................")
# print(Tyz[:,:,1]- 273.15)

Pg, P = initialize.pressure()
# print("Pg")
# print(".................")
# print(Pg[:,:,1]/1e5)

# print(".................")
# print("P")
# print(".................")
# print(P[:,:,1]/1e5)

Pg, P = update.pressure(P,Pg,Tg)
# print("Pg")
# print(".................")
# print(Pg[:,:,1]/1e5)

# print(".................")
# print("P")
# print(".................")
# print(P[:,:,1]/1e5)

R = update.resistance(P,Tyz,Pg,Tg)

# print(R[0,1,1,:,:])
# print(R[3,1,1,:,:])
# print(R[7,1,1,:,:])
# print("...........")
# print(R[0,3,1,:,:])
# print(R[3,3,1,:,:])
# print(R[7,3,1,:,:])
A,C = update.linear_equation_system(R,Tg,Tyz,Pg,P)
sol = np.linalg.inv(A)

print(sol)
# target = np.matmul(sol,C)

# print(target)
# print(np.size(A))
# print(np.size(C))

           
       