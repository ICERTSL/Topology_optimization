# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 15:46:55 2020

@author: HP
"""
import numpy as np
import matplotlib.pyplot as plt 
from CoolProp.CoolProp import PropsSI
import inputs as I
import find

def temperature():
    a = I.a
    b = I.b
    c = I.c
    Tg = np.zeros([a,b,c])
    Txy = np.zeros([a,b,c +1])
    Tyz = np.zeros([a +1 ,b,c])
    Txz = np.zeros([a,b +1 ,c])
    
    for i in range(0,a +1): # initialize face temperatures
        for j in range(0,b):
            for k in range(0,c):
                if(np.mod(j,4) == 1 and np.mod(k,2) == 1): # Hot channel
                    Tyz[i,j,k] = I.T_hot_in #- i*1e2
               
                if(np.mod(j,4) == 1 and np.mod(k,2) == 0): # Alternate wall nodes between Hot channels
                    Tyz[i,j,k] = I.T_hot_in
                    
                    
                if(np.mod(j,4) == 3 and np.mod(k,2) == 1): # Cold channel
                    Tyz[a-i,j,k] = I.T_cold_in #- i*1e-2     

                if(np.mod(j,4) == 3 and np.mod(k,2) == 0): # Alternate wall nodes between cold channels
                    Tyz[i,j,k] = I.T_cold_in                     
                
                if(np.mod(j,2) == 0): # Walls between a layer of hot and cold channels
                    Tyz[i,j,k] = (I.T_cold_in + I.T_hot_in)/2
                    if(j ==0):
                        Tyz[i,j,k] ==I.T_hot_in # Initialize top wall(Y = Ymax) as Thot

                    if(j ==b -1):
                        Tyz[i,j,k] == I.T_cold_in# Initialize bottom wall(Y = 0) as Tcold
                        
    for i in range(0,I.a):# initialize node temperatures as mean of face temperatures
        for j in range(0,I.b):
            for k in range(0,I.c):
                Tg[i,j,k] = (Tyz[i,j,k] + Tyz[i + 1,j,k])/2
                

    return(Tg,Txy,Tyz,Txz)

def pressure():
    a = I.a
    b = I.b
    c = I.c
    Pg = np.zeros([a,b,c])
    P = np.zeros([a +1 ,b,c])

    
    for i in range(0,a +1): # initialize face pressure
        for j in range(0,b):
            for k in range(0,c):
                if(np.mod(j,4) == 1 and np.mod(k,2) == 1): # Hot channel
                    P[i,j,k] = I.P_hot_in
               
                if(np.mod(j,4) == 1 and np.mod(k,2) == 0): # Alternate wall nodes between Hot channels
                    P[i,j,k] = 0
                    
                    
                if(np.mod(j,4) == 3 and np.mod(k,2) == 1): # Cold channel
                    P[i,j,k] = I.P_cold_in     

                if(np.mod(j,4) == 3 and np.mod(k,2) == 0): # Alternate wall nodes between cold channels
                    P[i,j,k] = 0
                
                if(np.mod(j,2) == 0): # Walls between a layer of hot and cold channels
                    P[i,j,k] = 0                        
    for i in range(0,I.a):# initialize node Pressures as mean of face temperatures
        for j in range(0,I.b):
            for k in range(0,I.c):
                Pg[i,j,k] = (P[i,j,k] + P[i + 1,j,k])/2
                

    return(Pg,P,Pg,P)
  