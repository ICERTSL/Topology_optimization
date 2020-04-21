# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 02:26:42 2020

@author: HP
"""
import numpy as np
import matplotlib.pyplot as plt 
from CoolProp.CoolProp import PropsSI
import inputs as I
import initialize
import find
import BC

def pressure(P,Pg,Pgstat,Pstat,Tg):# Calculates pressure drop using pressure and temperature and updates pressure matrices
    a = I.a
    b = I.b
    c = I.c
    
    for i in range(0,a): # initialize face pressure
        for j in range(0,b):
            for k in range(0,c):
                if(np.mod(j,4) == 1 and np.mod(k,2) == 1): # Hot channel
                    deltap = find.pressure_drop(Tg[i,j,k],Pg[i,j,k],I.hot_fluid)
                    P[i+1,j,k] = Pstat[i,j,k] - deltap
                                               
                if(np.mod(j,4) == 3 and np.mod(k,2) == 1): # Cold channel

                    P[a-i-1,j,k] = Pstat[a-i,j,k] - find.pressure_drop(Tg[a-i-1,j,k],Pg[a-i-1,j,k],I.hot_fluid)
               
    for i in range(0,I.a):# initialize node Pressures as mean of face temperatures
        for j in range(0,I.b):
            for k in range(0,I.c):
                Pg[i,j,k] = (P[i,j,k] + P[i + 1,j,k])/2

             
    return(Pg, P)

def temperature(target,Tg,Tyz):
    T = np.zeros([I.a,I.b,I.c])
    for index in range(0,np.size(target)):
        i,j,k = find.find_ijk(index)
        i = int(i)
        j = int(j)
        k = int(k)
        if(np.mod(j,2) == 1 and np.mod(k,2) == 1): #liquid channel
            Tyz[i,j,k] = target[index]

        else:
            Tg[i,j,k] = target[index]

    for i in range(0,I.a):
        for j in range(0,I.b):
            for k in range(0,I.c):
                if(np.mod(j,2) == 1 and np.mod(k,2) == 1): #liquid channel
                    Tg[i,j,k] = 0.5*(Tyz[i,j,k] + Tyz[i+1,j,k])
    return(Tg, Tyz)

def Rc(x,y,i,j,k):# finds convective resistance x- temp, y - pressure
    nu = find.Nu(x,y,i,j,k)
    kliq = find.k_liquid(x,y,i,j,k)
    htc = (kliq*nu)/I.d
    area = I.d*I.delta_x
    resistance = 1/(htc*area)
    return(resistance)

   
    
    
        
def linear_equation_system(Rwx,Rwy,Rwz,Tg,Tyz,Pg,P):# Generates coefficient matrix A and constant matrix C 
    # WALL BOUNDARY CONDITIONS ARE NOT IMPOSED IN THIS FUNCTION
    a = I.a
    b = I.b
    c = I.c
    A = np.zeros([a*b*c,a*b*c])
    C = np.zeros([a*b*c,1])
    for index in range(0,np.size(C)):
        i,j,k = find.find_equation_id(index)
        # print("i j k",i,j,k)
        i = int(i)
        j = int(j)
        k = int(k)

        if(j != 0 and j!= b-1 and k != 0 and k != c-1): # except boundary walls
            
            if(np.mod(j,4) == 1 and np.mod(k,2) == 1): # Hot channel
                
                # print("Hot channel")
                # print("Node of balance",i,j,k)
                # print(index)
                C11 = 1/(Rwy + Rc(Tg[i,j,k],Pg[i,j,k],i,j,k))
                C12 = 1/(Rwy + Rc(Tg[i,j,k],Pg[i,j,k],i,j,k)) 
                C21 = 1/(Rwz + Rc(Tg[i,j,k],Pg[i,j,k],i,j,k))
                C22 = 1/(Rwz + Rc(Tg[i,j,k],Pg[i,j,k],i,j,k))
                C01 = -C11 - C21 + I.m_channel*find.Cp(Tg[i,j,k],Pg[i,j,k],i,j,k)   
                C02 = -C11 - C21 - I.m_channel*find.Cp(Tg[i,j,k],Pg[i,j,k],i,j,k)                                                        

                A[index,find.find_index(i,j-1,k)] = C11# Tg j-1 term
                A[index,find.find_index(i,j+1,k)] = C12# Tg j+1 term
                A[index,find.find_index(i,j,k-1)] = C21# Tg k-1 term
                A[index,find.find_index(i,j,k+1)] = C22# Tg k+1 term
                                                                               
                if(i == 0): # Ti term is hot in BC
                    C[index,0] = -I.T_hot_in*C01# BC = Thot in
                    A[index,find.find_index(i+1,j,k)] = C02# Tyz i+1 term
                    
                else: # interior
                    A[index,find.find_index(i,j,k)] = C01 #Tyz i term
                    A[index,find.find_index(i+1,j,k)] = C02# Tyz i+1 term        

#end of Hot channel
#############################################################################################
#Cold Channel
                    
            if(np.mod(j,4) == 3 and np.mod(k,2) == 1): # Cold channel
                
                C11 = 1/(Rwy + Rc(Tg[i,j,k],Pg[i,j,k],i,j,k))
                C12 = 1/(Rwy + Rc(Tg[i,j,k],Pg[i,j,k],i,j,k))
                C21 = 1/(Rwz + Rc(Tg[i,j,k],Pg[i,j,k],i,j,k))
                C22 = 1/(Rwz + Rc(Tg[i,j,k],Pg[i,j,k],i,j,k))
                C01 = -C11 - C21 - I.m_channel*find.Cp(Tg[i,j,k],Pg[i,j,k],i,j,k)   
                C02 = -C11 - C21 + I.m_channel*find.Cp(Tg[i,j,k],Pg[i,j,k],i,j,k)                                                        

                A[index,find.find_index(i,j-1,k)] = C11# Tg j-1 term
                A[index,find.find_index(i,j+1,k)] = C12# Tg j+1 term
                A[index,find.find_index(i,j,k-1)] = C21# Tg k-1 term
                A[index,find.find_index(i,j,k+1)] = C22# Tg k+1 term
                
                if(i == a-1): # Ti-1 term is adjusted - Cold outlet
                    A[index,find.find_index(i,j,k)] = C01 #Tyz i term
                    C[index,0] = -I.T_cold_in*C02#Cold in BC-> Tyz i+1 term
                    
                else: # interior
                    A[index,find.find_index(i,j,k)] = C01 #Tyz i term
                    A[index,find.find_index(i+1,j,k)] = C02# Tyz i+1 term  

#end of cold channel
#############################################################################################
#Wall nodes                    
                
        
            if(np.mod(j,2) == 0 and np.mod(k,2) == 1): # Type 1 Solid node
                C01 = 1/(2*Rwx)
                C02 = 1/(2*Rwx)
                C11 = 0.5/(Rwy + Rc(Tg[i,j-1,k],Pg[i,j-1,k],i,j-1,k))
                C12 = 0.5/(Rwy + Rc(Tg[i,j+1,k],Pg[i,j+1,k],i,j+1,k))
                C21 = 1/(2*Rwz)
                C22 = 1/(2*Rwz)
 
                if(i == 0):# Hot in BC

                    if(np.mod(j+1,4) == 1):# hot channet at j+1 and cold at j-1 Hot in BC at i,j+1
                        C0 = -C02 - 2*C11 - 2*C12 - C21 - C22
                        A[index,find.find_index(i,j,k)] = C0 #i,j,k
                        A[index,find.find_index(i+1,j,k)] = C02 #i+1,j,k
                        A[index,find.find_index(i,j,k-1)] = C21 #i,j,k-1
                        A[index,find.find_index(i,j,k+1)] = C22 #i,j,k+1
                        A[index,find.find_index(i,j-1,k)] = C11 #i,j-1,k
                        A[index,find.find_index(i+1,j-1,k)] = C11
                        A[index,find.find_index(i+1,j+1,k)] = C12
                        C[index,0] = -I.T_hot_in*C12 
                                              
                        
                    elif(np.mod(j-1,4) ==1):# Hot channel at j-1 Hot in BC at i,j-1
                        C0 = -C02 - 2*C11 - 2*C12 - C21 - C22
                        A[index,find.find_index(i,j,k)] = C0
                        A[index,find.find_index(i+1,j,k)] = C02
                        A[index,find.find_index(i,j,k-1)] = C21
                        A[index,find.find_index(i,j,k+1)] = C22
                        A[index,find.find_index(i+1,j-1,k)] = C11
                        A[index,find.find_index(i,j+1,k)] = C12
                        A[index,find.find_index(i+1,j+1,k)] = C12
                        C[index,0] = -I.T_hot_in*C11

                        
                elif(i == a-1):# Cold in BC
                    
                    if(np.mod(j+1,4) == 1):# hot channet at j+1 and cold at j-1 Cold in BC at i+1,j-1
                        C0 = -C02 - 2*C11 - 2*C12 - C21 - C22
                        A[index,find.find_index(i,j,k)] = C0
                        A[index,find.find_index(i-1,j,k)] = C01
                        A[index,find.find_index(i,j,k-1)] = C21
                        A[index,find.find_index(i,j,k+1)] = C22
                        A[index,find.find_index(i,j+1,k)] = C12
                        A[index,find.find_index(i+1,j+1,k)] = C12
                        A[index,find.find_index(i,j-1,k)] = C11
                        C[index,0] = -I.T_cold_in*C11 
                                          
                    elif(np.mod(j-1,4) ==1):# Hot channel at j-1 Cold in BC at i+1,j+1
                        C0 = -C02 - 2*C11 - 2*C12 - C21 - C22
                        A[index,find.find_index(i,j,k)] = C0
                        A[index,find.find_index(i-1,j,k)] = C01
                        A[index,find.find_index(i,j,k-1)] = C21
                        A[index,find.find_index(i,j,k+1)] = C22
                        A[index,find.find_index(i,j+1,k)] = C12
                        A[index,find.find_index(i,j-1,k)] = C11
                        A[index,find.find_index(i+1,j-1,k)] = C11
                        C[index,0] = -I.T_cold_in*C12                 
                else:
                    C0 = -C02 - C01 - 2*C11 - 2*C12 - C21 - C22
                    A[index,find.find_index(i,j,k)] = C0
                    A[index,find.find_index(i-1,j,k)] = C01
                    A[index,find.find_index(i+1,j,k)] = C02
                    A[index,find.find_index(i,j,k-1)] = C21
                    A[index,find.find_index(i,j,k+1)] = C22
                    A[index,find.find_index(i,j+1,k)] = C12
                    A[index,find.find_index(i+1,j+1,k)] = C12
                    A[index,find.find_index(i,j-1,k)] = C11
                    A[index,find.find_index(i+1,j-1,k)] = C11

                    
            if(np.mod(j,2) == 0 and np.mod(k,2) == 0 ): # Type 2 Solid node
                C01 = 1/(2*Rwx)
                C02 = 1/(2*Rwx)
                C11 = 1/(2*Rwy)
                C12 = 1/(2*Rwy)
                C21 = 1/(2*Rwz)
                C22 = 1/(2*Rwz)
                if(i == 0):
                        C0 = -C02 - C11 - C12 - C21 - C22
                        A[index,find.find_index(i,j,k)] = C0
                        A[index,find.find_index(i+1,j,k)] = C02
                        A[index,find.find_index(i,j,k-1)] = C21
                        A[index,find.find_index(i,j,k+1)] = C22
                        A[index,find.find_index(i,j+1,k)] = C12
                        A[index,find.find_index(i,j-1,k)] = C11

                    
                elif(i == a-1):
                        C0 = -C01 - C11 - C12 - C21 - C22
                        A[index,find.find_index(i,j,k)] = C0
                        A[index,find.find_index(i-1,j,k)] = C01
                        A[index,find.find_index(i,j,k-1)] = C21
                        A[index,find.find_index(i,j,k+1)] = C22
                        A[index,find.find_index(i,j+1,k)] = C12
                        A[index,find.find_index(i,j-1,k)] = C11
                    
                else:
                        C0 = -C01 - C02 - C11 - C12 - C21 - C22
                        A[index,find.find_index(i,j,k)] = C0
                        A[index,find.find_index(i-1,j,k)] = C01
                        A[index,find.find_index(i+1,j,k)] = C02
                        A[index,find.find_index(i,j,k-1)] = C21
                        A[index,find.find_index(i,j,k+1)] = C22
                        A[index,find.find_index(i,j+1,k)] = C12
                        A[index,find.find_index(i,j-1,k)] = C11                
                
                
                
                
            if(np.mod(j,2) == 1 and np.mod(k,2) == 0 ): # Type 3 Solid node

                C01 = 1/(2*Rwx)
                C02 = 1/(2*Rwx)
                C11 = 1/(2*Rwy)
                C12 = 1/(2*Rwy)
                C21 = 0.5/(Rwz + Rc(Tg[i,j,k-1],Pg[i,j,k-1],i,j,k-1))
                C22 = 0.5/(Rwz + Rc(Tg[i,j,k+1],Pg[i,j,k+1],i,j,k+1))
        
               
                if(i == 0):# Hot in BC

                    if(np.mod(j,4) == 1):# hot channel at j and cold at j-1, j+1-> Hot in BC at i,k-1 k+1
                        C0 = -C02 - C11 - C12 - 2*C21 - 2*C22
                        A[index,find.find_index(i,j,k)] = C0 #i,j,k
                        A[index,find.find_index(i+1,j,k)] = C02 #i+1,j,k
                        A[index,find.find_index(i,j-1,k)] = C11 #i,j-1,k
                        A[index,find.find_index(i,j+1,k)] = C12 #i,j+1,k
                        A[index,find.find_index(i+1,j,k-1)] = C21 #i+1,j,k-1
                        A[index,find.find_index(i+1,j,k+1)] = C22 #i+1,j,k+1
                        
                        C[index,0] = -I.T_hot_in*(C21 + C22) 
                                              
                        
                    elif(np.mod(j,4) ==3):# Cold channet at j and hot at j-1, j+1
                        C0 = -C02 - C11 - C12 - 2*C21 - 2*C22
                        A[index,find.find_index(i,j,k)] = C0 #i,j,k
                        A[index,find.find_index(i+1,j,k)] = C02 #i+1,j,k
                        A[index,find.find_index(i,j-1,k)] = C11 #i,j-1,k
                        A[index,find.find_index(i,j+1,k)] = C12 #i,j+1,k
                        A[index,find.find_index(i+1,j,k-1)] = C21 #Tyz i+1,j,k-1
                        A[index,find.find_index(i+1,j,k+1)] = C22 #Tyz i+1,j,k+1
                        A[index,find.find_index(i,j,k-1)] = C21 #Tyz i,j,k-1
                        A[index,find.find_index(i,j,k+1)] = C22 #Tyz i+,j,k+1
                        
                        
                elif(i == a-1):# Cold in BC
                    
                    if(np.mod(j,4) == 1):# hot channet at j
                        C0 = -C01 - C11 - C12 - 2*C21 - 2*C22
                        A[index,find.find_index(i,j,k)] = C0#ijk
                        A[index,find.find_index(i-1,j,k)] = C01 #i-1,j,k
                        A[index,find.find_index(i,j-1,k)] = C11#i,j-1,k
                        A[index,find.find_index(i,j+1,k)] = C12#i,j+1,k
                        A[index,find.find_index(i,j,k-1)] = C21#i,j,k-1
                        A[index,find.find_index(i,j,k+1)] = C22#i,j,k+1
                        A[index,find.find_index(i+1,j,k-1)] = C21#i+1,j,k-1
                        A[index,find.find_index(i+1,j,k+1)] = C22#i+1,j,k+1

                        

                                          
                    elif(np.mod(j,4) ==3):# Cold channel at j Cold in BC at i+1,k-1,+1
                        C0 = -C01 - C11 - C12 - 2*C21 - 2*C22
                        A[index,find.find_index(i,j,k)] = C0
                        A[index,find.find_index(i-1,j,k)] = C01
                        A[index,find.find_index(i,j+1,k)] = C12
                        A[index,find.find_index(i,j-1,k)] = C11
                        A[index,find.find_index(i,j,k-1)] = C21
                        A[index,find.find_index(i,j,k+1)] = C22
                        C[index,0] = -I.T_cold_in*(C21 + C22)                 
                else:
                    C0 = -C02 - C01 - C11 - C12 - 2*C21 - 2*C22
                    A[index,find.find_index(i,j,k)] = C0
                    A[index,find.find_index(i-1,j,k)] = C01
                    A[index,find.find_index(i+1,j,k)] = C02
                    A[index,find.find_index(i,j,k-1)] = C21
                    A[index,find.find_index(i,j,k+1)] = C22
                    A[index,find.find_index(i+1,j,k+1)] = C22
                    A[index,find.find_index(i+1,j,k-1)] = C21
                    A[index,find.find_index(i,j+1,k)] = C12
                    A[index,find.find_index(i,j-1,k)] = C11
               
  
    A,C = BC.adiabatic_wall_BC(A,C,Rwx,Rwy,Rwz,Tg,Tyz,Pg,P)                
    return(A,C)


