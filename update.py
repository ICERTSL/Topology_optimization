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

def pressure(P,Pg,Tg):# Calculates pressure drop using pressure and temperature and updates pressure matrices
    a = I.a
    b = I.b
    c = I.c
    
    for i in range(0,a): # initialize face pressure
        for j in range(0,b):
            for k in range(0,c):
                if(np.mod(j,4) == 1 and np.mod(k,2) == 1): # Hot channel
                    deltap = find.pressure_drop(Tg[i,j,k],Pg[i,j,k],I.hot_fluid)
                    P[i+1,j,k] = P[i,j,k] - deltap
                                               
                if(np.mod(j,4) == 3 and np.mod(k,2) == 1): # Cold channel
                    P[a-i-1,j,k] = P[a-i,j,k] - find.pressure_drop(Tg[a-i-1,j,k],Pg[a-i-1,j,k],I.hot_fluid)
               
    for i in range(0,I.a):# initialize node Pressures as mean of face temperatures
        for j in range(0,I.b):
            for k in range(0,I.c):
                Pg[i,j,k] = (P[i,j,k] + P[i + 1,j,k])/2
                
    return(Pg, P)

def resistance(P,Tyz,Pg,Tg): #generates resistance matrix
    a = I.a
    b = I.b
    c = I.c
    R = np.zeros([a,b,c,3,4])


    for i in range(0,a):

        for j in range(0,b):
            for k in range(0,c):

                if(np.mod(j,4) == 1 and np.mod(k,2) == 1): # Hot channel

                    htc = find.k_liquid(Tg[i,j,k],Pg[i,j,k],i,j,k)
                    htc = htc*find.Nu(Tg[i,j,k],Pg[i,j,k],i,j,k)
                    htc = htc/I.d # Convective Heat transfer coefficient - phi
                    R[i,j,k,2,1] = a/(htc*I.d*I.L)#Rfz2
                    R[i,j,k,2,3] = a/(htc*I.d*I.L)#Rfz4                                               
                    R[i,j,k,1,1] = a/(htc*I.d*I.L)#Rfy2
                    R[i,j,k,1,3] = a/(htc*I.d*I.L)#Rfy4
                    if(i == 0): # Hot inlet
                        print("if 1 i",i)
                        dh1 = PropsSI('H','T',Tg[i+1,j,k],'P',Pg[i+1,j,k],I.hot_fluid) - PropsSI('H','T',Tg[i,j,k],'P',Pg[i,j,k],I.hot_fluid) 
                        dT1 = Tg[i+1,j,k] - Tg[i,j,k]                       
                        dh2 = PropsSI('H','T',Tyz[i,j,k],'P',P[i,j,k],I.hot_fluid) - PropsSI('H','T',Tg[i,j,k],'P',Pg[i,j,k],I.hot_fluid) 
                        dT2 = Tyz[i,j,k] - Tg[i,j,k]
                        
                    elif(i == a-1): # Hot out
                        print("if 2 i",i)
                        dh1 = PropsSI('H','T',Tyz[i+1,j,k],'P',P[i+1,j,k],I.hot_fluid) - PropsSI('H','T',Tg[i,j,k],'P',Pg[i,j,k],I.hot_fluid) 
                        dT1 = Tyz[i+1,j,k] - Tg[i,j,k]                       
                        dh2 = PropsSI('H','T',Tg[i-1,j,k],'P',Pg[i-1,j,k],I.hot_fluid) - PropsSI('H','T',Tg[i,j,k],'P',Pg[i,j,k],I.hot_fluid) 
                        dT2 = Tg[i-1,j,k] - Tg[i,j,k]
                     
                    else: # Hot interior
                        print("else i",i)
                        dh1 = PropsSI('H','T',Tg[i+1,j,k],'P',Pg[i+1,j,k],I.hot_fluid) - PropsSI('H','T',Tg[i,j,k],'P',Pg[i,j,k],I.hot_fluid) 
                        dT1 = Tyz[i+1,j,k] - Tg[i,j,k]                       
                        dh2 = PropsSI('H','T',Tg[i-1,j,k],'P',Pg[i-1,j,k],I.hot_fluid) - PropsSI('H','T',Tg[i,j,k],'P',Pg[i,j,k],I.hot_fluid) 
                        dT2 = Tg[i-1,j,k] - Tg[i,j,k]
                        
                    R[i,j,k,0,0] = dT1/(I.m_channel*dh1)#Rfx1
                    R[i,j,k,0,1] = dT2/(I.m_channel*dh2)#Rfx2
                    
                        
                        
                    
                    
                if(np.mod(j,4) == 3 and np.mod(k,2) == 1): # Cold channel

                    htc = find.k_liquid(Tg[i,j,k],Pg[i,j,k],i,j,k)*find.Nu(Tg[i,j,k],Pg[i,j,k],i,j,k)
                    htc = htc/I.d # Convective Heat transfer coefficient - phi
                    R[i,j,k,2,1] = a/(htc*I.d*I.L)#Rfz2
                    R[i,j,k,2,3] = a/(htc*I.d*I.L)#Rfz4                                               
                    R[i,j,k,1,1] = a/(htc*I.d*I.L)#Rfy2
                    R[i,j,k,1,3] = a/(htc*I.d*I.L)#Rfy4
                    if(i == 0): # Cold out
                        dh1 = PropsSI('H','T',Tg[i+1,j,k],'P',Pg[i+1,j,k],I.cold_fluid) - PropsSI('H','T',Tg[i,j,k],'P',Pg[i,j,k],I.cold_fluid) 
                        dT1 = Tg[i+1,j,k] - Tg[i,j,k]                       
                        dh2 = PropsSI('H','T',Tyz[i,j,k],'P',P[i,j,k],I.cold_fluid) - PropsSI('H','T',Tg[i,j,k],'P',Pg[i,j,k],I.cold_fluid) 
                        dT2 = Tyz[i,j,k] - Tg[i,j,k]
                        
                    elif(i == a-1): # Cold in
                        dh1 = PropsSI('H','T',Tyz[i+1,j,k],'P',P[i+1,j,k],I.cold_fluid) - PropsSI('H','T',Tg[i,j,k],'P',Pg[i,j,k],I.cold_fluid) 
                        dT1 = Tyz[i+1,j,k] - Tg[i,j,k]                       
                        dh2 = PropsSI('H','T',Tg[i-1,j,k],'P',Pg[i-1,j,k],I.cold_fluid) - PropsSI('H','T',Tg[i,j,k],'P',Pg[i,j,k],I.cold_fluid) 
                        dT2 = Tg[i-1,j,k] - Tg[i,j,k]
                     
                    else: # Cold interior
                        dh1 = PropsSI('H','T',Tg[i+1,j,k],'P',Pg[i+1,j,k],I.cold_fluid) - PropsSI('H','T',Tg[i,j,k],'P',Pg[i,j,k],I.cold_fluid) 
                        dT1 = Tyz[i+1,j,k] - Tg[i,j,k]                       
                        dh2 = PropsSI('H','T',Tg[i-1,j,k],'P',Pg[i-1,j,k],I.cold_fluid) - PropsSI('H','T',Tg[i,j,k],'P',Pg[i,j,k],I.cold_fluid) 
                        dT2 = Tg[i-1,j,k] - Tg[i,j,k]
                        
                    R[i,j,k,0,0] = -dT1/(I.m_channel*dh1)#Rfx1
                    R[i,j,k,0,1] = -dT2/(I.m_channel*dh2)#Rfx2
                    
            
                if(np.mod(j,2) == 0 and np.mod(k,2) == 1): # Type 1 Solid node
                    deltax = I.L/a
                    deltay = I.t
                    deltaz = I.d
                    R[i,j,k,0,0] = deltax/(2*I.k_solid*deltay*deltaz)# Rwx1
                    R[i,j,k,0,1] = deltax/(2*I.k_solid*deltay*deltaz)# Rwx2
                    R[i,j,k,1,0] = deltay/(2*I.k_solid*deltax*deltaz)# Rwy1
                    R[i,j,k,1,1] = deltay/(2*I.k_solid*deltax*deltaz)# Rwy2
                    R[i,j,k,2,0] = deltaz/(2*I.k_solid*deltay*deltax)# Rwz1
                    R[i,j,k,2,1] = deltaz/(2*I.k_solid*deltay*deltax)# Rwz2
                        
                if(np.mod(j,2) == 0 and np.mod(k,2) == 0): # Type 2 Solid node
                    deltax = I.L/a
                    deltay = I.t
                    deltaz = I.w
                    R[i,j,k,0,0] = deltax/(2*I.k_solid*deltay*deltaz)# Rwx1
                    R[i,j,k,0,1] = deltax/(2*I.k_solid*deltay*deltaz)# Rwx2
                    R[i,j,k,1,0] = deltay/(2*I.k_solid*deltax*deltaz)# Rwy1
                    R[i,j,k,1,1] = deltay/(2*I.k_solid*deltax*deltaz)# Rwy2
                    R[i,j,k,2,0] = deltaz/(2*I.k_solid*deltay*deltax)# Rwz1
                    R[i,j,k,2,1] = deltaz/(2*I.k_solid*deltay*deltax)# Rwz2
                    
                if(np.mod(j,2) == 1 and np.mod(k,2) == 0): # Type 3 Solid node
                    deltax = I.L/a
                    deltay = I.d
                    deltaz = I.w
                    R[i,j,k,0,0] = deltax/(2*I.k_solid*deltay*deltaz)# Rwx1
                    R[i,j,k,0,1] = deltax/(2*I.k_solid*deltay*deltaz)# Rwx2
                    R[i,j,k,1,0] = deltay/(2*I.k_solid*deltax*deltaz)# Rwy1
                    R[i,j,k,1,1] = deltay/(2*I.k_solid*deltax*deltaz)# Rwy2
                    R[i,j,k,2,0] = deltaz/(2*I.k_solid*deltay*deltax)# Rwz1
                    R[i,j,k,2,1] = deltaz/(2*I.k_solid*deltay*deltax)# Rwz2  
        
    return(R)
    
def linear_equation_system(R,Tg,Tyz,Pg,P):# Generates coefficient matrix A and constant matrix C 
    # WALL BOUNDARY CONDITIONS ARE NOT IMPOSED IN THIS FUNCTION
    a = I.a
    b = I.b
    c = I.c
    A = np.zeros([a*b*c,a*b*c])
    C = np.zeros([a*b*c,1])
    for index in range(0,np.size(C)):
        i,j,k = find.find_ijk(index)
        # print("i j k",i,j,k)
        i = int(i)
        j = int(j)
        k = int(k)
        if(j != 0 and j!= b-1 and k != 0 and k != c-1): # except boundary walls
            
            if(np.mod(j,2) == 1 and np.mod(k,2) == 1): # fluid channel

                if(i == 0): # Ti-1 term is adjusted
                    C[index,0] = -Tyz[i,j,k]/(R[i,j,k,0,1])
                    A[index,find.find_index(i+1,j,k)] = 1/R[i,j,k,0,0]
                    
                elif(i == a-1): # Ti+1 term is adjusted
                    C[index,0] = -Tyz[i+1,j,k]/(R[i,j,k,0,0])
                    A[index,find.find_index(i-1,j,k)] = 1/R[i,j,k,0,1]                     
                else: # interior
                    A[index,find.find_index(i-1,j,k)] = 1/R[i,j,k,0,1]# i-1 term
                    A[index,find.find_index(i+1,j,k)] = 1/R[i,j,k,0,0]# i+1 term
                    

                    

                A[index,find.find_index(i,j+1,k)] = 1/(R[i,j,k,1,0] + R[i,j,k,1,1] + R[i,j+1,k,1,1])#j+1 term
                A[index,find.find_index(i,j-1,k)] = 1/(R[i,j,k,1,2] + R[i,j,k,1,3] + R[i,j-1,k,1,0])#j-1 term
                A[index,find.find_index(i,j,k+1)] = 1/(R[i,j,k,2,0] + R[i,j,k,2,1] + R[i,j,k+1,2,1])# k+1 term
                A[index,find.find_index(i,j,k-1)] = 1/(R[i,j,k,2,2] + R[i,j,k,2,3] + R[i,j,k-1,2,0])# k-1 term
                
                coeff = 1/(R[i,j,k,1,0] + R[i,j,k,1,1] + R[i,j+1,k,1,1]) + 1/(R[i,j,k,1,2] + R[i,j,k,1,3] +  R[i,j-1,k,1,0]) + 1/(R[i,j,k,2,0] + R[i,j,k,2,1] + R[i,j,k+1,2,1]) + 1/(R[i,j,k,2,2] + R[i,j,k,2,3] + R[i,j,k-1,2,0])
                coeff = coeff + 1/R[i,j,k,0,1] + 1/R[i,j,k,0,0]
                
                A[index,find.find_index(i,j,k)] = -coeff #i,j,k term
                    
                
        
            if(np.mod(j,2) == 0 and np.mod(k,2) == 1): # Type 1 Solid node

                if(i == 0):
                    A[index,find.find_index(i+1,j,k)] = 1/(R[i,j,k,0,0] + R[i+1,j,k,0,1])# i+1 term
                    temp = A[index,find.find_index(i+1,j,k)]
                    
                elif(i == a-1):
                    A[index,find.find_index(i-1,j,k)] = 1/(R[i,j,k,0,1] + R[i-1,j,k,0,0])# i-1 term
                    temp = A[index,find.find_index(i-1,j,k)]
                    
                else:
                    A[index,find.find_index(i-1,j,k)] = 1/(R[i,j,k,0,1] + R[i-1,j,k,0,0])# i-1 term
                    A[index,find.find_index(i+1,j,k)] = 1/(R[i,j,k,0,0] + R[i+1,j,k,0,1])# i+1 term
                    temp = A[index,find.find_index(i+1,j,k)] + A[index,find.find_index(i-1,j,k)]
                    

                A[index,find.find_index(i,j-1,k)] = 1/(R[i,j,k,1,1] + R[i,j-1,k,1,0] + R[i,j-1,k,1,1])# j-1 term
                A[index,find.find_index(i,j+1,k)] = 1/(R[i,j,k,1,0] + R[i,j+1,k,1,2] + R[i,j+1,k,1,3])# j+1 term
                A[index,find.find_index(i,j,k-1)] = 1/(R[i,j,k,2,1] + R[i,j,k-1,2,0])# k-1 term
                A[index,find.find_index(i,j,k+1)] = 1/(R[i,j,k,2,0] + R[i,j,k+1,2,1])# k+1 term

                

                coeff = temp + A[index,find.find_index(i,j-1,k)] + A[index,find.find_index(i,j+1,k)] + A[index,find.find_index(i,j,k-1)] +A[index,find.find_index(i,j,k+1)]
                A[index,find.find_index(i,j,k)] = -coeff
                    
            if(np.mod(j,2) == 0 and np.mod(k,2) == 0 ): # Type 2 Solid node

                if(i == 0):
                    A[index,find.find_index(i+1,j,k)] = 1/(R[i,j,k,0,0] + R[i+1,j,k,0,1])# i+1 term
                    temp = A[index,find.find_index(i+1,j,k)]
                    
                elif(i == a-1):
                    A[index,find.find_index(i-1,j,k)] = 1/(R[i,j,k,0,1] + R[i-1,j,k,0,0])# i-1 term
                    temp = A[index,find.find_index(i-1,j,k)]
                    
                else:
                    A[index,find.find_index(i-1,j,k)] = 1/(R[i,j,k,0,1] + R[i-1,j,k,0,0])# i-1 term
                    A[index,find.find_index(i+1,j,k)] = 1/(R[i,j,k,0,0] + R[i+1,j,k,0,1])# i+1 term
                    temp = A[index,find.find_index(i+1,j,k)] + A[index,find.find_index(i-1,j,k)]                
                
                
                
                
                A[index,find.find_index(i,j-1,k)] = 1/(R[i,j,k,1,1] + R[i,j-1,k,1,0])# j-1 term
                A[index,find.find_index(i,j+1,k)] = 1/(R[i,j,k,1,0] + R[i,j+1,k,1,1])# j+1 term
                A[index,find.find_index(i,j,k-1)] = 1/(R[i,j,k,2,1] + R[i,j,k-1,2,0])# k-1 term
                A[index,find.find_index(i,j,k+1)] = 1/(R[i,j,k,2,0] + R[i,j,k+1,2,1])# k+1 term                   

                coeff = temp + A[index,find.find_index(i,j-1,k)] + A[index,find.find_index(i,j+1,k)] + A[index,find.find_index(i,j,k-1)] +A[index,find.find_index(i,j,k+1)]
                A[index,find.find_index(i,j,k)] = -coeff

                
            if(np.mod(j,2) == 1 and np.mod(k,2) == 0 ): # Type 3 Solid node
                
                if(i == 0):
                    A[index,find.find_index(i+1,j,k)] = 1/(R[i,j,k,0,0] + R[i+1,j,k,0,1])# i+1 term
                    temp = A[index,find.find_index(i+1,j,k)]
                    
                elif(i == a-1):
                    A[index,find.find_index(i-1,j,k)] = 1/(R[i,j,k,0,1] + R[i-1,j,k,0,0])# i-1 term
                    temp = A[index,find.find_index(i-1,j,k)]
                    
                else:
                    A[index,find.find_index(i-1,j,k)] = 1/(R[i,j,k,0,1] + R[i-1,j,k,0,0])# i-1 term
                    A[index,find.find_index(i+1,j,k)] = 1/(R[i,j,k,0,0] + R[i+1,j,k,0,1])# i+1 term
                    temp = A[index,find.find_index(i+1,j,k)] + A[index,find.find_index(i-1,j,k)]  
                    

                A[index,find.find_index(i,j-1,k)] = 1/(R[i,j,k,1,1] + R[i,j-1,k,1,0])# j-1 term
                A[index,find.find_index(i,j+1,k)] = 1/(R[i,j,k,1,0] + R[i,j+1,k,1,1])# j+1 term
                A[index,find.find_index(i,j,k-1)] = 1/(R[i,j,k,2,1] + R[i,j,k-1,2,0] + R[i,j,k-1,2,1])# k-1 term
                A[index,find.find_index(i,j,k+1)] = 1/(R[i,j,k,2,0] + R[i,j,k+1,2,2] + R[i,j,k+1,2,3])# k+1 term                   
        
                coeff = temp + A[index,find.find_index(i,j-1,k)] + A[index,find.find_index(i,j+1,k)] + A[index,find.find_index(i,j,k-1)] +A[index,find.find_index(i,j,k+1)]
                A[index,find.find_index(i,j,k)] = -coeff
                
        
    A,C = adiabatic_wall_BC(A, C, R)                
    return(A,C)
def adiabatic_wall_BC(A,C,R):
    a = I.a
    b = I.b
    c = I.c
    
    for i in range(0,a):
###############################################################################################################
#side Walls        
        for j in range(0,b):# Side Walls
            index0 = find.find_index(i,j,0)
            index1 = find.find_index(i,j,c-1)       

            if(np.mod(j,2) == 0): # Type 2 Solid node

            
                if(i == 0): # Nodes along y axis at x = 0
                    A[index0,find.find_index(i+1,j,0)] = 1/(R[i,j,0,0,0] + R[i+1,j,0,0,0])# i+1 term
                    tempi0 = A[index0,find.find_index(i+1,j,0)]
                    
                    A[index1,find.find_index(i+1,j,c-1)] = 1/(R[i,j,c-1,0,0] + R[i+1,j,c-1,0,0])
                    tempi1 = A[index0,find.find_index(i+1,j,c-1)]
                    
                elif(i == a-1):#Nodes along y axis at x = L
                    A[index0,find.find_index(i-1,j,0)] = 1/(R[i,j,0,0,1] + R[i-1,j,0,0,0])# i-1 term
                    tempi0 = A[index0,find.find_index(i-1,j,0)]

                    A[index1,find.find_index(i-1,j,c-1)] = 1/(R[i,j,c-1,0,1] + R[i-1,j,c-1,0,0])
                    tempi1 = A[index0,find.find_index(i-1,j,c-1)]
                    
                else: # Nodes in interior of XY plane 
                    A[index0,find.find_index(i-1,j,0)] = 1/(R[i,j,0,0,1] + R[i-1,j,0,0,0])# i-1 term
                    A[index0,find.find_index(i+1,j,0)] = 1/(R[i,j,0,0,0] + R[i+1,j,0,0,1])# i+1 term
                    tempi0 = A[index0,find.find_index(i+1,j,0)] + A[index0,find.find_index(i-1,j,0)]

                    A[index1,find.find_index(i-1,j,c-1)] = 1/(R[i,j,c-1,0,1] + R[i-1,j,c-1,0,0])# i-1 term
                    A[index1,find.find_index(i+1,j,c-1)] = 1/(R[i,j,c-1,0,0] + R[i+1,j,c-1,0,1])# i+1 term
                    tempi1 = A[index1,find.find_index(i+1,j,c-1)] + A[index1,find.find_index(i-1,j,c-1)]
                
                if(j == 0):#Nodes along x axis at y = 0
                    A[index0,find.find_index(i,j+1,0)] = 1/(R[i,j,0,1,0] + R[i,j+1,0,1,1])# j+1 term
                    tempj0 = A[index0,find.find_index(i,j+1,0)]

                    A[index1,find.find_index(i,j+1,c-1)] = 1/(R[i,j,c-1,1,0] + R[i,j+1,c-1,1,1])# j+1 term
                    tempj1 = A[index1,find.find_index(i,j+1,c-1)]
                    
                elif(j == b-1):#Nodes along x axis at y = ymax
                    A[index0,find.find_index(i,j-1,0)] = 1/(R[i,j,0,1,1] + R[i,j-1,0,1,0])# j-1 term
                    tempj0 =  A[index0,find.find_index(i,j-1,0)]
                    
                    A[index1,find.find_index(i,j-1,c-1)] = 1/(R[i,j,c-1,1,1] + R[i,j-1,c-1,1,0])
                    tempj1 =  A[index1,find.find_index(i,j-1,c-1)]
                
                else: #interior of XY plane
                    A[index0,find.find_index(i,j-1,0)] = 1/(R[i,j,0,1,1] + R[i,j-1,0,1,0])# j-1 term                   
                    A[index0,find.find_index(i,j+1,0)] = 1/(R[i,j,0,1,0] + R[i,j+1,0,1,1])# j+1 term
                    tempj0 =  A[index0,find.find_index(i,j-1,0)] +  A[index0,find.find_index(i,j+1,0)]
                    
                    A[index1,find.find_index(i,j-1,c-1)] = 1/(R[i,j,c-1,1,1] + R[i,j-1,c-1,1,0])# j-1 term                   
                    A[index1,find.find_index(i,j+1,c-1)] = 1/(R[i,j,c-1,1,0] + R[i,j+1,c-1,1,1])# j+1 term
                    tempj1 =  A[index0,find.find_index(i,j-1,c-1)] +  A[index1,find.find_index(i,j+1,c-1)]                
                  
             
            
                A[index0,find.find_index(i,j,1)] = 1/(R[i,j,0,2,0] + R[i,j,1,2,1])# k+1 term                   
                A[index1,find.find_index(i,j,c-2)] = 1/(R[i,j,c-1,2,1] + R[i,j,c-2,2,0])# k-1 term
                
                coeff0 = tempi0 + tempj0 + A[index0,find.find_index(i,j,1)]
                coeff1 = tempi1 + tempj1 + A[index1,find.find_index(i,j,c-2)]
                A[index0,find.find_index(i,j,0)] = -coeff0
                A[index1,find.find_index(i,j,c-1)] = -coeff1
            
####################################################################################################                
            if(np.mod(j,2) == 1): # Type 3 Solid node
                
                if(i == 0):
                    A[index0,find.find_index(i+1,j,0)] = 1/(R[i,j,0,0,0] + R[i+1,j,0,0,1])# i+1 term
                    tempi0 = A[index0,find.find_index(i+1,j,0)]

                    A[index1,find.find_index(i+1,j,c-1)] = 1/(R[i,j,c-1,0,0] + R[i+1,j,c-1,0,1])# i+1 term
                    tempi1 = A[index1,find.find_index(i+1,j,c-1)]
                    
                elif(i == a-1):
                    A[index0,find.find_index(i-1,j,0)] = 1/(R[i,j,0,0,1] + R[i-1,j,0,0,0])# i-1 term
                    tempi0 = A[index0,find.find_index(i-1,j,0)]

                    A[index1,find.find_index(i-1,j,c-1)] = 1/(R[i,j,c-1,0,1] + R[i-1,j,c-1,0,0])# i-1 term
                    tempi1 = A[index1,find.find_index(i-1,j,c-1)]
                    
                else:
                    A[index0,find.find_index(i-1,j,0)] = 1/(R[i,j,0,0,1] + R[i-1,j,0,0,0])# i-1 term
                    A[index0,find.find_index(i+1,j,0)] = 1/(R[i,j,0,0,0] + R[i+1,j,0,0,1])# i+1 term
                    tempi0 = A[index0,find.find_index(i+1,j,0)] + A[index0,find.find_index(i-1,j,0)]  
                    
                    A[index1,find.find_index(i-1,j,c-1)] = 1/(R[i,j,c-1,0,1] + R[i-1,j,c-1,0,0])# i-1 term
                    A[index1,find.find_index(i+1,j,c-1)] = 1/(R[i,j,c-1,0,0] + R[i+1,j,c-1,0,1])# i+1 term
                    tempi1 = A[index1,find.find_index(i+1,j,c-1)] + A[index1,find.find_index(i-1,j,c-1)]               

                if(j == 0):
                    A[index0,find.find_index(i,j+1,0)] = 1/(R[i,j,0,1,0] + R[i,j+1,0,1,1])# j+1 term
                    tempj0 = A[index0,find.find_index(i,j+1,0)]

                    A[index1,find.find_index(i,j+1,c-1)] = 1/(R[i,j,c-1,1,0] + R[i,j+1,c-1,1,1])# j+1 term
                    tempj1 = A[index1,find.find_index(i,j+1,c-1)]
                    
                elif(j == b-1):
                    A[index0,find.find_index(i,j-1,0)] = 1/(R[i,j,0,1,1] + R[i,j-1,0,1,0])# j-1 term
                    tempj0 = A[index0,find.find_index(i,j-1,0)]

                    A[index1,find.find_index(i,j-1,c-1)] = 1/(R[i,j,c-1,1,1] + R[i,j-1,c-1,1,0])# j-1 term
                    tempj1 = A[index1,find.find_index(i,j-1,c-1)]
                    
                else:
                    A[index0,find.find_index(i,j-1,0)] = 1/(R[i,j,0,1,1] + R[i,j-1,0,1,0])# j-1 term
                    A[index0,find.find_index(i,j+1,0)] = 1/(R[i,j,0,1,0] + R[i,j+1,0,1,1])# j+1 term 
                    tempj0 = A[index0,find.find_index(i,j-1,0)] + A[index0,find.find_index(i,j+1,0)]
                                                                    
                    A[index1,find.find_index(i,j-1,c-1)] = 1/(R[i,j,c-1,1,1] + R[i,j-1,c-1,1,0])# j-1 term
                    A[index1,find.find_index(i,j+1,c-1)] = 1/(R[i,j,c-1,1,0] + R[i,j+1,c-1,1,1])# j+1 term 
                    tempj1 = A[index1,find.find_index(i,j-1,c-1)] + A[index0,find.find_index(i,j+1,c-1)]
                    
            
                A[index0,find.find_index(i,j,1)] = 1/(R[i,j,0,2,0] + R[i,j,1,2,2] + R[i,j,1,2,3])# k+1 term
                A[index1,find.find_index(i,j,c-2)] = 1/(R[i,j,c-1,2,1] + R[i,j,c-2,2,0] + R[i,j,c-2,2,1])# k-1 term
                   
            
                coeff0 = tempi0 + tempj0 + A[index0,find.find_index(i,j,1)]
                A[index0,find.find_index(i,j,0)] = -coeff0
        
                coeff1 = tempi1 + tempj1 + A[index1,find.find_index(i,j,c-2)]
                A[index1,find.find_index(i,j,c-1)] = -coeff1        
###############################################################################################################        
# End of Side wall
###############################################################################################################
# Top and bottom wall
        for k in range(0,c):# Top wall
            index0 = find.find_index(i,0,k)
            index1 = find.find_index(i,b-1,k)
            
            if(np.mod(k,2) == 1): # Type 1 Solid node

                if(i == 0):
                    A[index0,find.find_index(i+1,0,k)] = 1/(R[i,0,k,0,0] + R[i+1,0,k,0,1])# i+1 term
                    tempi0 = A[index0,find.find_index(i+1,0,k)]

                    A[index1,find.find_index(i+1,b-1,k)] = 1/(R[i,b-1,k,0,0] + R[i+1,b-1,k,0,1])# i+1 term
                    tempi1 = A[index1,find.find_index(i+1,b-1,k)]
                    
                elif(i == a-1):
                    A[index0,find.find_index(i-1,0,k)] = 1/(R[i,0,k,0,1] + R[i-1,0,k,0,0])# i-1 term
                    tempi0 = A[index0,find.find_index(i-1,0,k)]

                    A[index1,find.find_index(i-1,b-1,k)] = 1/(R[i,b-1,k,0,1] + R[i-1,b-1,k,0,0])# i-1 term
                    tempi1 = A[index1,find.find_index(i-1,b-1,k)]                    
                    
                else:
                    A[index0,find.find_index(i-1,0,k)] = 1/(R[i,0,k,0,1] + R[i-1,0,k,0,0])# i-1 term
                    A[index0,find.find_index(i+1,0,k)] = 1/(R[i,0,k,0,0] + R[i+1,0,k,0,1])# i+1 term
                    tempi0 = A[index0,find.find_index(i+1,0,k)] + A[index0,find.find_index(i-1,0,k)]
                    
                    A[index1,find.find_index(i-1,b-1,k)] = 1/(R[i,b-1,k,0,1] + R[i-1,b-1,k,0,0])# i-1 term
                    A[index1,find.find_index(i+1,b-1,k)] = 1/(R[i,b-1,k,0,0] + R[i+1,b-1,k,0,1])# i+1 term
                    tempi1 = A[index1,find.find_index(i+1,b-1,k)] + A[index1,find.find_index(i-1,b-1,k)]

                if(k == 0):
                    A[index0,find.find_index(i,0,k+1)] = 1/(R[i,0,k,2,0] + R[i,0,k+1,2,1])# k+1 term
                    tempk0 = A[index0,find.find_index(i,0,k+1)]
                    
                    A[index1,find.find_index(i,b-1,k+1)] = 1/(R[i,b-1,k,2,0] + R[i,b-1,k+1,2,1])# k+1 term
                    tempk1 = A[index1,find.find_index(i,b-1,k+1)]
                    
                                     
                elif(k == c-1):
                    A[index0,find.find_index(i,0,k-1)] = 1/(R[i,0,k,2,1] + R[i,0,k-1,2,0])# k-1 term
                    tempk0 = A[index0,find.find_index(i,0,k-1)]
                    
                    A[index1,find.find_index(i,b-1,k-1)] = 1/(R[i,b-1,k,2,1] + R[i,b-1,k-1,2,0])# k-1 term
                    tempk1 = A[index1,find.find_index(i,b-1,k-1)]
                
                    
                else:
                    A[index0,find.find_index(i,0,k-1)] = 1/(R[i,0,k,2,1] + R[i,0,k-1,2,0])# k-1 term
                    A[index0,find.find_index(i,0,k+1)] = 1/(R[i,0,k,2,0] + R[i,0,k+1,2,1])# k+1 term
                    tempk0 = A[index0,find.find_index(i,0,k-1)] + A[index0,find.find_index(i,0,k+1)]
                    
                    A[index1,find.find_index(i,b-1,k-1)] = 1/(R[i,b-1,k,2,1] + R[i,b-1,k-1,2,0])# k-1 term
                    A[index1,find.find_index(i,b-1,k+1)] = 1/(R[i,b-1,k,2,0] + R[i,b-1,k+1,2,1])# k+1 term
                    tempk1 = A[index1,find.find_index(i,b-1,k-1)] + A[index1,find.find_index(i,b-1,k+1)]
                    

                A[index0,find.find_index(i,1,k)] = 1/(R[i,0,k,1,0] + R[i,1,k,1,2] + R[i,1,k,1,3])# j+1 term
                A[index1,find.find_index(i,b-2,k)] = 1/(R[i,b-1,k,1,1] + R[i,b-2,k,1,0] + R[i,b-2,k,1,1])# j-1 term

                

                coeff0 = tempi0 +tempk0 + A[index0,find.find_index(i,1,k)]
                A[index0,find.find_index(i,0,k)] = -coeff0     
                
                coeff1 = tempi1 +tempk1 + A[index1,find.find_index(i,b-2,k)]
                A[index1,find.find_index(i,b-1,k)] = -coeff1
            
            
            
    return(A,C)                