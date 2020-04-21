# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 20:25:21 2020

@author: HP
"""
import numpy as np
import matplotlib.pyplot as plt 
from CoolProp.CoolProp import PropsSI
import inputs as I
import initialize
import find
import BC

def Rc(x,y,i,j,k):# finds convective resistance x- temp, y - pressure
    nu = find.Nu(x,y,i,j,k)
    kliq = find.k_liquid(x,y,i,j,k)
    htc = (kliq*nu)/I.d
    area = I.d*I.delta_x
    resistance = 1/(htc*area)
    return(resistance)

def adiabatic_wall_BC(A,C,Rwx,Rwy,Rwz,Tg,Tyz,Pg,P):
    a = I.a
    b = I.b
    c = I.c
    
    for i in range(0,a):
###############################################################################################################
#side Walls
                
        for j in range(0,b):# Side Walls
            index0 = find.find_equation_number(i,j,0)
            index1 = find.find_equation_number(i,j,c-1)       
            index0 = int(index0)
            index1 = int(index1)

            
            if(np.mod(j,2) == 0): # Type 2 Solid node
                C0_0 = 0
                C0_1 = 0
                C01 = 1/(2*Rwx)
                C02 = 1/(2*Rwx)
                C11 = 1/(2*Rwy)
                C12 = 1/(2*Rwy)
                C21 = 1/(2*Rwz)
                C22 = 1/(2*Rwz)
                
                if(i == 0): # Nodes along y axis at x = 0
                    C0_0 = C0_0 -C02 -C22 
                    A[index0,find.find_index(i+1,j,0)] = C02# i+1 term
                    A[index0,find.find_index(i,j,1)] = C22# k+1 term

                    C0_1 = C0_1 -C02 -C21 
                    A[index1,find.find_index(i+1,j,c-1)] = C02# i+1 term
                    A[index1,find.find_index(i,j,c-2)] = C21# k-1 term

                    
                elif(i == a-1):#Nodes along y axis at x = L
                    C0_0 = C0_0 -C01 -C22 
                    A[index0,find.find_index(i-1,j,0)] = C02# i-1 term
                    A[index0,find.find_index(i,j,1)] = C22# k+1 term

                    C0_1 = C0_1 -C01 -C21 
                    A[index1,find.find_index(i-1,j,c-1)] = C02# i-1 term
                    A[index1,find.find_index(i,j,c-2)] = C21# k-1 term
                    
                else: # Nodes in interior of XY plane 
                    C0_0 = C0_0 -C01 - C02 -C22 
                    A[index0,find.find_index(i-1,j,0)] = C02# i-1 term
                    A[index0,find.find_index(i+1,j,0)] = C02# i+1 term
                    A[index0,find.find_index(i,j,1)] = C22# k+1 term

                    C0_1 = C0_1 -C01 -C02 -C21 
                    A[index1,find.find_index(i-1,j,c-1)] = C02# i-1 term
                    A[index1,find.find_index(i+1,j,c-1)] = C02# i+1 term
                    A[index1,find.find_index(i,j,c-2)] = C21# k-1 term
                    
                if(j == 0):#Nodes along x axis at y = 0
                    C0_0 = C0_0 - C12 
                    A[index0,find.find_index(i,j+1,0)] = C12# j+1 term
                    
                    C0_1 = C0_1 - C12
                    A[index1,find.find_index(i,j+1,c-1)] = C12# j+1 term

                    
                elif(j == b-1):#Nodes along x axis at y = ymax
                    C0_0 = C0_0 - C11 
                    A[index0,find.find_index(i,j-1,0)] = C11# j-1 term
                    
                    C0_1 = C0_1 - C11
                    A[index1,find.find_index(i,j-1,c-1)] = C11# j-1 term
                
                else: #interior of XY plane
               
                    C0_0 = C0_0 - C11 -C12
                    A[index0,find.find_index(i,j-1,0)] = C11# j-1 term
                    A[index0,find.find_index(i,j+1,0)] = C12# j+1 term
                    
                    C0_1 = C0_1 - C11 - C12
                    A[index1,find.find_index(i,j-1,c-1)] = C11# j-1 term
                    A[index1,find.find_index(i,j+1,c-1)] = C12# j+1 term
             
                                
                A[index0,find.find_index(i,j,0)] = C0_0
                A[index1,find.find_index(i,j,c-1)] = C0_1
        
####################################################################################################                
            if(np.mod(j,2) == 1): # Type 3 Solid node
                C0_0 = 0
                C0_1 = 0
                C01 = 1/(2*Rwx)
                C02 = 1/(2*Rwx)
                C11 = 1/(2*Rwy)
                C12 = 1/(2*Rwy)
                C21 = 0.5/(Rwz + Rc(Tg[i,j,c-2],Pg[i,j,c-2],i,j,c-2))
                C22 = 0.5/(Rwz + Rc(Tg[i,j,1],Pg[i,j,1],i,j,1))
        
               
                if(i == 0): 
                    C0_0 = C0_0 -C02 - 2*C22
                    C0_1 = C0_1 -C02 - 2*C21
                    A[index0,find.find_index(i+1,j,0)] = C02# Tg i+1,j,k term
                    A[index1,find.find_index(i+1,j,c-1)] = C02# Tg i+1,j,k term
                    
                    if(np.mod(j,4) == 1):#Along hot channel -> i,j,k+1 is hot inlet bc
                        A[index0,find.find_index(i+1,j,1)] = C22#i+1, k+1 term
                        C[index0,0] = -I.T_hot_in*C22 #Hot in BC

                        A[index1,find.find_index(i+1,j,c-2)] = C21# i+1,k-1 term
                        C[index1,0] = -I.T_hot_in*C21 #Hot in BC
                        
                    elif(np.mod(j,4) == 3):#Along Cold channel
                        A[index0,find.find_index(i+1,j,1)] = C22# k+1 term
                        A[index0,find.find_index(i,j,1)] = C22# k+1 term                  

                        A[index1,find.find_index(i+1,j,c-2)] = C21# k-1 term
                        A[index1,find.find_index(i,j,c-2)] = C21# k-1 term   
                   
                elif(i == a-1):#Nodes along y axis at x = L
                    C0_0 = C0_0 -C01 - 2*C22
                    C0_1 = C0_1 -C01 - 2*C21
                    A[index0,find.find_index(i-1,j,0)] = C01# Tg i-1,j,k term
                    A[index1,find.find_index(i-1,j,c-1)] = C01# Tg i+1,j,k term
                    
                    if(np.mod(j,4) == 1):#Along hot channel
                        
                        A[index0,find.find_index(i+1,j,1)] = C22# k+1 term
                        A[index0,find.find_index(i,j,1)] = C22# k+1 term                  

                        A[index1,find.find_index(i+1,j,c-2)] = C21# k-1 term
                        A[index1,find.find_index(i,j,c-2)] = C21# k-1 term                         

                        
                    elif(np.mod(j,4) == 3):#Along Cold channel
 
                    
                        A[index0,find.find_index(i,j,1)] = C22#, k+1 term
                        C[index0,0] = -I.T_cold_in*C22 #Cold in BC

                        A[index1,find.find_index(i,j,c-2)] = C21# i,k-1 term
                        C[index1,0] = -I.T_cold_in*C21 #Cold in BC                    
                    
                else: # Nodes in interior of XY plane 
                    C0_0 = C0_0 -C01 - C02 - 2*C22
                    A[index0,find.find_index(i+1,j,0)] = C02# Tg i+1,j,k term
                    A[index0,find.find_index(i-1,j,0)] = C01# Tg i-1,j,k term
                    A[index0,find.find_index(i+1,j,1)] = C22#Tyz i+1, k+1 term
                    A[index0,find.find_index(i,j,1)] = C22#Tyz i, k+1 term    

                    C0_1 = C0_1 -C01 - C02 - 2*C21
                    A[index1,find.find_index(i+1,j,c-1)] = C02# Tg i+1,j,k term
                    A[index1,find.find_index(i-1,j,c-1)] = C01# Tg i-1,j,k term
                    A[index1,find.find_index(i+1,j,c-2)] = C21#Tyz i+1, k-1 term
                    A[index1,find.find_index(i,j,c-2)] = C21  #Tyz i, k-1 term


                    
                if(j == 0):#Nodes along x axis at y = 0
                    C0_0 = C0_0 - C12 
                    A[index0,find.find_index(i,j+1,0)] = C12# j+1 term                   
                    C0_1 = C0_1 - C12
                    A[index1,find.find_index(i,j+1,c-1)] = C12# j+1 term

                    
                elif(j == b-1):#Nodes along x axis at y = ymax
                    C0_0 = C0_0 - C11 
                    A[index0,find.find_index(i,j-1,0)] = C11# j-1 term                    
                    C0_1 = C0_1 - C11
                    A[index1,find.find_index(i,j-1,c-1)] = C11# j-1 term
                
                else: #interior of XY plane
               
                    C0_0 = C0_0 - C11 -C12
                    A[index0,find.find_index(i,j-1,0)] = C11# j-1 term
                    A[index0,find.find_index(i,j+1,0)] = C12# j+1 term
                    
                    C0_1 = C0_1 - C11 - C12
                    A[index1,find.find_index(i,j-1,c-1)] = C11# j-1 term
                    A[index1,find.find_index(i,j+1,c-1)] = C12# j+1 term
             
            
                A[index0,find.find_index(i,j,0)] = C0_0
                A[index1,find.find_index(i,j,c-1)] = C0_1   
                
 
###############################################################################################################        
# End of Side wall
###############################################################################################################
# Top and bottom wall
        for k in range(0,c):# Top wall
            index0 = find.find_equation_number(i,0,k)
            index1 = find.find_equation_number(i,b-1,k)

        
            
            if(np.mod(k,2) == 1): # Type 1 Solid node

                C0_0 = 0
                C0_1 = 0
                C01 = 1/(2*Rwx)
                C02 = 1/(2*Rwx)
                C11 = 0.5/(Rwy + Rc(Tg[i,b-2,k],Pg[i,b-2,k],i,b-2,k))
                C12 = 0.5/(Rwy + Rc(Tg[i,1,k],Pg[i,1,k],i,1,k))
                C21 = 1/(2*Rwz)
                C22 = 1/(2*Rwz)
                
                if(i == 0):
                    C0_0 = C0_0 - C02 - 2*C12
                    A[index0,find.find_index(i+1,0,k)] = C02# Tg i+1,j,k term
                    A[index0,find.find_index(i+1,1,k)] = C12# Tyz i+1,j+1
                    C[index0,0] = -I.T_hot_in*C12# Hot in BC

                    C0_1 = C0_1 - C02 - 2*C11
                    A[index1,find.find_index(i+1,b-1,k)] = C02# Tg i+1,j,k term
                    A[index1,find.find_index(i,b-2,k)] = C11# Tyz i,j-1
                    A[index1,find.find_index(i+1,b-2,k)] = C11# Tyz i+1,j-1
               
                    
                elif(i == a-1):
                    C0_0 = C0_0 - C01 - 2*C12
                    A[index0,find.find_index(i-1,0,k)] = C01# Tg i-1,j,k term
                    A[index0,find.find_index(i+1,1,k)] = C12# Tyz i+1,j+1
                    A[index0,find.find_index(i,1,k)] = C12# Tyz i,j+1


                    C0_1 = C0_1 - C01 - 2*C11
                    A[index1,find.find_index(i-1,b-1,k)] = C01# Tg i-1,j,k term
                    A[index1,find.find_index(i,b-2,k)] = C11# Tyz i,j-1
                    C[index1,0] = -I.T_cold_in*C11# Cold in BC
                    
                else:
                    C0_0 = C0_0 - C01 - C02 - 2*C12
                    A[index0,find.find_index(i-1,0,k)] = C01# Tg i-1,j,k term
                    A[index0,find.find_index(i+1,0,k)] = C02# Tg i+1,j,k term
                    A[index0,find.find_index(i+1,1,k)] = C12# Tyz i+1,j+1
                    A[index0,find.find_index(i,1,k)] = C12# Tyz i,j+1

                    C0_1 = C0_1 - C01 - C02 - 2*C11
                    A[index1,find.find_index(i+1,b-1,k)] = C02# Tg i+1,j,k term
                    A[index1,find.find_index(i,b-2,k)] = C11# Tyz i,j-1
                    A[index1,find.find_index(i+1,b-2,k)] = C11# Tyz i+1,j-1
                    A[index1,find.find_index(i-1,b-1,k)] = C01# Tg i-1,j,k term
                    
                if(k == 0):
                    C0_0 = C0_0 - C22
                    A[index0,find.find_index(i,0,k+1)] = C22
                    
                    C0_1 = C0_1 - C22
                    A[index1,find.find_index(i,b-1,k+1)] = C22
                    
                    
                elif(k == c-1):
                    C0_0 = C0_0 - C21
                    A[index0,find.find_index(i,0,k-1)] = C21
                    
                    C0_1 = C0_1 - C21
                    A[index1,find.find_index(i,b-1,k-1)] = C21                
                    
                else:
                    C0_0 = C0_0 - C21 - C22
                    A[index0,find.find_index(i,0,k-1)] = C21
                    A[index0,find.find_index(i,0,k+1)] = C22
                    
                    C0_1 = C0_1 - C21 - C22
                    A[index1,find.find_index(i,b-1,k-1)] = C21
                    A[index1,find.find_index(i,b-1,k+1)] = C22

                A[index0,find.find_index(i,0,k)] = C0_0
                A[index1,find.find_index(i,b-1,k)] = C0_1    
                
            if(np.mod(k,2) == 0): # Type 2 Solid node

                C0_0 = 0
                C0_1 = 0
                C01 = 1/(2*Rwx)
                C02 = 1/(2*Rwx)
                C11 = 1/(2*Rwy)
                C12 = 1/(2*Rwy)
                C21 = 1/(2*Rwz)
                C22 = 1/(2*Rwz)
                
                if(i == 0):
                    C0_0 = C0_0 - C02 - C12
                    A[index0,find.find_index(i+1,0,k)] = C02# Tg i+1,j,k
                    A[index0,find.find_index(i,1,k)] = C12 # Tg i, j+1,k
                    
                    C0_1 = C0_1 - C02 - C11
                    A[index1,find.find_index(i+1,b-1,k)] = C02# Tg i+1,j,k
                    A[index1,find.find_index(i,b-2,k)] = C11 # Tg i, j-1,k                    
                    
                 
                    
                elif(i == a-1):
                    C0_0 = C0_0 - C01 - C12
                    A[index0,find.find_index(i-1,0,k)] = C01# Tg i-1,j,k
                    A[index0,find.find_index(i,1,k)] = C12 # Tg i, j+1,k
                    
                    C0_1 = C0_1 - C01 - C11
                    A[index1,find.find_index(i-1,b-1,k)] = C01# Tg i-1,j,k
                    A[index1,find.find_index(i,b-2,k)] = C11 # Tg i, j-1,k 
                    
                else:
                    C0_0 = C0_0 - C01 -C02 - C12
                    A[index0,find.find_index(i+1,0,k)] = C02# Tg i+1,j,k
                    A[index0,find.find_index(i-1,0,k)] = C01# Tg i-1,j,k
                    A[index0,find.find_index(i,1,k)] = C12 # Tg i, j+1,k
                    
                    C0_1 = C0_1 - C01 - C02 - C11
                    A[index1,find.find_index(i+1,b-1,k)] = C02# Tg i+1,j,k
                    A[index1,find.find_index(i-1,b-1,k)] = C01# Tg i-1,j,k
                    A[index1,find.find_index(i,b-2,k)] = C11 # Tg i, j-1,k 


                if(k == 0):
                    C0_0 = C0_0 - C22
                    A[index0,find.find_index(i,0,k+1)] = C22# Tg i,j,k+1
                    
                    C0_1 = C0_1 - C22
                    A[index1,find.find_index(i,b-1,k+1)] = C22# Tg i,j,k+1                    
                                                         
                elif(k == c-1):
                    C0_0 = C0_0 - C21
                    A[index0,find.find_index(i,0,k-1)] = C21# Tg i,j,k-1
                    
                    C0_1 = C0_1 - C21
                    A[index1,find.find_index(i,b-1,k-1)] = C21# Tg i,j,k-1                  
                    
                else:
                    C0_0 = C0_0 - C21 -C22
                    A[index0,find.find_index(i,0,k-1)] = C21# Tg i,j,k-1
                    A[index0,find.find_index(i,0,k+1)] = C22# Tg i,j,k+1
                    
                    C0_1 = C0_1 - C21-C22
                    A[index1,find.find_index(i,b-1,k-1)] = C21# Tg i,j,k-1
                    A[index1,find.find_index(i,b-1,k+1)] = C22# Tg i,j,k+1

                A[index0,find.find_index(i,0,k)] = C0_0# Tg i,j,k
                A[index1,find.find_index(i,b-1,k)] = C0_1# Tg i,j,k
            
    return(A,C)                