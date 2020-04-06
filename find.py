# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 01:42:50 2020

@author: HP
"""
import numpy as np
import matplotlib.pyplot as plt 
from CoolProp.CoolProp import PropsSI
import inputs as I
import initialize


def enthalpy(Pg,Tg):# generates Nodal enthalpy matrix from Pg and Tg
    a = I.a
    b = I.b
    c = I.c
    hg = np.zeros([a,b,c])
    
    for i in range(0,a +1): 
        for j in range(0,b):
            for k in range(0,c):
                if(np.mod(j,4) == 1 and np.mod(k,2) == 1): # Hot channel
                    hg[i,j,k] = PropsSI('H','T',Tg[i,j,k],'P',Pg[i,j,k], I.hot_fluid)
                                
                if(np.mod(j,4) == 3 and np.mod(k,2) == 1): # Cold channel
                    hg[i,j,k] = PropsSI('H','T',Tg[i,j,k],'P',Pg[i,j,k], I.cold_fluid) 

    return(hg)

def entropy(Pg,Tg):# generates Nodal entropy matrix from Pg and Tg
    a = I.a
    b = I.b
    c = I.c
    sg = np.zeros([a,b,c])
    
    for i in range(0,a +1): 
        for j in range(0,b):
            for k in range(0,c):
                if(np.mod(j,4) == 1 and np.mod(k,2) == 1): # Hot channel
                    sg[i,j,k] = PropsSI('S','T',Tg[i,j,k],'P',Pg[i,j,k], I.hot_fluid)
                                
                if(np.mod(j,4) == 3 and np.mod(k,2) == 1): # Cold channel
                    sg[i,j,k] = PropsSI('S','T',Tg[i,j,k],'P',Pg[i,j,k], I.cold_fluid)    

    return(sg)

def pressure_drop(x,y,fluid):

    mu = PropsSI('V','T',x,'P',y, fluid)
    density=PropsSI('D','T',x,'P',y, fluid)
    Re=I.m_channel/(I.d*mu)
    deltap=friction_factor(Re)*(0.5*I.m_channel*I.m_channel*I.delta_x/(density*(np.power(I.d,5))))  
          
    return(deltap)

def friction_factor(x): # x is Reynolds number
	f=np.power((0.79*(np.log(x))-1.64),-2) # use blassius eqn?
	return(f)

#Finds the Nusselt number given the Reynolds number
def Nu(x,y,i,j,k): #x - temperature, y - pressure ,i,j,k are node index
  
    if(np.mod(j,4) == 1 and np.mod(k,2) == 1): # Hot channel
        
        mu=PropsSI('V','T',x,'P',y,I.hot_fluid)      
        cp=PropsSI('C','T',x,'P',y, I.hot_fluid)
        k=PropsSI('L','T',x,'P',y,I.hot_fluid)
        Re=I.m_channel/(I.d*mu)
        f = friction_factor(Re)
        # Pr=refpropm('^','T',T,'P',P,I.hot_fluid);
        Pr = (mu*cp)/k
        Pr2 = np.power(Pr, 0.66666)
        f2 = f/8
        f2 = np.power(f2,0.5)                                
    elif(np.mod(j,4) == 3 and np.mod(k,2) == 1): # Cold channel


        mu=PropsSI('V','T',x,'P',y, I.cold_fluid)
        cp=PropsSI('C','T',x,'P',y, I.cold_fluid)
        k=PropsSI('L','T',x,'P',y, I.cold_fluid)
        Re=I.m_channel/(I.d*mu)
        f = friction_factor(Re)
        # Pr=refpropm('^','T',T,'P',P, I.hot_fluid);
        Pr = (mu*cp)/k
        Pr2 = np.power(Pr, 0.66666)
        f2 = f/8
        f2 = np.power(f2,0.5)
        
    return(((f/8)*(Re-1000)*Pr)/(1+12.7*f2*(Pr2-1)))


def Cp(x,y,i,j,k): #finds Cp, x -Temperature , y -pressure i,j,k are node index
    if(np.mod(j,4) == 1 and np.mod(k,2) == 1): # Hot channel
        C=PropsSI('C','T',x,'P',y, I.hot_fluid)                               
    if(np.mod(j,4) == 3 and np.mod(k,2) == 1): # Cold channel
        C=PropsSI('C','T',x,'P',y,I.cold_fluid)
    return(C)

def k_liquid(x,y,i,j,k):#finds thermal conductivity of fluid, x -Temperature , y -pressure i,j,k are node index
    if(np.mod(j,4) == 1 and np.mod(k,2) == 1): # Hot channel
        k=PropsSI('L','T',x,'P',y,I.hot_fluid)                               
    if(np.mod(j,4) == 3 and np.mod(k,2) == 1): # Cold channel
        k=PropsSI('L','T',x,'P',y, I.cold_fluid)

    return(k)

def find_ijk(index):# finds the indices i,j,k corresponding to an index in the constants matrix C
    i = np.mod(index,I.a)
    j = (index - i)/I.a
    j = np.mod(j,I.b)
    k= (index - I.a*j - i)/(I.a*I.b)

    return(i,j,k)

def find_index(i,j,k): # finds index in Constant Matrix C corresponding to i,j,k
    index = I.a*I.b*k + I.a*j + i

    return(index)
    