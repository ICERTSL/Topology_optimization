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
        
    return(((f/8)*(Re-1000)*Pr)/(1+12.7*f2*(Pr2-1)))#return(100)


def Cp(x,y,i,j,k): #finds Cp, x -Temperature , y -pressure i,j,k are node index
    if(np.mod(j,4) == 1 and np.mod(k,2) == 1): # Hot channel
        C=PropsSI('C','T',x,'P',y, I.hot_fluid)                               
    if(np.mod(j,4) == 3 and np.mod(k,2) == 1): # Cold channel
        C=PropsSI('C','T',x,'P',y,I.cold_fluid)
        
    return(C)#return(1211)

def k_liquid(x,y,i,j,k):#finds thermal conductivity of fluid, x -Temperature , y -pressure i,j,k are node index
    if(np.mod(j,4) == 1 and np.mod(k,2) == 1): # Hot channel
        k=PropsSI('L','T',x,'P',y,I.hot_fluid)                               
    if(np.mod(j,4) == 3 and np.mod(k,2) == 1): # Cold channel
        k=PropsSI('L','T',x,'P',y, I.cold_fluid)

    return(k)#return(0.01)

def find_equation_id(index):# finds the indices i,j,k corresponding to an index in the constants matrix C
    k = np.mod(index,I.c)
    j = (index - k)/I.c
    j = np.mod(j,I.b)
    i= (index - I.c*j - k)/(I.c*I.b)

    return(i,j,k)

def find_equation_number(i,j,k): # finds index in Constant Matrix C corresponding to i,j,k where i,j,k is node at which we are doing energy balance
    index = I.b*I.c*i + I.c*j + k
    if(i>I.a-1 or j>I.b-1 or k>I.c-1):
        print("error")
    return(index)

def find_ijk(index):# finds the indices i,j,k corresponding to an index in the constants matrix C
    c = I.c
    b = I.b
    a = I.a
    k = np.mod(index,c)
    j = np.mod((index - k)/c,b)
    if(np.mod(j,4) == 3 and np.mod(k,2) == 1):# cold channel
        info = 2 #yz
        i = (index - c*j - k)/(b*c)
        
    elif(np.mod(j,4) == 1 and np.mod(k,2) == 1):# hot channel
        info = 2 #yz
        i = (index - c*j - k)/(b*c) + 1
    
    else:
        info = 1 #g
        i = (index - c*j - k)/(b*c)    
        

    return(i,j,k)# info = 1 => ijk are for a node in Tg and info = 2 => ijk are for a face in Tyz

def find_index(i,j,k): # finds index in Constant Matrix C corresponding to i,j,k
    c = I.c
    b = I.b
    a = I.a
    for index in range(0,a*b*c):
        x,y,z = find_ijk(index)
        if(x==i and y== j and z==k):
            return(index)


def Rwall():
    Rwx = I.delta_x/(2*I.k_solid*I.d*I.d)
    Rwy = 1/(2*I.k_solid*I.delta_x)
    Rwz= 1/(2*I.k_solid*I.delta_x)
    return(Rwx,Rwy,Rwz)    
