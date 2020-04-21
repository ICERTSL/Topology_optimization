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
# print(Tg[:,:,1]- 273.15)
# print(".................")
# print(Tg[:,:,2]- 273.15)

# print(".................")
# print("Tyz")
# print(".................")
# print(Tyz[:,:,1]- 273.15)
# print(".................")
# print(Tyz[:,:,1]- 273.15)
# print(".................")
# print(Tyz[:,:,2]- 273.15)

Pg, P,Pgstat,Pstat = initialize.pressure()
# print("Pg")
# print(".................")
# print(Pg[:,:,1]/1e5)

# print(".................")
# print("P")
# print(".................")
# print(P[:,:,1]/1e5)

Pg, P = update.pressure(P,Pg,Pgstat,Pstat,Tg)

# print("Pg")
# print(".................")
# print(Pg[:,:,1]/1e5)

# print(".................")
# print("P")
# print(".................")
# print(P[:,:,1]/1e5)

# print(".................")
# print("Rhot")
# R = update.resistance(P,Tyz,Pg,Tg)

ordinate = I.delta_x*np.arange(0,I.a+1)
ordinateg = I.delta_x*np.arange(0,I.a) + I.delta_x*0.5
Rwx,Rwy,Rwz = find.Rwall()
res = 10
count = 0
while(res>1e-3):
    print("iteration = ",count)
    A,C = update.linear_equation_system(Rwx,Rwy,Rwz,Tg,Tyz,Pg,P)
    sol = np.linalg.inv(A)
    target = np.matmul(sol,C)
    Tg,Tyz = update.temperature(target,Tg,Tyz)
    Pg, P = update.pressure(P,Pg,Tg,Pstat,Tg)
  
    
    plt.plot(ordinate,Tyz[:,1,1] -273)
    plt.plot(ordinate,Tyz[:,3,1] -273)


    
    A,C = update.linear_equation_system(Rwx,Rwy,Rwz,Tg,Tyz,Pg,P)
    sol = np.linalg.inv(A)
    target = np.matmul(sol,C)
    sumdhhot = 0
    sumdhcold =0
    for height in range(2,I.b -1,4):
        print("height = ",height)
        for width in range(1,I.c,2):
            print("width = ",width)
            sumdhhot = sumdhhot + abs(PropsSI("H","T",Tyz[0,height -1,width],"P",P[0,height-1,width],"co2") - PropsSI("H","T",Tyz[I.a,height -1,width],"P",P[I.a,height -1,width],"co2"))
            sumdhcold =sumdhcold+ abs(PropsSI("H","T",Tyz[0,height +1,width],"P",P[0,height +1,width],"co2") - PropsSI("H","T",Tyz[I.a,height +1,width],"P",P[I.a,height +1,width],"co2"))    


    res = abs((sumdhhot/sumdhcold) - 1)
    print("res = ",res)
    count = count + 1

plt.show()


for i in range(1,I.b,2):
    plt.plot(ordinate,Tyz[:,i,1] -273)
    print(i)
   
plt.show()    


   

# print("a-1,2,1")
# print(A[find.find_equation_number(I.a-1,3,0),:])
# print(C[find.find_equation_number(I.a-1,3,0)])



# x = 4
# count = 0
# for i in range(0,np.size(A[x,:])):
#     if(A[x,i] != 0):
#         print(i)
#         print(find.find_ijk(i))
#         print("coeff",A[x,i])
#         print("Val = ",target[i])
#         count = count + 1
# print(C[x])
# print("count =",count)


# for i in range(0,np.size(target)):
#     if (target[i] == np.max(target)):
#         print(i)
#         print(find.find_ijk(i))
#         print(np.max(target) - 273)
    
    
# print(target -273)


