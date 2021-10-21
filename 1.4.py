import math as m
from random import*
import numpy as np
import matplotlib.pyplot as plt

time = 10
dt = 0.001
n = 50  #кол-во частиц
c = 0.2 #концетрация
T = 2 #температура
massa = 1
velocity = int(1.6*m.sqrt(3*n*T/massa) + 1) #предел начальной скорости
#константы взаимодействия
E = 1 
D = 1
K = 1
#сигма
radius = 1
#границы куба
wall = pow(n/c, 1./3.)
wallX1 = int(-wall/2 - 1)
wallY1 = int(-wall/2 - 1)
wallZ1 = int(-wall/2 - 1)
wallX2 = int(wall/2 + 1)
wallY2 = int(wall/2 + 1)
wallZ2 = int(wall/2 + 1)
#в функцию подается квадрат расстояния между частицами, а возвращается модуль силы
def force(r):
    f=(radius**2/r)**3
    return -24*E*f*(2*f-1)/m.sqrt(r)

def distancesquare(el1, el2):
    return (Xcoordinate[el1] - Xcoordinate[el2])**2 + (Ycoordinate[el1] - Ycoordinate[el2])**2 + (Zcoordinate[el1] - Zcoordinate[el2])**2 

#в функцию подается номер частицы, а возвращается ее ускорение
def Xacceleration(x):
    result = (9*D/radius)*(radius/(wallX1 - Xcoordinate[x]))**10 - (9*D/radius)*(radius/(wallX2 - Xcoordinate[x]))**10
    for el in range(n):
        if (el==x):
            result += 0
        else:
            d = distancesquare(el,x)
            result += force(d)*(Xcoordinate[x] - Xcoordinate[el])/m.sqrt(d) 
    return result/massa

def Yacceleration(x):
    result = (9*D/radius)*(radius/(wallY1 - Ycoordinate[x]))**10 - (9*D/radius)*(radius/(wallY2 - Ycoordinate[x]))**10
    for el in range(n):
        if (el==x):
            result += 0
        else:
            d = distancesquare(el,x)
            result += force(d)*(Ycoordinate[x] - Ycoordinate[el])/m.sqrt(d)
    return result/massa

def Zacceleration(x):
    result = (9*D/radius)*(radius/(wallZ1 - Zcoordinate[x]))**10 - (9*D/radius)*(radius/(wallZ2 - Zcoordinate[x]))**10
    for el in range(n):
        if (el==x):
            result += 0
        else:
            d = distancesquare(el,x)
            result += force(d)*(Zcoordinate[x] - Zcoordinate[el])/m.sqrt(d)
    return result/massa

def kineticenergy():
    result = 0.
    for el in range(n):
        result += Xvelocity[el]**2 + Yvelocity[el]**2
    return result*massa/2

def potentionalenergy():
    result = 0.
    for el in range(n):
        for ol in range(el+1, n):
            f = (radius**2/distancesquare(el, ol))**3
            result += f*(f - 1)
    result = 4*E*result

    for el in range(n):
        r1 = abs(radius/(Xcoordinate[el] - wallX1))
        r2 = abs(radius/(Xcoordinate[el] - wallX2))
        r3 = abs(radius/(Ycoordinate[el] - wallY1)) 
        r4 = abs(radius/(Ycoordinate[el] - wallY2))
        r5 = abs(radius/(Zcoordinate[el] - wallZ1))
        r6 = abs(radius/(Zcoordinate[el] - wallZ2))
        result += D*(r1**9 + r2**9 + r3**9 + r4**9 + r5**2 + r6**2) 
    return result

Xcoordinate = [randint(wallX1*100, wallX2*100)*0.009 for i in range(n)]
Ycoordinate = [randint(wallY1*100, wallY2*100)*0.009 for i in range(n)]
Zcoordinate  = [randint(wallZ1*100, wallZ2*100)*0.009 for i in range(n)]
Xvelocity = [randint(-velocity, velocity) for i in range(n)]
Yvelocity = [randint(-velocity, velocity) for i in range(n)]
Zvelocity = [randint(-velocity, velocity) for i in range(n)]
Xacceleration1 = [0]*n
Yacceleration1 = [0]*n
Zacceleration1 = [0]*n
Xacceleration2 = [0]*n
Yacceleration2 = [0]*n
Zacceleration2 = [0]*n

for el in range(n):
    Xacceleration1[el] = Xacceleration(el)
    Yacceleration1[el] = Yacceleration(el)
    Zacceleration1[el] = Zacceleration(el)

kinenergy = [0]*time
totalenergy = [0]*time
potenergy = [0]*time

for ul in range(time):
    for el in range(n):
        Xcoordinate[el] += Xvelocity[el]*dt + 0.5*Zacceleration1[el]*dt**2
        Ycoordinate[el] += Yvelocity[el]*dt + 0.5*Yacceleration1[el]*dt**2
        Zcoordinate[el] += Zvelocity[el]*dt + 0.5*Zacceleration1[el]*dt**2
    for el in range(n):
        Xacceleration2[el] = Xacceleration(el)
        Yacceleration2[el] = Yacceleration(el)
        Zacceleration2[el] = Zacceleration(el)
    for el in range(n):
        Xvelocity[el] += 0.5*(Xacceleration1[el] + Xacceleration2[el])*dt
        Yvelocity[el] += 0.5*(Yacceleration1[el] + Yacceleration2[el])*dt
        Zvelocity[el] += 0.5*(Zacceleration1[el] + Zacceleration2[el])*dt
    for el in range(n):
        Xacceleration1[el] = Xacceleration2[el]
        Yacceleration1[el] = Yacceleration2[el]
        Zacceleration1[el] = Zacceleration2[el]
    kinenergy[ul] = kineticenergy()
    potenergy[ul] = potentionalenergy()
    totalenergy[ul] = potenergy[ul] + kinenergy[ul]

kinenergy  = np.array(kinenergy)
potenergy  = np.array(potenergy)
totalenergy = np.array(totalenergy)
x = np.arange(0, time, 1)

plt.figure()
plt.plot(x, kinenergy[x])
plt.xlabel(r'$time$')
plt.ylabel(r'$Kinetic Energy$')
plt.grid(True)

plt.figure()
plt.plot(x, totalenergy[x])
plt.xlabel(r'$time$')
plt.ylabel(r'$Total Energy$')
plt.grid(True)

plt.figure()
plt.plot(x, potenergy[x])
plt.xlabel(r'$time$')
plt.ylabel(r'$Potential Energy$')
plt.grid(True)

plt.show()