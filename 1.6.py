from math import sqrt
from random import randint
import numpy as np
import matplotlib.pyplot as plt

time = 1000
dt = 0.001
n = 50  #количество частиц
c = 0.2 #концентрация
T = 2 #температура
massa = 1
velocity = int(0.7*sqrt(3*T/massa) + 0.5) #предел начальной скорости
#константы взаимодействия
E = 1 
radius = 1 #сигма
#границы куба
wall = int(0.5*(pow(n/c, 1/3)) + 0.5)


#в функцию подается квадрат расстояния между частицами, а возвращается модуль силы
def force(r): 
    f=(radius**2/r)**3
    result = -24*E*f*(2*f-1)
    return result

def distancesquare(coord1, coord2):
    return (coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2

def vector(coord1, coord2):
    result = [coord2[0] - coord1[0], coord2[1] - coord1[1], coord2[2] - coord1[2]]
    for cel in range(3):
        subres = abs(result[cel])
        if subres > wall:
            while subres > 2*wall:
                subres += -2*wall
            result[cel] *= -(2*wall - subres)/abs(result[cel])
    return result

#в функцию подается номер частицы, а возвращается ее ускорение
def acceleration(x):
    result = [0, 0, 0]
    for zel in range(n):
        if (zel==x):
            pass
        else:
            radius_vector = vector(coordinate[x], coordinate[zel])
            dist = radius_vector[0]**2 + radius_vector[1]**2 + radius_vector[2]**2
            force_all = force(dist)/(dist*massa)
            for zul in range(3):
                result[zul] += force_all*radius_vector[zul]
    return result

def kineticenergy():
    result = 0.
    for al in range(n):
        result += velocitys[al][0]**2 + velocitys[al][1]**2 + velocitys[al][2]**2
    return result*massa/2

def potentionalenergy():
    result = 0.
    for el in range(n):
        for ol in range(el+1, n):
            radius_vector = vector(coordinate[el], coordinate[ol])
            dist = radius_vector[0]**2 + radius_vector[1]**2 + radius_vector[2]**2
            f = (radius**2/dist)**3
            result += f*(f - 1)
    return 4*E*result

coordinate = np.array([(0, 0, 0)]*n, float)
for el in range(n):
    coordinate[el] = (0.001*randint(-wall*1000, wall*1000), 0.001*randint(-wall*1000, wall*1000), 0.001*randint(-wall*1000, wall*1000))
control = False
while not control:
    control = True
    for ul in range(n):
        for el in range(n):
            if ul==el:
                pass
            elif distancesquare(coordinate[ul], coordinate[el]) < radius**2:
                control = False
                for ol in range(3):
                    coordinate[el][ol] = 0.001*randint(-wall*1000, wall*1000)
             
velocitys = np.array([(0, 0, 0)]*n, float)
for el in range(n):
    velocitys[el] = (randint(-velocity, velocity), randint(-velocity, velocity), randint(-velocity, velocity))

acceleration1 = np.array([(0, 0, 0)]*n, float)
acceleration2 = np.array([(0, 0, 0)]*n, float)

for el in range(n):
    acceleration1[el] = acceleration(el)

kinenergy = np.array([0]*time, float)
totalenergy = np.array([0]*time, float)
potenergy = np.array([0]*time, float)

for ul in range(time):
    for el in range(n):
        for ol in range(3):
            coordinate[el][ol] += velocitys[el][ol]*dt + 0.5*acceleration1[el][ol]*dt**2
    for el in range(n):
         acceleration2[el] = acceleration(el)
    for el in range(n):
        for ol in range(3):
            velocitys[el][ol] += 0.5*(acceleration1[el][ol] + acceleration2[el][ol])*dt
    for el in range(n):
        acceleration1 = acceleration2
    kinenergy[ul] = kineticenergy()
    potenergy[ul] = potentionalenergy()
    totalenergy[ul] = potenergy[ul] + kinenergy[ul]

print(2*np.sum(kinenergy)/(3*n*time))
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