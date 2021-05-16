import numpy as np
from scipy.special import eval_legendre
import math
import random
from porazdelitve import nakljucna,fibonacci

polinom = eval_legendre
cos, acos = math.cos, math.acos
pi = math.pi

def S(data,l):
    N = len(data)
    vsota = 0

    for k in data:
        for t in data:
            produkt = np.dot(k,t)
            vsota += polinom(l,produkt)
    return vsota/N

def varianca(data):
    vsota = 0
    N = len(data)

    for l in range(1,101):
        produkt = polinom(l+1,cos(theta)) - polinom(l-1,cos(theta))
        nov = S(data,l)*produkt**2/(2*l+1)
        vsota += nov
        #print(nov)
    return vsota*N/4

def stevilo(n,kot,porazdelitev):
    m = 0
    for i in range(n):
        so_not = 0
        sredisce = nakljucna(1)[0]
        for vektor in porazdelitev:
            if acos(np.dot(sredisce,vektor)) <= kot: 
                so_not += 1
        m += so_not/n
    return m
    
datoteka = "85.xyz"
data = np.loadtxt(datoteka, skiprows=2, usecols=[1,2,3]) 
#data = nakljucna(85)
theta = pi/20


#### Izpisen sfericni faktor za razlicne l-je
if 0:
    for i in range(1,101):
        print(i,S(data,i),sep='\t')

### Izracunam varianco z metropolisom
if 1:
    n = 1000
    m = stevilo(n,theta,data)
    varianca = math.sqrt(m*(n-m))/n
    print(m,n)
    print(varianca)

### Izracunam varianco z uporabo strukturnega faktorja 
if 0:   # izracun variance.2
    print(varianca(data))

