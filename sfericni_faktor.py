import numpy as np
from numpy import linalg as LA
from scipy.special import eval_legendre
import math
import random
from porazdelitve import nakljucna,fibonacci
from threading import Thread
from time import sleep, perf_counter

legendre = eval_legendre
cos, acos = math.cos, math.acos
pi = math.pi

### Strukturni faktor kot vsota legendrovih polinomov
def S(data,l):
    N = len(data)
    vsota = 0

    for k in data:
        for t in data:
            produkt = np.dot(k,t)/(LA.norm(k)*LA.norm(t))
            vsota += legendre(l,produkt)
    return vsota/N

### Varianca kot vsota razlik legendrovih polinomov
def varianca(data,kot):
    vsota = 0
    N = len(data)

    for l in range(1,101):
        produkt = legendre(l+1,cos(kot)) - legendre(l-1,cos(kot))
        nov = S(data,l)*produkt**2/(2*l+1)
        vsota += nov
        #print(nov)
    print(kot/pi, vsota*N/4)
    return

### Steje ko je razdalja po krogu manjsa od kota kapice
def stevilo(n,kot,porazdelitev):
    m = np.zeros(n)
    for i in range(n):
        so_not = 0
        sredisce = nakljucna(1)[0]
        t = 0
        for vektor in porazdelitev:
            t += 1
            produkt = np.dot(sredisce,vektor)
            if produkt <= cos(kot) and produkt >= 0:
                so_not += 1
        m[i] += so_not
    return m

### Porazdelitev — random/fibonacci/Thomson ---> lahko bi se kej

#datoteka = "fibonacci/932.txt"
#data = np.loadtxt(datoteka, usecols=[0,1,2]) 

datoteka = "thomson/10.xyz"
data = np.loadtxt(datoteka, skiprows=2, usecols=[1,2,3]) 

#data = nakljucna(20)

### Interval dolg pi/2 na a delov v th
a = 32
th = np.linspace(0, pi/2, a)

### Izpise sfericni faktor za razlicne l-je
if 0:
    stopnja = np.linspace(1,100,100)
    for l in stopnja:
        print(l,S(data,l),sep='\t')

start_time = perf_counter()
### Izracuna varianco preko Legendrovih polinomov 
if 1:
    th = th.reshape(1,32)
    vse_niti = []
    for i in range(len(th[0])):
        for n in range(1):
            nit = (Thread(target=varianca, args=(data,th[n][i])))
            nit.start()
            vse_niti.append(nit)
        for nit in vse_niti:
            nit.join()
    end_time = perf_counter()
    print(f'Vzame {end_time- start_time: 0.2f} sekund.')


if 0:
    for theta in th:
        napaka = varianca(data, theta)
    end_time = perf_counter()
    print(f'Vzame {end_time- start_time: 0.2f} sekund.')


### Izracuna variance za veliko stevilo n
if 0:
    th = np.linspace(0,pi/2,101)
    for theta in th:
        n = 10000
        m = stevilo(n,theta,data)
        varianca = np.sum(np.square(m))/n - (np.sum(m)/n)**2
        print(theta/pi,varianca,sep='\t')
