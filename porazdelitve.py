import numpy as np
import math
import random

random.seed()
arange, pi, sin, cos, arccos = np.arange, np.pi, np.sin, np.cos, np.arccos

def nakljucna(N):
    u = np.array([random.uniform(0,1) for i in range(N)])
    v = np.array([random.uniform(0,1) for i in range(N)])

    fi = 2*pi*u
    theta = arccos(2*v-1)

    x,y,z = cos(fi)*sin(theta), sin(fi)*sin(theta), cos(theta)
    matrika = np.stack((x,y,z),axis=-1)

    return np.array(matrika)

def fibonacci(N):
    zlato_razmerje = (1+5**0.5)/2

    i = arange(0,N)
    theta = 2*pi*i/zlato_razmerje
    phi = arccos(1-2*(i+0.5)/N)

    x,y,z = cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)
    matrika = np.stack((x,y,z),axis=-1)

    return np.array(matrika)
