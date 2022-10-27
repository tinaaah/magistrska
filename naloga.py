import math
import random
import numpy as np
from matplotlib import patches
from scipy.spatial import distance

cos = math.cos
sin = math.sin
pi = math.pi

## A = R(theta) * A' * R^T(theta)
def Rotate(theta, matrix):
    R = np.matrix([ [cos(theta), -sin(theta)],
                    [sin(theta), cos(theta)]    ])
    R_t = R.transpose()

    return R@matrix@R_t


class distribution():
    def __init__(self, samples, grid, N):
        self.grid = grid
        self.width, self.height = grid
        self.N = N

        ## generira porazdelitev
        while len(samples) < self.N:
            k = len(samples)//2 + 1

            ## randomly generate k candidates 
            candidates = np.array( [
                                    [random.uniform(0,1)*self.width, 
                                    random.uniform(0,1)*self.height] for i in range(k)
                                    ])

            ## calculate distance matrix between candidates and existing points
            distances = distance.cdist(candidates, samples, self.periodic_metric)

            ## select the candidate with the largest minimum distance
            r_min = np.min(distances, axis=1)
            i = np.argmax(r_min)

            generated = candidates[i]
            samples = np.vstack((samples, candidates[i]))

        self.samples = samples
    
    ### define pbc metric
    def check_period(point1, point2, width, height):
        dx = point1[0] - point2[0]
        dy = point1[1] - point2[1]

        if dx > 0.5*width:
            dx = dx - width
        elif dx < -0.5*width:
            dx = dx + width  
        if dy > 0.5*height:
            dy = dy - height
        elif dy < -0.5*height:
            dy = dy + height
        return np.array([dx, dy])

    def periodic_metric(self, point1, point2, vector=False):
        dx,dy = distribution.check_period(point1, point2, self.width, self.height)
        return math.sqrt( (dx)**2 + (dy)**2)


class ellipse():

    def __init__(self, center, angle):
        self.center = center
        self.angle = angle

    def convert_to_patches(self, a, b):
        return patches.Ellipse(self.center, a, b, angle=self.angle*180/math.pi)

class ellipses():
    def __init__(self, distribution, a, b):
        self.grid = distribution.grid
        self.width, self.height = self.grid
        self.N = distribution.N
        self.a, self.b = a, b

        ## v centrih naredim elipse
        self.ell = np.apply_along_axis(
            lambda center: ellipse(center, random.uniform(0,1)*pi*2), 
            1, distribution.samples   )
    
    ##  sprejme ellipse: center, širino(a), višino(b) in kot
    ## f = lambda*(1-lambda)*r^T * A^(-1) * r
    def f(self, x, E1, E2):
        A1 = np.diag( [self.a**2, self.b**2] )
        theta1 = E1.angle
        A1 = Rotate(theta1, A1)

        A2 = np.diag( [self.a**2, self.b**2] )
        theta2 = E2.angle
        A2 = Rotate(theta2, A2)

        dr = distribution.check_period(E1.center, E2.center, self.width, self.height)
        C = np.linalg.inv( (1-x)*A1 + x*A2 )
        return (x*(1-x)*dr@C@dr.T)[0,0]
    
    def mu(self, E1, E2):
        x = np.linspace(0, 1, endpoint=True, num=100)
        y = np.vectorize(self.f)
        Y = y(x, E1, E2)
        return max(Y)