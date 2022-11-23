import math
import random
import copy
import numpy as np
from matplotlib import patches
from scipy.spatial.distance import cdist 

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

        ## generate distribution
        while len(samples) < self.N:
            k = len(samples)//2 + 1

            ## randomly generate k candidates
            candidates = np.array( [
                                    [random.uniform(0,1)*self.width, 
                                    random.uniform(0,1)*self.height] for i in range(k)
                                    ])

            ## calculate distance matrix between candidates and existing points
            distances = cdist(candidates, samples, self.periodic_metric)

            ## select the candidate with the largest minimum distance
            r_min = np.min(distances, axis=1)
            i = np.argmax(r_min)

            generated = candidates[i]
            samples = np.vstack((samples, candidates[i]))

        self.samples = samples
    
    ## define pbc metric
    def check_period(point1, point2, width, height):
        dx = point1[0] - point2[0]
        dy = point1[1] - point2[1]

        if dx > 0.5*width:      dx = dx - width
        elif dx < -0.5*width:   dx = dx + width  
        if dy > 0.5*height:     dy = dy - height
        elif dy < -0.5*height:  dy = dy + height
        return np.array([dx, dy])

    def periodic_metric(self, point1, point2, vector=False):
        dx,dy = distribution.check_period(point1, point2, self.width, self.height)
        return math.sqrt( (dx)**2 + (dy)**2)


class ellipse():

    def __init__(self, center, angle):
        self.center = center
        self.angle = angle

    ## see matplotlib.patches.Ellipse
    def convert_to_patches(self, a, b):
        return patches.Ellipse(self.center, a, b, angle=self.angle*180/math.pi)

class ellipses():
    def __init__(self, distribution, a, b):
        self.grid = distribution.grid
        self.width, self.height = self.grid
        self.N = distribution.N
        self.a, self.b = a, b

        ## create ellipses from centres
        self.ell = np.apply_along_axis(
            lambda center: ellipse(center, random.uniform(0,1)*pi*2), 
            1, distribution.samples   )

    ## input is an ellipse: centre, width(a), height(b) and angle
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
    
    ## mu is maximum value of f
    def mu(self, E1, E2):
        x = np.linspace(0, 1, endpoint=True, num=100)
        y = np.vectorize(self.f)
        Y = y(x, E1, E2)
        return max(Y)
    
    ## energy is contact function between all the neighbours
    def energy(self, j_ellipse, proximity):
        new_E = 0
        if len(proximity) != 0:
            all_mu = np.array( [self.mu(j_ellipse, point_b) for point_b in proximity])
            new_E = np.sum(-np.log(all_mu))
        return new_E
            
    def metropolis(self, proximity, n=100, T=0):
        count, E = 0, np.zeros( shape=(n, self.N) )
        accepted_theta, rejected_theta = np.nan*np.empty(n), np.nan*np.empty(n)

        ## energy of the starting system    
        E[0] = [self.energy(self.ell[i], self.ell[proximity[i]]) for i in range(self.N)]

        for t in range(1, n):
            E[t] = E[t-1]
            ## randomly select an ellipse and copy it
            j = np.random.randint(0, self.N-1)
            j_ellipse = copy.deepcopy(self.ell[j])

            ## generate random angle theta
            j_ellipse.angle = np.random.vonmises(j_ellipse.angle, kappa=1)
            
            ## decide if you want to accept new step or not
            new_E = self.energy(j_ellipse, self.ell[proximity[j]])
            delta_E = new_E - E[t,j]

            ## boltzmann probability distribution
            u = random.uniform(0,1)
            if (delta_E<0 and u<1) or (delta_E>0 and delta_E<=T*math.log(1)/u):
                ##accept this step
                self.ell[j].angle = j_ellipse.angle
                ## change energy
                E[t,j] = new_E
                ## save parametres for later
                accepted_theta[count] = j_ellipse.angle
                count += 1
            else:
                ##save parametres for later
                rejected_theta[t-count] = j_ellipse.angle
                continue
        return [E, accepted_theta, rejected_theta]