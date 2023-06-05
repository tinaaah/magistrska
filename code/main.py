import math
import random
import copy
import numpy as np
from matplotlib import patches
from scipy.spatial.distance import cdist
from scipy.optimize import brute

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
        self.accepted_rot = 0
        self.rejected_rot = 0

    ## see matplotlib.patches.Ellipse
    def convert_to_patches(self, a, b):
        return patches.Ellipse(self.center, a, b, angle=self.angle*180/math.pi)

class ellipses():
    def __init__(self, distribution, a, b):
        self.grid = distribution.grid
        self.periodic_metric = distribution.periodic_metric
        self.width, self.height = self.grid
        self.N = distribution.N
        self.a, self.b = a, b
        self.x = np.linspace(0, 1, endpoint=True, num=100)

        ## create ellipses from centres
        self.ell = np.apply_along_axis(
            lambda center: ellipse(center, random.uniform(0,1)*pi*2), 
            1, distribution.samples   )

        ## define ellipse matrices for all ellipses
        for sample in self.ell:
            sample.A = Rotate(sample.angle, np.diag([self.a**2, self.b**2]))
        
    ## fix ellipse matrices when changing size
    @property
    def fix_A(self):
        for sample in self.ell:
            sample.A = Rotate( sample.angle, np.diag( [self.a**2, self.b**2] ))
    

    ## orientational correlation function
    def psi2(self, ellipse, proximity):
        if len(proximity) < 2:
            psi = 1
        else:
            psi = np.sum([np.exp(2*1j*(ellipse.angle-ellipse_b.angle)) for ellipse_b in proximity])/len(proximity)
        return psi
    
    def angle_correlation(self, radius_range):
            g = np.zeros(len(radius_range))

            for i in range(len(radius_range)):
                d_a = radius_range[i]

                ## distance matrix between all the centres
                S = [x.center for x in self.ell]
                distance = cdist(S, S, self.periodic_metric)

                ## matrix of neighbouring ellipses, which are less than d_a away
                neighbours = [np.where((line<=d_a) & (line!=0)) for line in distance]

                ## calculate all psi_2 at current r
                psi_2 = [self.psi2(self.ell[j], self.ell[neighbours[j]]) for j in range(self.N)] 

                ## g(r) for each j
                g[i] = np.average( [[np.abs(psi_2[j]*np.conjugate(psi_2[k])) for k in range(self.N)] for j in range(self.N)])
            return g

    ## input is an ellipse: centre, width(a), height(b) and angle
    ## f = lambda*(1-lambda)*r^T * A^(-1) * r
    def f(self, x, E1, E2):
        A1 = E1.A
        A2 = E2.A

        dr = distribution.check_period(E1.center, E2.center, self.width, self.height)
        C = np.linalg.inv( (1-x)*A1 + x*A2 )
        # C = np.linalg.inv( (1-x[0])*A1 + x[0]*A2 )
        return (x*(1-x)*dr@C@dr.T)[0,0]
    
    ## mu is maximum value of f
    def mu(self, E1, E2):
        x = self.x
        y = np.vectorize(self.f)
        Y = y(x, E1, E2)
        return max(Y)
    
    ## energy is contact function between all the neighbours
    def energy(self, j_ellipse, proximity):
        new_E, new_Z = 0, 0
        if len(proximity) != 0:
            all_mu = np.array( [self.mu(j_ellipse, point_b) for point_b in proximity])

            ## energy is zero unless in contact
            E = -np.log(all_mu[np.where(all_mu<=1)])
            new_Z = len(E)
            new_E = np.sum(E)
        return np.asarray([new_E, new_Z])
            
    def metropolis(self, proximity, n=100, T=0):
        t, E, Z = 0, np.zeros( shape=(n, self.N) ), np.zeros( shape=(n, self.N) )
        accepted_theta, rejected_theta = np.nan*np.ones(n), np.nan*np.ones(n)

        ## energy of the starting system    
        temp = np.asarray([self.energy(self.ell[i], self.ell[proximity[i]]) for i in range(self.N)])
        E[0] = temp[:,0]
        Z[0] = temp[:,1]

        for i in range(1, n):
            E[i] = E[i-1]
            Z[i] = Z[i-1]

            ## randomly select an ellipse and copy it
            j = np.random.randint(0, self.N-1)
            j_ellipse = copy.deepcopy(self.ell[j])
            
            ## generate random angle theta
            j_ellipse.angle = np.random.vonmises(j_ellipse.angle, kappa=3)
            j_ellipse.A = Rotate(j_ellipse.angle, np.diag( [self.a**2, self.b**2] ))

            ## because delta_theta can get to 2*pi
            def transform(x):
                if x > math.pi:     x = x - 2*math.pi
                elif x < -math.pi:  x = x + 2*math.pi
                return x
            delta_theta = transform(self.ell[j].angle - j_ellipse.angle)
           
            ## decide if you want to accept new step or not
            new_E, new_Z = self.energy(j_ellipse, self.ell[proximity[j]])
            delta_E = new_E - E[i,j]

            ## boltzmann probability distribution
            u = random.uniform(0,1)

            if (delta_E<=0 and u<1) or (delta_E>0 and delta_E<=T*math.log(1/u)):
                ## accept this step
                self.ell[j].angle = j_ellipse.angle
                self.ell[j].A = Rotate(j_ellipse.angle, np.diag( [self.a**2, self.b**2] ))

                ## save parameters for later
                ## change energy
                E[i,j] = new_E
                Z[i,j] = new_Z

                ## change energy of all neighbors as well
                neighbourhood = self.ell[proximity[j]]
                for neighbour in neighbourhood:
                    k = np.where(self.ell == neighbour)[0][0]
                    E[i,k], Z[i,k] = self.energy(self.ell[k], self.ell[proximity[k]])

                ## delta_theta
                accepted_theta[t] = delta_theta
                self.ell[j].accepted_rot += 1
                t += 1

            else:
                ## save parameters for later
                rejected_theta[i-t] = delta_theta
                self.ell[j].rejected_rot += 1
                continue

            ## If energy reaches zero stop iterating
            if sum(E[i]) == 0:
                break
        
        ## calculate energy for all pairs after the final rotation
        temp = np.asarray([self.energy(self.ell[i], self.ell[proximity[i]]) for i in range(self.N)])
        E[i] = temp[:,0]
        Z[i] = temp[:,1]
        return [E, Z, accepted_theta, rejected_theta]