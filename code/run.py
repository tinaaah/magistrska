import main
import math
import random
import tikzplotlib
import timeit
import threading
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import patches, cm, colors
from scipy.spatial import distance
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

def draw_state(file_name, Energy, coord_number, variable="energy"):
    ## let's have fun with colours --> map energy values to color palette
    my_palette = LinearSegmentedColormap.from_list('custom red', ["#577590", "#f94144", "#f8961e"], N=256)
    # my_palette = LinearSegmentedColormap.from_list('custom red', ["#ffffff", "#ff0a54"], N=256)

    if variable == "energy": 
        ## heatmap of energy
        # minima = np.round(np.min(Energy), 2)
        minima = 0
        # maxima = np.round(np.max(Energy), 2)
        maxima = 2.6
    
    elif variable == "coord_number":
        ## heatmap of coordination number
        minima = np.round(np.min(coord_number), 2)
        maxima = np.round(np.max(coord_number), 2)
    
    norm = colors.Normalize(vmin=minima, vmax=maxima, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=my_palette)

    plt.figure(figsize=(8,8))

    ## dimensions for patches
    a_width, b_height = 2*new_dist.a, 2*new_dist.b

    for center in new_dist.ell:
        index = np.where(new_dist.ell == center)
        
        if variable == "energy":
            v = Energy[-1,index]

        elif variable == "coord_number":
            v = coord_number[-1,index]

        ellipse = main.ellipse.convert_to_patches(center, a_width, b_height)
        plt.gca().add_patch(ellipse)
        ellipse.set_facecolor(mapper.to_rgba(v)[0])
        ellipse.set_label('_nolegend_')

    cbar = np.round( np.linspace(minima, maxima, 5, endpoint=True), 2)
    plt.colorbar(mapper, ticks=cbar, shrink=0.85, orientation='vertical')

    plt.axis('scaled')
    plt.xlim(0.0, width)
    plt.ylim(0.0, height)
    plt.savefig(file_name, bbox_inches='tight')
    plt.close()

def one_run(a_i, e=5/2, state=False):
    start_time = timeit.default_timer()

    r = np.linspace(1,3,100)
    pic_num = 0

    b_i = a_i/e
    new_dist.a = a_i
    new_dist.b = b_i

    a_width, b_height = 2*a_i, 2*b_i

    ## update a in matrices A
    new_dist.fix_A

    ## distance matrix between all the centres
    S = [x.center for x in new_dist.ell]
    distance_matrix = distance.cdist(S, S, dist.periodic_metric)
    
    ## matrix of indices of neighbouring ellipses, which are less than 2*a apart
    in_proximity = [np.where((line<=2*new_dist.a) & (line!=0)) for line in distance_matrix]

    ## find new state
    Energy, coord_number, accepted, rejected = new_dist.metropolis(in_proximity, n=10)
    correlation = new_dist.angle_correlation(r*new_dist.a)

    if state == True:
        # output_name = "../figures/gifs/delta_a_05/pic" + f"{pic_num:02d}.png"
        output_name = "test.png"
        draw_state(output_name, Energy, coord_number)
        pic_num += 1

    elapsed = timeit.default_timer() - start_time

    return [Energy, coord_number, accepted, rejected, correlation], elapsed

if __name__ == "__main__":
    ## define parametres
    grid = [100, 100]
    width, height = grid
    m, N = 2, 128                   #  number of starting points and number of all points

    ## Selects two random points on a grid
    initial = np.array( [[random.uniform(0,1)*width, random.uniform(0,1)*height] for i in range(m)] )
    ## and and generate distribution with Mitchell algorithm
    dist = main.distribution(initial, grid, N)

    ## change circles into ellipses
    a0, b0 = 5, 2                       # big and small semi-axis
    eps = math.sqrt(1 - (b0/a0)**2)     # eccentricity
    a_width, b_height = 2*a0, 2*b0      # values for matplotlib.pathces

    ## create ellipses
    new_dist = main.ellipses(dist, a0, b0)

    ## distance matrix between all the centres
    S = [x.center for x in new_dist.ell]
    distance_matrix = distance.cdist(S, S, dist.periodic_metric)

    ## run
    parmeters, time = one_run(a_i=7, e=5/2, state=True)
    print(time)

    # threads = []

    # for i in np.arange(5,8,step=0.1):
    #     threads.append(threading.Thread(target=one_run, args=(i, )))
    
    # for thread in threads:
    #     thread.start()
    
    # for thread in threads:
    #     thread.join()
    