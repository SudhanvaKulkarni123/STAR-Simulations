import numpy as np
import matplotlib.pyplot as plt
from sklearn import neighbors


DOMAIN_WIDTH = 500
DOMAIN_HEIGHT = 800
BASE_DENSITY = 1
FIGURE_SIZE = (4, 6)
PLOT_EVERY = 6
SCATTER_DOT_SIZE = 200

class Kernels:
    def general_kernel(h,r):
        return 315.0*((h**2 - r**2)**3)/(64*np.pi*(h**9))

    def grad_gen(h,r):
        return (945.0*((h**2 - r**2)**2)/(64*np.pi*(h**6)))
    
    def lapl_gen(h,r):
        return (315.0*((h**2 - r**2))/(64*np.pi*(h**3)))

    def __pressure_kernel(h,r):
        return 15.0*((h - r)**3)/(np.pi*(h**6))
    
    def __viscosity_kernel(h,r):
        return -2*(r/h)**3 + (r/h)**2 + h/(2*r) - 1

def main(T,N, total, smoothing, iso, base, rho, bounce, viscosity):    #T total time, N number of teps per second

    x = np.zeros((total, 2))  #position
    x = np.zeros((total, 2))   #positions at previous step for verlet integration
    v = np.zeros_like(x)   #velocity
    f = np.zeros_like(x)   #force
    m = np.zeros(total)      #mass
    p = np.zeros(total)       #particles
    step = 1/N                      #for verlet
    astep = step**2

    DOMAIN_X_LIM = np.array([smoothing, DOMAIN_WIDTH - smoothing,])
    DOMAIN_Y_LIM = np.array([smoothing,DOMAIN_HEIGHT - smoothing,])

    for i in range(total):
        new_positions = np.array([
                [np.random.uniform(0.0,50.0), DOMAIN_Y_LIM[1]],
            ])
        new_velocities = np.array([
                [0.0, np.random.normal(10.0,2.0)] 
            ])
        new_mass = np.array([np.random.normal(1.0,0.5)])
        x = np.concatenate((x, new_positions), axis=0)
        v = np.concatenate((v, new_velocities), axis=0)
        m = np.concatenate((m,new_mass), axis= 0)

    for i in range(int(T*N)):
        
        tree = neighbors.KDTree(x,8)
        ids, dist = tree.query_radius(x,smoothing,return_distance=True,sort_results=True,)
        rho = np.zeros(total)     #density
        for k in range(total):
            for j_in_list, j in enumerate(ids[k]):
                rho[k] += m[k]*Kernels.general_kernel(smoothing, dist[k][j_in_list]) 
            p = (rho - base)*iso
        f = np.zeros(total)     #forces
        for k in range(total):
            for j_in_list, j in enumerate(ids[k]):
                f[k] += m[k]*p[k]*Kernels.grad_gen(smoothing, dist[k][j_in_list])
                f[k] += m[k]*v[k]*Kernels.lapl_gen(smoothing, dist[k][j_in_list])*viscosity
        if i == 0:
            x = x + (step)*v + 0.5*(astep)*f/rho[:, np.newaxis]
            
        else :
            tmp = x
            x = 2*x - prev_x + (astep)*f/rho[:, np.newaxis]
            prev_x = tmp
        v = v + (step)*f/rho[:, np.newaxis]
        out_of_left_boundary = x[:, 0] < DOMAIN_X_LIM[0]
        out_of_right_boundary = x[:, 0] > DOMAIN_X_LIM[1]
        out_of_bottom_boundary = x[:, 1] < DOMAIN_Y_LIM[0]
        out_of_top_boundary = x[:, 1] > DOMAIN_Y_LIM[1]

        v[out_of_left_boundary, 0]     *= bounce
        x [out_of_left_boundary, 0]      = DOMAIN_X_LIM[0]

        v[out_of_right_boundary, 0]    *= bounce
        x[out_of_right_boundary, 0]     = DOMAIN_X_LIM[1]

        v[out_of_bottom_boundary, 1]   *= bounce
        x[out_of_bottom_boundary, 1]    = DOMAIN_Y_LIM[0]

        v[out_of_top_boundary, 1]      *= bounce
        x[out_of_top_boundary, 1]       = DOMAIN_Y_LIM[1]


        plt.style.use("dark_background")
        plt.figure(figsize=FIGURE_SIZE, dpi=160)
        if iter % PLOT_EVERY == 0:
            plt.scatter(
                x[:, 0],
                x[:, 1],
                s=SCATTER_DOT_SIZE,
                c=x[:, 1],
                cmap="Wistia_r",
            )
            plt.ylim(DOMAIN_Y_LIM)
            plt.xlim(DOMAIN_X_LIM)
            plt.xticks([], [])
            plt.yticks([], [])
            plt.tight_layout()
            plt.draw()
            plt.pause(0.0001)
            plt.clf()


if __name__ == "__main__":
    time = 10
    N = 10
    total = 50
    smoothing_length = 9
    bounce = -0.9
    base = 1
    iso = 20
    rho = 1
    viscosity = 0.001
    main(time,N, total=total, smoothing=smoothing_length, iso=iso, base=base, rho = 1, bounce=bounce,viscosity=viscosity)