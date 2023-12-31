from stl import mesh
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from inspect import signature

#the purpose of this folder is to ensurre consistency of numerical schemes used for simulating our system,
#the probability of our simulatino and other design features

#sample of how to write ur sim is below
def euler(h : float, a : int):
    result = []
    guy = 1  #-- this is the initial conditon
    t = np.linspace(0,1,int(1/h))
    result.append(guy)
    for i in t:
        guy = guy + h*(np.exp(i))
        result.append(guy)
    return result

def check_consistency(func, h0, order):
    N = np.array([h0,h0/10,h0/100,h0/1000, h0/100000])
    M = 2*N
    p = []
    tester = True
    for i in range(5):
        eul = func(N[i],0)
        eul2 =  func(M[i],0)
        p.append( abs(eul[-1] - ((2**order)*eul2[-1] - eul[-1])/(2**order - 1)))
    for i in range(4):
        print(abs((np.log(p[i+1]) - np.log(p[i]))/(np.log(N[i]) - np.log(N[i+1]))))
        tester = tester and abs(abs((np.log(p[i+1]) - np.log(p[i]))/(np.log(N[i]) - np.log(N[i+1]))) - order) < 0.1
    plt.loglog(p,N)
    plt.show()
    if not tester:
        print("failed, check implementation of scheme")
        return
    print("passed")
    return

check_consistency(euler, 0.05, 1)

def check_stability(func):
    pass


def check_probability(func):
    num_inputs = signature(func)
    print(type(num_inputs.parameters['h'].annotation))
    
    return

def expected_performance(func, ranges):
    for i in range(len(ranges)):
        




    

