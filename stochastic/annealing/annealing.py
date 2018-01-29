## GLOBAL OPTIMIZATION


import numpy as np
from Parabolic_System import Parabolic_System
from scipy.stats import beta
import matplotlib.pyplot as plt
from scipy.optimize import basinhopping, differential_evolution
from eifs import compute_eifs


np.random.seed(323)

# Define system parameters

t = np.linspace(0,3,50)
n = t.shape[0]
N = 16
p = n//3

q1 = 1.1
q2 = 9.75

ds = 1 # downsampling rate
tau = (t[1]-t[0])*ds

z = Parabolic_System(q1, q2, n, t, N)

M=compute_eifs(p,n,tau,'linear').T

# Define artificial input signal

f = lambda s: 0.3*beta.pdf(s,12,7) + 0.6*beta.pdf(s,4,11) if s<=1 else 0

u = np.abs(np.array([f(s) for s in t]) + np.random.normal(scale=0.02,size=(50,)))
y = z.fwd(u)

def J_dJ(qcreg):
    from Parabolic_System import Parabolic_System
    from numpy.linalg import norm
   
    q1,q2 = qcreg[:2]
    c = qcreg[2:(2+p)]
    reg = qcreg[1:]

    z = Parabolic_System(q1,q2,n,t,N)
    model = z.fwd(M @ c)
    
    
    
    return norm(model-y,2) + reg*norm(M @ c,2)

initial = np.array([1,1]+[0]*(p)+[0.1])
minimizer_kwargs = {"method": "L-BFGS-B", "jac":True}

fit = basinhopping(J_dJ,initial, minimizer_kwargs = minimizer_kwargs,
                   niter = 10, callback=print_fun)
