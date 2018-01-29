# -*- coding: utf-8 -*-
"""
eifs.py

Creates file eifs.stan containing an array whose rows are the 
first p AVE approximations of some subspace of an RKHS.

The argument 'rkhs' specifies a basis of the subspace.
'h1' uses the eigenfunctions of the first order Green's kernel.
'pwc' uses piecewise constants
'rbf' uses non-normalized Gaussians
'linear' uses linear splines
"""

def compute_eifs(p,n,tau,rkhs='h1',rbfwidth=1):
    import numpy as np
    from scipy import integrate
    eifs = np.zeros((p,n),dtype=float)
    T = n*tau

    if rkhs == 'h1':
        for j in range(p):
            e = lambda s: np.sqrt(2. / T) * np.sin((s / T) * (j * np.pi - np.pi / 2))
            for i in range(1,n):
                eifs[j,i] = tau**(-1)*integrate.quad(e,tau*(i-1),tau*i)[0]
    if rkhs == 'pwc':
        assert(p<n and n%p==0)
        k = int(n/p)
        for i in range(p):
            eifs[i,(i*k):((i+1)*k)] = np.repeat(np.array([1/n]),k)
    if rkhs == 'rbf':
        # Not normalized
        for j in range(p):
            e = lambda s: (2*np.pi*rbfwidth)**(-2)*np.exp(-(0.5*rbfwidth**(-2))*(s-T*j/p)**2)
            for i in range(1,n):
                eifs[j,i] = tau**(-1)*integrate.quad(e,tau*(i-1),tau*i)[0]
    if rkhs == 'linear':
        w = (n-1)*tau / p
        for j in range(p):
            e = lambda x: (j*w < x <= (j+1)*w)*(x-j*w)/w + ((j+1)*w< x <= (j+2)*w)*(-(x-(j+2)*w))/w
            for i in range(1,n):
                eifs[j,i] = tau**(-1)*integrate.quad(e,tau*(i-1),tau*i)[0]
    return np.array(eifs)

def save_eifs(eifs,filename):
    import numpy
    numpy.set_printoptions(threshold = 1e6)
    with open(filename,'w') as file:
        file.write('eifs = ' + numpy.array2string(eifs,separator=',') + ';')
    
def eifs(p,n,tau,rkhs='h1',rbfwidth=1):
    save_eifs(compute_eifs(p,n,tau,rkhs),'../eifs.stan')
