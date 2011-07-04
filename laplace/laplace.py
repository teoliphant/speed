from numpy import zeros, abs
from scipy import weave
import numexpr as ne
import pyximport
import numpy as np

dx = 0.1
dy = 0.1
dx2 = dx*dx
dy2 = dy*dy

def dummy_update(u):
    pass

def py_update(u):
    nx, ny = u.shape
    for i in xrange(1,nx-1):
        for j in xrange(1, ny-1):
            u[i,j] = ((u[i+1, j] + u[i-1, j]) * dy2 +
                      (u[i, j+1] + u[i, j-1]) * dx2) / (2*(dx2+dy2))

def num_update(u):
    u[1:-1,1:-1] = ((u[2:,1:-1]+u[:-2,1:-1])*dy2 + 
                    (u[1:-1,2:] + u[1:-1,:-2])*dx2) / (2*(dx2+dy2))

def expr_update(u):
    bottom = u[2:,1:-1]
    top = u[:-2,1:-1]
    left = u[1:-1,2:]
    right = u[1:-1,:-2]
    u[1:-1,1:-1] = ne.evaluate("((bottom + top)*dy2 + "\
                    "(left + right)*dx2) / (2*(dx2+dy2))")
    
def weave_update(u):
    code = """
    int i, j;
    for (i=1; i<Nu[0]-1; i++) {
       for (j=1; j<Nu[1]-1; j++) {
           U2(i,j) = ((U2(i+1, j) + U2(i-1, j))*dy2 + \
                       (U2(i, j+1) + U2(i, j-1))*dx2) / (2*(dx2+dy2));
       }
    }
    """
    weave.inline(code, ['u', 'dx2', 'dy2'])

from laplace import cy_update

import pyximport
import numpy as np
pyximport.install(setup_args={'include_dirs':[np.get_include()]})
from _laplace import cy_update as cy_update2
from _laplace_for import for_update1, for_update2

def calc(N, Niter=100, func=py_update, args=()):
    order = 'C' if func not in [for_update1, for_update2] else 'F'
    u = zeros([N, N], order=order)
    u[0] = 1
    for i in range(Niter):
        func(u,*args)
    return u
    
