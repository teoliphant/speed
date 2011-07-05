from numpy import zeros, abs
from scipy import weave
import numexpr as ne
import pyximport
import numpy as np
from _laplace_for import for_update1, for_update2
pyximport.install(setup_args={'include_dirs':[np.get_include()]})
from _laplace import cy_update 


dx = 0.1
dy = 0.1
dx2 = dx*dx
dy2 = dy*dy


def num_update(u, dx2, dy2):
    u[1:-1,1:-1] = ((u[2:,1:-1]+u[:-2,1:-1])*dy2 + 
                    (u[1:-1,2:] + u[1:-1,:-2])*dx2) / (2*(dx2+dy2))

def expr_update(u, dx2, dy2):
    bottom = u[2:,1:-1]
    top = u[:-2,1:-1]
    left = u[1:-1,2:]
    right = u[1:-1,:-2]
    u[1:-1,1:-1] = ne.evaluate("((bottom + top)*dy2 + "\
                    "(left + right)*dx2) / (2*(dx2+dy2))")
    
def weave_update(u, dx, dy2):
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


def calc(N, Niter=100, func=num_update, args=()):
    order = 'C' if func not in [for_update1, for_update2] else 'F'
    u = zeros([N, N], order=order)
    u[0] = 1
    for i in range(Niter):
        func(u,*args)
    return u



def main():
    import time
    N = 150
    Niter = 8000
    modes = [['NumPy', num_update, (dx2, dy2)],
             ['Numexpr', expr_update, (dx2, dy2)],
             ['Cython', cy_update, (dx2, dy2)],
             ['Weave', weave_update, (dx2, dy2)],
             ['Looped Fortran', for_update1, (dx2, dy2, N, N)],
             ['Vectorized Fortran', for_update2, (dx2, dy2, N, N)]
             ]
    for mode, update, args in modes:
        start = time.time()
        calc(N, Niter, func=update, args=args)
        elapsed = time.time() - start
        print "%s: %f seconds" % (mode, elapsed)

if __name__ == '__main__':
    main()
    
