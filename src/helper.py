import numpy as np

def func1(x):
    '''
    This is the Goldstein-Price Function
    Bound X1=[-2,2], X2=[-2,2]
    Global Optimum: 3.0, at (0,-1)
    '''
    x1 = x[0]
    x2 = x[1]
    u1 = (x1 + x2 + 1.0)**2
    u2 = 19. - 14.*x1 + 3.*x1**2 - 14.*x2 + 6.*x1*x2 +3.*x2**2
    u3 = (2.*x1 - 3.*x2)**2
    u4 = 18. - 32.*x1 + 12.*x1**2 + 48.*x2 -36.*x1*x2 + 27.*x2**2
    u5 = u1 * u2
    u6 = u3 * u4
    f  = (1. + u5) * (30. + u6)
    return f


def func2(x):
    '''
    %  This is the Rosenbrock Function
    %  Bound: X1=[-5,5], X2=[-2,8]
    %  Global Optimum: 0, at (1,1)
    '''
    x1 = x[0]
    x2 = x[1]
    a  = 100.0
    f  = a * (x2 - x1**2)**2 + (1 - x1)**2
    return f


def func3(x):
    '''
    %  This is the Six-hump Camelback Function.
    %  Bound: X1=[-3,2], X2=[-3,2]
    %  True Optima: -1.031628453489877, at (-0.08983,0.7126), (0.08983,-0.7126)
    '''
    x1 = x[0]
    x2 = x[1]
    f  = (4 - 2.1*x1**2 + x1**4/3)*x1**2 + x1*x2 + (-4 + 4*x2**2)*x2**2
    return f


def func4(x):
    '''
    %  This is the Rastrigin Function
    %  Bound: X1=[-1,1], X2=[-1,1]
    %  Global Optimum: -2, at origin
    '''
    x1 = x[0]
    x2 = x[1]
    f  = x1**2 + x2**2 - np.cos(18.0*x1) - np.cos(18.0*x2)
    return f


def func5(x, nopt=2):
    '''
    This is the Griewank Function (2-D or 10-D)
    Bound: X(i)=[-600,600], for i=1,2,...,10
    Global Optimum: 0, at origin
    '''
    if nopt==2:    d = 200.0
    else:          d = 4000.0
    u1 = 0.0
    u2 = 1.0
    for j in range(nopt):
        u1 = u1 + x[j]**2/d
        u2 = u2 * np.cos(x[j]/np.sqrt(float(j+1)))
    f = u1 - u2 + 1
    return f