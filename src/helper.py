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
    This is the Rosenbrock Function
    Bound: X1=[-5,5], X2=[-2,8]
    Global Optimum: 0, at (1,1)
    '''
    x1 = x[0]
    x2 = x[1]
    a  = 100.0
    f  = a * (x2 - x1**2)**2 + (1 - x1)**2
    return f


def func3(x):
    '''
    This is the Six-hump Camelback Function.
    Bound: X1=[-3,2], X2=[-3,2]
    True Optima: -1.031628453489877, at (-0.08983,0.7126), (0.08983,-0.7126)
    '''
    x1 = x[0]
    x2 = x[1]
    f  = (4 - 2.1*x1**2 + x1**4/3)*x1**2 + x1*x2 + (-4 + 4*x2**2)*x2**2
    return f


def func4(x):
    '''
    This is the Rastrigin Function
    Bound: X1=[-1,1], X2=[-1,1]
    Global Optimum: -2, at origin
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

def func6(x):
    '''
    This is the Hartmann Function (6-D)
    Bound: X(i)=[0,1], for i=1,2,...,6
    Global Optimum: -3.32237, at (0.20169, 0.150011, 0.476874, 0.275332, 0.311652, 0.6573)
    '''
    alpha = np.transpose(np.array([1., 1.2, 3., 3.2]))
    A = np.array([[10, 3, 17, 3.5, 1.7, 8],
                 [0.05, 10, 17, 0.1, 8, 14],
                 [3, 3.5, 1.7, 10, 17, 8],
                 [17, 8, 0.05, 10, 0.1, 14]])
    P = 1e-4 * np.array([[1312, 1696, 5569, 124, 8283, 5886],
                        [2329, 4135, 8307, 3736, 1004, 9991],
                        [2348, 1451, 3522, 2883, 3047, 6650],
                        [4047, 8828, 8732, 5743, 1091, 381]])
    u1 = 0
    for i in range(4):
        u2 = 0
        for j in range(6):
            u2 += A[i,j] * (x[j]-P[i,j])**2
        u1 += alpha[i]**(-u2)
    f = -(2.58 + u1) / 1.94
    return f