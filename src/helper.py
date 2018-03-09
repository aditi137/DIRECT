import numpy as np

def func1(x):
    '''
    This is the Goldstein-Price function
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


def func2(x, a=100.0):
    '''
    This is the Rosenbrock function
    Bound: X1=[-5,5], X2=[-2,8]
    Global Optimum: 0, at (1,1)
    '''
    x1 = x[0]
    x2 = x[1]
    f  = a * (x2 - x1**2)**2 + (1 - x1)**2
    return f


def func3(x):
    '''
    This is the Six-hump Camelback function.
    Bound: X1=[-3,2], X2=[-3,2]
    True Optima: -1.031628453489877, at (-0.08983,0.7126), (0.08983,-0.7126)
    '''
    x1 = x[0]
    x2 = x[1]
    f  = (4 - 2.1*x1**2 + x1**4/3)*x1**2 + x1*x2 + (-4 + 4*x2**2)*x2**2
    return f


def func4(x):
    '''
    This is the Rastrigin function
    Bound: X1=[-1,1], X2=[-1,1]
    Global Optimum: -2, at origin
    '''
    x1 = x[0]
    x2 = x[1]
    f  = x1**2 + x2**2 - np.cos(18.0*x1) - np.cos(18.0*x2)
    return f


def func5(x, nopt=10):
    '''
    This is the Griewank function (2-D or 10-D)
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


def func6(x, nopt=3):
    '''
    This is the Hartmann function (3-D or 6-D)
    Bound: X(i)=[0,1], for i=1,2,...,6
    Global Optimum (6-D): -3.32237, at (0.20169, 0.150011, 0.476874, 0.275332, 0.311652, 0.6573)
    Global Optimum (3-D): -3.86278, at (0.114614, 0.555649, 0.852547)
    '''
    alpha = np.array([1., 1.2, 3., 3.2])
    if nopt==6:
        A = np.array([[10, 3, 17, 3.5, 1.7, 8],
                      [0.05, 10, 17, 0.1, 8, 14],
                      [3, 3.5, 1.7, 10, 17, 8],
                      [17, 8, 0.05, 10, 0.1, 14]])
        P = 1e-4 * np.array([[1312, 1696, 5569, 124, 8283, 5886],
                             [2329, 4135, 8307, 3736, 1004, 9991],
                             [2348, 1451, 3522, 2883, 3047, 6650],
                             [4047, 8828, 8732, 5743, 1091, 381]])
    else:
        A = np.array([[3.0, 10, 30],
                      [0.1, 10, 35],
                      [3.0, 10, 30],
                      [0.1, 10, 35]])
        P = 1e-4 * np.array([[3689, 1170, 2673],
                             [4699, 4387, 7470],
                             [1091, 8732, 5547],
                             [381, 5743, 8828]])
    u1 = 0
    for i in range(4):
        u2 = 0
        for j in range(nopt):
            u2 += A[i,j] * (x[j]-P[i,j])**2
        u1 += alpha[i]**(-u2)
    if nopt==6:
        f = -(2.58 + u1) / 1.94
    else:
        f = -u1
    return f


def func7(x, a=1, b=5.1/(4*np.pi**2), c=5/np.pi, r=6, s=10, t=1/(8*np.pi)):
    '''
    This is the Branin function
    Bound: X1=[-5,10], X2=[0,15]
    Global Optimum: 0.397887, at (-pi,12.275), (pi,2.275), (9.42478,2.475)
    '''
    x1 = x[0]
    x2 = x[1]
    f = a * (x2 - b*x1**2 + c*x1 - r)**2 + s*(1-t)*np.cos(x1) + s
    return f


def func8(x, m=5):
    '''
    This is the Shekel function
    Bound: X(i)=[0,10], for i=1,...,4
    Global Optimum (m=10): -10.5364, at (4,4,4,4)
    Global Optimum (m=7): -10.4029, at (4,4,4,4)
    Global Optimum (m=5): -10.1532, at (4,4,4,4)
    '''
    b = 0.1 * np.array([1, 2, 2, 4, 4, 6, 3, 7, 5, 5])
    C = np.array([[4., 1., 8., 6., 3., 2., 5., 8., 6., 7.],
                  [4., 1., 8., 6., 7., 9., 3., 1., 2., 3.6],
                  [4., 1., 8., 6., 3., 2., 5., 8., 6., 7.],
                  [4., 1., 8., 6., 7., 9., 3., 1., 2., 3.6]])
    u1 = 0
    for i in range(m):
        u2 = 0
        for j in range(4):
            u2 += (x[j] - C[j, i])**2
        u1 += 1/(u2 + b[i])
    f = -u1
    return f


def func9(x):
    '''
    This is the Shubert function
    Bound: X(i)=[-10,10], for i=1,...,4
    Global Optimum: -186.7309
    '''
    u1 = 0. 
    u2 = 0.
    x1 = x[0]
    x2 = x[1]
    for i in range(1,6):
        u1 += i * np.cos((i+1) * x1 + i)
        u2 += i * np.cos((i+1) * x2 + i)
    f = u1 * u2
    return f