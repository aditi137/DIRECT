import numpy as np
import pylab as plt

from direct import Direct

def prog1():
    f = lambda x: (x**2)[0]
    b = np.array([[0,1]])
    d = Direct(f, b, max_iter=4)
    plt.figure()
    x = np.linspace(b[0,0], b[0,1], 100).reshape((1,-1))
    y = f(x)
    plt.plot(x.flatten(),y, label='true function')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(loc='best')
    i = 0
    for hist in d.l_hist:
        plt.plot(hist[0], hist[1], 'r.')
        plt.text(hist[0], hist[1]+0.05, i)
        i += 1
    d.run()
    plt.show()

def func_1d(x):
    return ((x-0.5)**2 + np.sin(x*2.1)* np.sin(x*12)*1.7)
    
def prog2():
    b = np.array([[-1,2]])
    d = Direct(func_1d, b, max_feval=12)
    d.run()
    plt.figure()
    x = np.linspace(b[0,0], b[0,1], 100)
    y = np.array([func_1d(xx) for xx in x])
    plt.plot(x.flatten(),y, label='true function')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(loc='best')
    i = 0
    for hist in d.l_hist:
        plt.plot(hist[0], hist[1], 'r.')
        plt.text(hist[0], hist[1]+0.05, i)
        i += 1
    plt.show()

def prog3():
    b = np.array([[-1,2]])
    d = Direct(func_1d, b, max_feval=12)
    d.run()
    plt.figure()
    x = np.linspace(b[0,0], b[0,1], 100)
    y = np.array([func_1d(xx) for xx in x])
    h = np.array([yy for xx, yy in d.l_hist])
    search = np.minimum.accumulate(h)
    nRand = 100
    rand_search = np.zeros((search.shape[0],))
    for i in range(nRand):
        np.random.shuffle(y)
        rand_search += np.minimum.accumulate(y[:search.shape[0]])
    rand_search /= float(nRand)
    plt.plot(search, label='DiRect')
    plt.plot(rand_search, label='random search avg. 100 times')
    plt.legend(loc='best')
    plt.xlabel('iteration')
    plt.ylabel('found optimum')
    plt.show()
    
if __name__ == "__main__":
    prog1()
    prog2()
    prog3()