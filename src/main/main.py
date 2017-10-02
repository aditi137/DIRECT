import numpy as np

from direct.direct import Direct
from problem.helper import call_obj_func

def callfunc(bound, x, testnr):
    f = call_obj_func(x, testnr)
    d = Direct(f, bound)
    d.run()
    print()

if __name__ == "__main__":
    bound = np.array([[-2,2],[-2,2]])
    x = np.array([2,2])
    callfunc(bound, x, testnr=1)
    
    bound = np.array([[-5,5],[-2,8]])
    x = np.array([-2.,7.])
    callfunc(bound, x, testnr=2)
    
    bound = np.array([[-5,5],[-5,5]])
    x = np.array([-1.,1.])
    callfunc(bound, x, testnr=3)
    
    bound = np.array([[-1,1],[-1,1]])
    x = np.array([-1.,1.])
    callfunc(bound, x, testnr=4)

    bound = np.array([[-600,600],[-600,600]])
    x = np.zeros(2)
    callfunc(bound, x, testnr=5)