from direct.direct import Direct
from problem.helper import *

def result(d):
    curr_opt, x_at_opt, l_hist = d.run()
    print("curr_opt =", curr_opt, ",  x_at_opt =", x_at_opt)
    print()


if __name__ == "__main__":
    print('test Goldstein-Price results:')
#    x=np.array([2,2])
    d = Direct(func1, bounds=np.array([[-2,2,],[-2,2,]]))
    result(d)
  
    print('test Rosenbrock results:')
#    x = np.array([-2.,7.])
    d = Direct(func2, bounds=np.array([[-5,5],[-2,8]]))
    result(d)
   
    print('test Six-hump Camelback results:')
#    x = np.array([-1.,1.])
    d = Direct(func3, bounds=np.array([[-5,5],[-5,5]]))
    result(d)
     
    print('test Rastrigin results:')
#    x = np.array([-1.,1.])
    d = Direct(func4, bounds=np.array([[-1,1],[-1,1]]))
    result(d)
 
    print('test Griewank results:')
#    x = np.zeros(2)
    d = Direct(func5, bounds=np.array([[-600,600],[-600,600]]))
    result(d)