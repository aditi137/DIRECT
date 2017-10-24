from direct import Direct, GlobalMin
from helper import *

def result(f, bounds, **kwargs):
    curr_opt, x_at_opt, l_hist = Direct(f,bounds, **kwargs).run()
    print("curr_opt =", curr_opt, ",  x_at_opt =", x_at_opt)
    print()


if __name__ == "__main__":
    print('test Goldstein-Price results:')
    result(func1, bounds=np.array([[-2,2,],[-2,2,]]), max_feval=10000, max_iter=10000)
  
    print('test Rosenbrock results:')
    result(func2, bounds=np.array([[-5,5],[-2,8]]), max_feval=10000, max_iter=10000)
   
    print('test Six-hump Camelback results:')
    result(func3, bounds=np.array([[-5,5],[-5,5]]))
      
    print('test Rastrigin results:')
    result(func4, bounds=np.array([[-1,1],[-1,1]]))
  
    print('test Griewank results:')
    result(func5, bounds=np.array([[-600,600],[-600,600]]))