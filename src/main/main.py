from direct.direct import Direct
from problem.helper import *

def result(f, bounds):  #TODO: add kwargs for Direct(max_feval,max_iter)
    curr_opt, x_at_opt, l_hist = Direct(f,bounds).run()
    print("curr_opt =", curr_opt, ",  x_at_opt =", x_at_opt)
    print()


if __name__ == "__main__":
    print('test Goldstein-Price results:')
    result(func1, bounds=np.array([[-2,2,],[-2,2,]]))
  
    print('test Rosenbrock results:')
    result(func2, bounds=np.array([[-5,5],[-2,8]]))
   
    print('test Six-hump Camelback results:')
    result(func3, bounds=np.array([[-5,5],[-5,5]]))
      
    print('test Rastrigin results:')
    result(func4, bounds=np.array([[-1,1],[-1,1]]))
  
    print('test Griewank results:')
    result(func5, bounds=np.array([[-600,600],[-600,600]]))