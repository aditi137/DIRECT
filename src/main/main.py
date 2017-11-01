from direct import Direct, GlobalMin
from helper import *

def result(f, bounds, **kwargs):
    curr_opt, x_at_opt, l_hist = Direct(f,bounds, **kwargs).run()
    print("curr_opt =", curr_opt, ",  x_at_opt =", x_at_opt)
    print()


if __name__ == "__main__":
# Direct parameter list:
# f, bounds, epsilon=1e-4, max_feval=20, max_iter=10, max_rectdiv=100,
# globalmin=GlobalMin(minimize=True, known=False, val=0.), tol=1e-2

    print('test Goldstein-Price results:')
    result(func1, bounds=np.array([[-2,2,],[-2,2,]]), globalmin=GlobalMin(known=True, val=3.))
  
    print('test Rosenbrock results:')
    result(func2, bounds=np.array([[-5,5],[-2,8]]))
   
    print('test Six-hump Camelback results:')
    result(func3, bounds=np.array([[-5,5],[-5,5]]))
      
    print('test Rastrigin results:')
    result(func4, bounds=np.array([[-1,1],[-1,1]]))
  
    print('test Griewank results:')
    result(func5, bounds=np.array([[-600,600],[-600,600]]))