from direct import Direct, GlobalMin
from helper import *
import datetime

if __name__ == "__main__":
# Direct parameter list:
# f, bounds, epsilon=1e-4, max_feval=200, max_iter=10, max_rectdiv=200,
# globalmin=GlobalMin(minimize=True, known=False, val=0.), tol=1e-2

    with open("direct-run.log", 'a') as file:
        file.write("=================================================\n")
        file.write("RUN MAIN.PY "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M")+"\n")

        print('test Goldstein-Price results:')
        file.write('test Goldstein-Price results:\n')
        Direct(func1, bounds=np.array([[-2,2,],[-2,2,]]), globalmin=GlobalMin(known=True, val=3.)).run(file)

        print('test Rosenbrock results:')
        file.write('test Rosenbrock results:\n')
        Direct(func2, bounds=np.array([[-5,5],[-2,8]])).run(file)
   
        print('test Six-hump Camelback results:\n')
        file.write('test Six-hump Camelback results:\n')
        Direct(func3, bounds=np.array([[-5,5],[-5,5]])).run(file)

        print('test Rastrigin results:')
        file.write('test Rastrigin results:\n')
        Direct(func4, bounds=np.array([[-1,1],[-1,1]]), max_feval=5).run(file)

        print('test Griewank results:\n')
        file.write('test Griewank results:')
        Direct(func5, bounds=np.array([[-600,600],[-600,600]]), max_feval=5).run(file)
    file.close()