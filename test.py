import sys
sys.path.append("build")

import Hilbert
print(dir(Hilbert))

try:
    X = Hilbert.r2t(5, 3, [5, 10, 20])  # Hilbert transpose for 5 bits and 3 dimensions
    print("return type", type(X))
    print("return value", X)    # check [10, 14, 27]
    # check = 7865 or 001111010111001
    print("Hilbert integer = " + \
    str(X[0] >> 4 & 1) + str(X[1] >> 4 & 1) + str(X[2] >> 4 & 1) + str(X[0] >> 3 & 1) + str(X[1] >> 3 & 1) + \
    str(X[2] >> 3 & 1) + str(X[0] >> 2 & 1) + str(X[1] >> 2 & 1) + str(X[2] >> 2 & 1) + str(X[0] >> 1 & 1) + \
    str(X[1] >> 1 & 1) + str(X[2] >> 1 & 1) + str(X[0] >> 0 & 1) + str(X[1] >> 0 & 1) + str(X[2] >> 0 & 1))
    
    print(Hilbert.t2r(5, 3, X)) # check [5, 10, 20]
    print("ok.")
    sys.exit(1)
    
except Exception as e:
    print(e)