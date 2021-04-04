from numpy.linalg import eig, inv, norm, cond
import numpy as np
from tabulate import tabulate

myA=[
 [-0.81417, -0.01937, 0.41372],
 [-0.01937, 0.54414, 0.00590],
 [0.41372, 0.00590, -0.81445]
]

def Jacobi():
 E = np.eye(3)
 

print("Вариант 1:")
# 1)
