from tabulate import tabulate
import math
import numpy as np

class Printer:
    @staticmethod
    # Печать матрицы
    def printMatrix(A):
       for i in range (len(A)):
          for j in range ( len(A[i]) ):
              print ( "{:4f}".format(A[i][j]), end = " " )
          print ()

    @staticmethod
    # Печать таблицы решения на "крупной" сетке 6х6
    def print_grid(U, N, M):
        table = []
        header = ['x/t']
        h = 1 / N
        tau = 0.1 / M
        hx = math.ceil(N / 5)
        ht = math.ceil(M / 5)
        Us = np.zeros((6, 6))
        for i in range(6):
            for k in range(6):
                Us[i][k] = U[i * hx, k * ht]

        for k in range(6):
            header.append(k * ht * tau)

        for i in range(6):
            line = []
            line.append(i * hx * h)
            for k in range(6):
                line.append(Us[i][k])
            table.append(line)
        print(tabulate(table, headers=header, tablefmt="pipe"))

