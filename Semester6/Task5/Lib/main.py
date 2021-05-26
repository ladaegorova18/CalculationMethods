import numpy as np
import tabulate
import sympy
import scipy.integrate as integrate
import scipy.misc as misc
from mpmath import chebyt, chop, taylor



def printMatrix(A, B, selected):
    for row in range(len(B)):
        print("(", end='')
        for col in range(len(A[row])):
             print("\t{1:10.2f}{0}".format(" " if (selected is None
or selected != (row, col)) else "*", A[row][col]), end='')
        print("\t) * (\tX{0}) = (\t{1:10.2f})".format(row + 1, B[row]))


class Composable(object):
    def __init__(self, function):
        self.function = function
    def __call__(self, *args, **kwargs):
        return self.function(*args, **kwargs)
    def __mul__(self, other):
        @Composable
        def composed(*args, **kwargs):
            return self.function(other(*args, **kwargs))
        return composed
    def __rmul__(self, other):
        @Composable
        def composed(*args, **kwargs):
            return other(self.function(*args, **kwargs))
        return composed


def infnorm(A):
    return max(A, key=lambda i: abs(i))


def scalar(alpha, beta):
   # d = sympy.Mul(alpha, beta)
    #return integrate.quad(d, -1, 1)
    return 1


def Lu(u):
    return -1 * misc.derivative((u.diff(x)).diff(x)) /(x - 3) + (1 + x/2) * u.diff(x) + np.power(np.e, x / 2) * u


def f(x):
    return 2 - x


print("Проекционные методы решения краевой задачи для обыкновенного дифференциального уравнения второго порядка")
print("Вариант 3:")
print("Метод Галёркина:")
#n = int(input("Введите число координатных функций:"))
n = 3

x = sympy.symbols('x')

polynomials = []
for i in range(n):
    polynomials.append(sympy.jacobi(i, 1, 1, x))

b = (lambda x: polynomials[0])(4)


def polynom(j):
    return polynomials[j]

A = np.eye(n)
B = np.ones((n, 1))

for i in range(n):
    B[i] = scalar(f, polynomials[i])
    for j in range(n):
       # A[i][j] = scalar(Lu(polynomials[j]), polynomials[i])
         A[i][j] = 1

print("Расширенная матрица системы:")
#printMatrix(A, B)
#C = np.linalg.solve(A, B)


print("Число обусловленности матрицы A:")
#mu = infnorm(A) * infnorm(np.linalg.inv(A))


print("Коэффициенты разложения С:")
#print(C)

solution = []

Cond = []
Cond.append(["n", "mu(A)", "y^n(x)", "y*(x) - y^n(x)"])

print("Метод коллокации:")
t = []
for i in range(n):
    a = np.roots(chop(taylor(lambda x: chebyt(i, x), 0, n))[::-1])
    t.append(a)


A1 = np.eye(n)
B1 = np.ones((n, 1))
for i in range(n):
    B[i] = f(t[i])
    for j in range(n):
        A1[i][j] = Lu(polynomials[j])(t[i])

C1 = np.linalg.solve(A1, B1)
print("Расширенная матрица системы:")
printMatrix(A1, B1)


print("Число обусловленности матрицы A:")
mu = infnorm(A1) * infnorm(np.linalg.inv(A1))


print("Коэффициенты разложения С:")
print(C1)
f = 9