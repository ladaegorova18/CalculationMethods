import numpy as np
from tabulate import tabulate
import sympy
import scipy.integrate as integrate
from mpmath import chebyt, chop, taylor

x = sympy.symbols('x')


def f(x):
    return 1 + x / 3

# вывод матрицы на экран
def printMatrix(A, B):
    selected = None
    for row in range(len(B)):
        print("(", end='')
        for col in range(len(A[row])):
            print("\t{1:10.2f}{0}".format(" " if (selected is None
                                                  or selected != (row, col)) else "*", A[row][col]), end='')
        print("\t) * (\tX{0}) = (\t{1:10.2f})".format(row + 1, B[row][0]))

# норма матрицы, где p = infinity
def infnorm(A):
    return max(map(max, A))

# скалярное произведение (Lwj, wi) для метода Галёркина
def scalar(i, j):
    x = sympy.symbols('x')
    wj = sympy.jacobi(j, 1, 1, x) * (1 - np.power(x, 2))
    d = sympy.lambdify(x, sympy.diff(wj, x, 1))
    d2 = sympy.lambdify(x, sympy.diff(wj, x, 2))
    l = integrate.quad(lambda x: (-1 * d2(x) * ((4 - x) / (5 - 2 * x)) + ((1 - x) / 2) * d(x)
                                 + 0.5 * sympy.ln(3 + x) * sympy.jacobi(j, 1, 1, x) * (1 - np.power(x, 2))) * sympy.jacobi(i, 1, 1, x) * (1 - np.power(x, 2)), -1, 1)[0]
    return l

# Lu - для метода коллокации
def Lu(wj):
    x = sympy.symbols('x')
    d = sympy.lambdify(x, sympy.diff(wj, x, 1))
    d2 = sympy.lambdify(x, sympy.diff(wj, x, 2))
    return -1 * d2(x) * ((4 - x) / (5 - 2 * x)) + ((1 - x) / 2) * d(x) + 0.5 * sympy.ln(3 + x) * wj

# формирование решения
def solution(x0, C): 
    result = 0
    for i in range(len(C)):
        result += C[i] * sympy.lambdify(x, w[i])(x0)
    return result

# печать результата для нескольких n
def printResult(Cond):
    headers = ["n", "mu(A)", "y^n(-0.5)", "y^n(0)", "y^n(0.5)", "y*(x) - y^n(x)"]
    print(tabulate(Cond, headers, tablefmt='grid'))


# метод Галёркина
def galerkin(n, w):
    A = np.eye(n)
    B = np.ones((n, 1))

    for i in range(n):
        B[i] = integrate.quad(lambda x: (1 + x / 3) * sympy.jacobi(i, 1, 1, x) * (1 - np.power(x, 2)), -1, 1)[0]
        for j in range(n):
            A[i][j] = scalar(i, j)
    C = np.linalg.solve(A, B)
    mu = infnorm(A) * infnorm(np.linalg.inv(A))
    return A, B, C, mu


# Метод коллокации
def collocation(n, w):
    t = sorted(
        np.squeeze(np.asarray(np.roots(chop(taylor(lambda x: chebyt(n, x), 0, n))[::-1]))))  # Корни многочлена Чебышева

    A = np.eye(n)
    B = np.ones((n, 1))
    for i in range(n):
        B[i] = f(t[i])
        for j in range(n):
            A[i][j] = sympy.lambdify(x, Lu(w[j]))(t[i])
    C = np.linalg.solve(A, B)
    mu = infnorm(A) * infnorm(np.linalg.inv(A))
    return A, B, C, mu


print("Проекционные методы решения краевой задачи для обыкновенного дифференциального уравнения второго порядка")
print("Вариант 8")

Galerkin = []
Collocation = []

for n in range(3, 11):
    w = [] # формируем семейство ортогональных функций — здесь это многочлены Якоби
    for i in range(n):
        w.append(sympy.jacobi(i, 1, 1, x) * (1 - np.power(x, 2)))

    A, B, C, mu = galerkin(n, w)
    A1, B1, C1, mu1 = collocation(n, w)

    Galerkin.append([n, mu, solution(-0.5, C), solution(0, C), solution(0.5, C), np.abs(solution(-0.5, C) - solution(-0.5, C1))])
    Collocation.append([n, mu1, solution(-0.5, C1), solution(0, C1), solution(0.5, C1), np.abs(solution(-0.5, C) - solution(-0.5, C1))])

print("Метод Галёркина:")
printResult(Galerkin)

print("Метод коллокации:")
printResult(Collocation)