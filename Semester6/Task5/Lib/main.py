import numpy as np
from tabulate import tabulate
import sympy
import scipy.integrate as integrate
from mpmath import chebyt, chop, taylor

x = sympy.symbols('x')


def printMatrix(A, B):
    selected = None
    for row in range(len(B)):
        print("(", end='')
        for col in range(len(A[row])):
            print("\t{1:10.2f}{0}".format(" " if (selected is None
                                                  or selected != (row, col)) else "*", A[row][col]), end='')
        print("\t) * (\tX{0}) = (\t{1:10.2f})".format(row + 1, B[row][0]))


def infnorm(A):
    return max(map(max, A))


def scalar(i, j):
    x = sympy.symbols('x')
    wj = sympy.jacobi(j, 1, 1, x)
    d = sympy.lambdify(x, sympy.diff(wj, x, 1))
    d2 = sympy.lambdify(x, sympy.diff(wj, x, 2))
    l = integrate.quad(lambda x: -1 * d2(x) / (x - 3) + (1 + x / 2) * d(x)
                                 + np.power(np.e, x / 2) * sympy.jacobi(j, 1, 1, x) * sympy.jacobi(i, 1, 1, x), -1, 1)[
        0]
    return l


def Lu(wj):
    x = sympy.symbols('x')
    d = sympy.lambdify(x, sympy.diff(wj, x, 1))
    d2 = sympy.lambdify(x, sympy.diff(wj, x, 2))
    return -1 * d2(x) / (x - 3) + (1 + x / 2) * d(x) + np.power(np.e, x / 2) * wj


def f(x):
    return 2 - x


def solution(x0, C):
    result = 0
    for i in range(len(C)):
        result += C[i] * sympy.lambdify(x, w[i])(x0)
    return result


def printresult(Cond):
    headers = ["n", "mu(A)", "y^n(-0.5)", "y^n(0)", "y^n(0.5)", "y*(x) - y^n(x)"]
    print(tabulate(Cond, headers, tablefmt='grid'))


def printSolon(A, B, mu, C, n):
    print("Расширенная матрица:")
    printMatrix(A, B)

    print("Число обусловленности матрицы A:")
    print(mu)

    print("Коэффициенты разложения С:")
    print(C)


def galerkin(n, w):
    A = np.eye(n)
    B = np.ones((n, 1))

    for i in range(n):
        B[i] = integrate.quad(lambda x: (2 - x) * sympy.jacobi(i, 1, 1, x), -1, 1)[0]
        for j in range(n):
            A[i][j] = scalar(i, j)
    C = np.linalg.solve(A, B)
    mu = infnorm(A) * infnorm(np.linalg.inv(A))
    return A, B, C, mu


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
print("Вариант 3:")
v = int(input("Введите число координатных функций или нажмите 0, чтобы оставить значения от 3 до 7:"))

CondGalerkin = []
CondColloc = []
if v == 0:
    for n in range(3, 8):
        w = []
        for i in range(n):
            w.append(sympy.jacobi(i, 1, 1, x))

        A, B, C, mu = galerkin(n, w)
        CondGalerkin.append([n, mu, solution(-0.5, C), solution(0, C), solution(0.5, C), 0])

        A1, B1, C1, mu1 = collocation(n, w)
        CondColloc.append([n, mu1, solution(-0.5, C1), solution(0, C1), solution(0.5, C1), 0])

    print("Метод Галёркина:")
    printresult(CondGalerkin)

    print("Метод коллокации:")
    printresult(CondColloc)
else:
    n = v
    w = []
    for i in range(n):
        w.append(sympy.jacobi(i, 1, 1, x))
    print("Метод Галёркина:")
    A, B, C, mu = galerkin(n, w)
    CondGalerkin.append([n, mu, solution(-0.5, C), solution(0, C), solution(0.5, C), 0])

    printSolon(A, B, mu, C, n)

    print("Метод коллокации:")
    A1, B1, C1, mu1 = collocation(n, w)
    CondColloc.append([n, mu1, solution(-0.5, C1), solution(0, C1), solution(0.5, C1), 0])
    printSolon(A1, B1, mu1, C1, n)
